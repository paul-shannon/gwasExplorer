#----------------------------------------------------------------------------------------------------
library(gwasExplorer)
library(ghdb)
library(EndophenotypeExplorer)
library(ADvariantExplorer)
library(TrenaProjectAD)
library(motifbreakR)
library(MotifDb)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)

source("~/github/TrenaMultiScore/tools/runner/v2/tmsCore.R")
#----------------------------------------------------------------------------------------------------
targetGene <- "TFAM"
#----------------------------------------------------------------------------------------------------
if(!exists("etx")){  # let etx be the proxy for all initializations
    etx <- EndophenotypeExplorer$new(targetGene, "hg38", vcf.project="ADNI", initialize.snpLocs=FALSE)
    trenaProject <- TrenaProjectAD()
    tbl.fimo <- get(load("data/tbl.fimo.TFAM.RData"))
       # 16k straddling TSS: chr10:58,377,275-58,393,364
    tbl.fimo.small <- subset(tbl.fimo, chrom=="chr10" & start >= 58377275 & end <= 58393364)
    dim(tbl.fimo.small)   # 27538

    data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
    filename <- "boca-hg38-consensus-ATAC.RData"
    tbl.boca <- get(load(file.path(data.dir, filename)))

       #tbl.snpLocs <- get(load("data/snpLocs-TFAM.RData"))
    print(load("data/ampad.eqtls.TFAM.RData"))
    dim(tbl.ampad.eqtls)
    table(tbl.ampad.eqtls$study)   #   ampad-mayo ampad-rosmap         GTEx
                               #         8740         4429           28
    eqtl.pval.max <- 1e-4
    tbl.eqtls.rosmap <- subset(tbl.ampad.eqtls, study=="ampad-rosmap" & pvalue <= eqtl.pval.max)
    dim(tbl.eqtls.rosmap)   # 151 10

    genomic.shoulder <- 1000
    gene.expression.absCor.minimum <- 0.2
    tbl.eqtl <- get(load("data/tbl.eqtls.6.tissues-chr10:57385407-59385407.2022-06-08.08:03:52"))
    tissues <- unique(tbl.eqtl$id)
    tbl.eqtls.gtex <- subset(tbl.eqtl, pvalue <= eqtl.pval.max)
    table(tbl.eqtls.gtex$id)
   # GTEx_V8.Brain_Cerebellar_Hemisphere            GTEx_V8.Brain_Cerebellum
   #                                  72                                 362
   #                GTEx_V8.Brain_Cortex    GTEx_V8.Brain_Frontal_Cortex_BA9
   #                                 634                                  58
   #           GTEx_V8.Brain_Hippocampus          GTEx_V8.Brain_Hypothalamus
   #                                  21                                  27
   } # initialize widely used variables
#----------------------------------------------------------------------------------------------------
do.tms <- function(brain.tissue, tbl.fimo)
{
    tms <- TMS$new(trenaProject, targetGene, tbl.fimo, tbl.boca, quiet=FALSE)
    tms$scoreFimoTFBS()   # chip, conservation, genehancer, genic annotations, distance to tss

    names(etx$get.rna.matrix.codes())
    mtx.rna <- etx$get.rna.matrix(brain.tissue)
    dim(mtx.rna)
    tms$add.tf.mrna.correlations(mtx.rna, featureName="cor.all")

    tbl.tms <- tms$getTfTable()
    gr.tms <- GRanges(tbl.tms)

    #------------------------------------------------------
    # add the ampd eqtls to each fimo region they intersect
    #------------------------------------------------------
    gr.ampad <- with(tbl.ampad.eqtls, GRanges(seqnames=chrom[1], IRanges(hg38)))
    tbl.ov <- as.data.frame(findOverlaps(gr.ampad, gr.tms))
    dim(tbl.ov)
    ampad.eqtl <- rep(FALSE, nrow(tbl.tms))
    if(nrow(tbl.ov) > 0)
        ampad.eqtl[unique(tbl.ov[,2])] <- TRUE
    tbl.tms$ampad.eqtl <- ampad.eqtl

    #-----------------------------------------------------------------------
    # add the tissue-specific gtex eqtls to each fimo region they intersect
    #--------------------------------------------------------------------
    tbl.gtex.tissue.eqtls <- subset(tbl.eqtls.gtex, id==brain.tissue)
    dim(tbl.gtex.tissue.eqtls)
    gr.gtex <- with(tbl.gtex.tissue.eqtls, GRanges(seqnames=sprintf("chr%s", chrom[1]), IRanges(hg38)))
    tbl.ov <- as.data.frame(findOverlaps(gr.gtex, gr.tms))
    dim(tbl.ov)
    gtex.eqtl <- rep(FALSE, nrow(tbl.tms))
    if(nrow(tbl.ov) > 0)
        gtex.eqtl[unique(tbl.ov[,2])] <- TRUE
    tbl.tms$gtex.eqtl <- gtex.eqtl


    tbl.tms.filtered <- subset(tbl.tms, ampad.eqtl & gtex.eqtl & abs(cor.all) > 0.4)
    tf.candidates <- unique(tbl.tms.filtered)$tf
    if(length(tf.candidates) == 0)
        return(list(trena=data.frame(), tms=data.frame()))

    solvers=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest")

    solver <- EnsembleSolver(mtx.rna,
                             targetGene=targetGene,
                             candidateRegulators=tf.candidates,
                             geneCutoff=1.0,
                             solverNames=solvers)
    tbl.trena <- trena::run(solver)
    new.order <- order(tbl.trena$rfScore, decreasing=TRUE)
    tbl.trena <- tbl.trena[new.order,]
    tbl.trena$rfNorm <- tbl.trena$rfScore / max(tbl.trena$rfScore)
    tbl.trena <- cbind(tbl.trena[1], as.data.frame(lapply(tbl.trena[-1], function(x) round(x, digits=2))))
    tbl.trena$tfbs <- unlist(lapply(tbl.trena$gene, function(tf) length(grep(tf, tbl.tms.filtered$tf))))

    return(list(trena=tbl.trena, tms=tbl.tms.filtered))

} # do.tms
#----------------------------------------------------------------------------------------------------
run <- function()
{
   brain.tissues <- unique(tbl.eqtls.gtex$id)
   x.small <- lapply(brain.tissues, function(tissue) do.tms(tissue, tbl.fimo.small))
   names(x.small) <- brain.tissues
   lapply(x.small, function(e) dim(e$trena))
   lapply(x.small, function(e) dim(e$tms))

   x <- lapply(brain.tissues, function(tissue) do.tms(tissue, tbl.fimo))
   names(x) <- brain.tissues
   lapply(x, function(e) dim(e$trena))
   lapply(x, function(e) dim(e$tms))
   save(x.small, x, file="tfam.trena.results.RData")

} # run
#----------------------------------------------------------------------------------------------------
viz <- function()
{
   if(!exists("igv")){
      igv <- start.igv(targetGene, "hg38")
      } #

   coi <- c("chrom", "hg38", "hg38", "pvalue")
   tbl.rosmap <-
       subset(tbl.ampad.eqtls, study=="ampad-rosmap" & pvalue < 0.05)[, coi]
   dim(tbl.rosmap)
   colnames(tbl.rosmap) <- c("chrom", "start", "end", "score")
   tbl.rosmap$start <- tbl.rosmap$start - 1
   tbl.rosmap$score <- -log10(tbl.rosmap$score)
   track <- DataFrameQuantitativeTrack("rosmap", tbl.rosmap, autoscale=FALSE,
                                       min=0, max=22, color="darkblue")
   displayTrack(igv, track)

   tbl.mayo <-
       subset(tbl.ampad.eqtls, study=="ampad-mayo" & pvalue < 0.05)[, coi]
   dim(tbl.mayo)
   colnames(tbl.mayo) <- c("chrom", "start", "end", "score")
   tbl.mayo$start <- tbl.mayo$start - 1
   tbl.mayo$score <- -log10(tbl.mayo$score)
   track <- DataFrameQuantitativeTrack("mayo", tbl.mayo, autoscale=FALSE,
                                       min=0, max=22,color="darkgreen")
   displayTrack(igv, track)

   ghdb <- GeneHancerDB()
   tbl.gh <- retrieveEnhancersFromDatabase(ghdb, targetGene, tissues="all")
   tbl.gh$score <- asinh(tbl.gh$combinedscore)
   track <- DataFrameQuantitativeTrack("GH", tbl.gh[, c("chrom", "start", "end", "score")],
                                       autoscale=TRUE, color="brown")
   displayTrack(igv, track)

   for(TISSUE in tissues){
       coi <- c("chrom", "hg38", "hg38", "pvalue", "beta")
       tbl.sub <- subset(tbl.eqtl, tissue==TISSUE & pvalue < 0.05 & !is.na(chrom))[, coi]
       #tbl.sub$score <- -log10(tbl.sub$pvalue) * tbl.sub$beta
       tbl.sub$score <- -log10(tbl.sub$pvalue)
       colnames(tbl.sub) <- c("chrom", "start", "end", "pvalue", "beta", "score")
       tbl.sub$start <- tbl.sub$start - 1
       tbl.track <- tbl.sub[, c("chrom", "start", "end", "score")]
       printf("%s: %d eqtls", TISSUE, nrow(tbl.sub))
       #title <- sprintf("%s score", TISSUE)
       title <- sprintf("%s pval", TISSUE)
       track <- DataFrameQuantitativeTrack(title, tbl.track, color="random", autoscale=FALSE,
                                           min=0, max=20)
                                           #min=-5, max=5)
       displayTrack(igv, track)
       } # for TISSUE


} # viz
#----------------------------------------------------------------------------------------------------
break.models <- function()
{
    print(load("tfam.trena.results.RData"))
    names(x)
    names(x.small)
    print(load("data/ampad.eqtls.TFAM.RData"))
    table(tbl.ampad.eqtls$study)   #   ampad-mayo ampad-rosmap         GTEx
                                   #         8740         4429           28
    eqtl.pval.max <- 1e-4
    tbl.eqtls <- subset(tbl.ampad.eqtls, study=="ampad-rosmap" & pvalue <= eqtl.pval.max)
    dim(tbl.eqtls)   # 151 10

    tissues <- names(x)
    for(tissue in tissues){
        tbl.trena <- x[[tissue]]$trena
        tbl.trena.tmp <- subset(tbl.trena, rfNorm > 0.2 | abs(betaLasso) > 0.05)
        rownames(tbl.trena.tmp) <- NULL
        tbl.tms   <- x[[tissue]]$tms
        tbl.tfbs <- subset(tbl.tms, tf %in% tbl.trena.tmp$gene)
        dim(tbl.tfbs)
        gr.tfbs <- GRanges(tbl.tfbs[, c("chrom", "start", "end")])
        tbl.eqtls.tmp <- tbl.eqtls[, c("chrom", "hg38", "hg38")]
        colnames(tbl.eqtls.tmp) <- c("chrom", "start", "end")
        gr.eqtls <- GRanges(tbl.eqtls.tmp)
        tbl.ov <- as.data.frame(findOverlaps(gr.eqtls, gr.tfbs))
        tbl.eqtls.ov <- tbl.eqtls[tbl.ov$queryHits,]
        tbl.tfbs.ov <- tbl.tfbs[tbl.ov$subjectHits,]

        tfs.oi <- unique(tbl.tfbs.ov$tf)
        tfs.oi <- tfs.oi[match(tbl.trena.tmp$gene,tfs.oi)]
        deleters <- which(is.na(tfs.oi))
        if(length(deleters) > 0)
            tfs.oi <- tfs.oi[-deleters]
        rsids.oi <- unique(tbl.eqtls.ov$rsid)
        mdb.human <- query(MotifDb, "sapiens", c("jaspar2022", "hocomoco-core-A"))
        motifs.selected <- query(mdb.human, andStrings="sapiens", orStrings=tfs.oi)

        x <- system.time(snps.gr <- snps.from.rsid(rsid = rsids.oi,
                                           dbSNP=SNPlocs.Hsapiens.dbSNP155.GRCh38,
                                           search.genome=BSgenome.Hsapiens.UCSC.hg38))
        printf("snps.gr for %d snps obtained, elapsed time: %f", length(rsids.oi), x[["elapsed"]])

        bpparam <- MulticoreParam(workers=3)
        results <- motifbreakR(snpList = snps.gr,
                               filterp = TRUE,
                               pwmList = motifs.selected,
                               show.neutral=FALSE,
                               method = c("ic", "log", "notrans")[1],
                               bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                               BPPARAM = bpparam,
                               verbose=TRUE)
        tbl.breaks <- as.data.frame(results, row.names=NULL)
        tbl.breaks <- subset(tbl.breaks, effect=="strong")
        tbl.breaks$start <- tbl.breaks$start - 1
        tbl.breaks$pctDelta <- with(tbl.breaks, pctAlt-pctRef)
        new.order <- order(abs(tbl.breaks$pctDelta), decreasing=TRUE)
        tbl.breaks <- tbl.breaks[new.order,]
        colnames(tbl.breaks)[1] <- "chrom"
        tbl.breaks$chrom <- as.character(tbl.breaks$chrom)
        motifbreakR.coi <- c("chrom", "start", "end", "SNP_id", "geneSymbol", "providerName", "pctRef", "pctAlt", "pctDelta")
        tbl.broken <- subset(tbl.breaks, pctDelta < -0.1 & pctRef > 0.6)[, motifbreakR.coi]
        tbl.broken <- subset(tbl.broken, geneSymbol %in% tfs.oi)
        tbl.broken$score <- with(tbl.broken, pctRef - pctDelta)
        rownames(tbl.broken) <- NULL
        colnames(tbl.broken)[4] <- "rsid"
        timestamp <- sub(" ", "", tolower(format(Sys.time(), "%Y.%b.%e-%H:%M")))
        filename <- sprintf("tbl.breaks.%s.%s", targetGene, timestamp)
        save(tbl.breaks, results, file=filename)
        if(FALSE){
           if(!exists("igv")) igv <- start.igv(targetGene, "hg38")
           shoulder <- 1000
           roi <- with(tbl.tfbs.ov, sprintf("%s:%d-%d", chrom[1], min(start)+shoulder, max(end)+shoulder))
           showGenomicRegion(igv,roi)
           for(tf.oi in tfs.oi){
               tbl.track <- subset(tbl.tfbs.ov, tf==tf.oi)[, 1:4]
               track <- DataFrameAnnotationTrack(tf.oi, tbl.track, color="random", trackHeight=25)
               displayTrack(igv, track)
               }
           tbl.track <- tbl.eqtls.ov[, c("chrom", "hg38", "hg38", "pvalue", "rsid")]
           colnames(tbl.track) <- c("chrom", "start", "end", "score", "rsid")
           tbl.track$start <- tbl.track$start - 1
           tbl.track$score <- -log10(tbl.track$score)
           track <- DataFrameQuantitativeTrack("rosmap eQTL", tbl.track, autoscale=TRUE, color="red")
           displayTrack(igv, track)
           tbl.track <- tbl.broken[, c("chrom", "start", "end", "score")]
           track <- DataFrameQuantitativeTrack("broken motifs", tbl.track, autoscale=TRUE, col="brown")
           displayTrack(igv, track)
           }
      } # for tissue

} # break.models
#----------------------------------------------------------------------------------------------------
