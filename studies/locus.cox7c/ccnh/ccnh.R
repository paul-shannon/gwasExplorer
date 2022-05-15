# get gtex tissue names:
#    tbl.summary <- gwex$getEQTLSummary()
#    checkTrue(nrow(tbl.summary) > 450)
#    brain.tissue.study.ids <- grep("GTEx_V8.Brain_[CH]", tbl.summary$unique_id, v=TRUE)
#----------------------------------------------------------------------------------------------------
library(gwasExplorer)
library(EndophenotypeExplorer)
library(RUnit)
library(GenomicRanges)
library(rtracklayer)
#----------------------------------------------------------------------------------------------------
eqtl.list <- get(load("../shared/tbl.eqtls.gtex.13brain.tissues.chr5:85,927,378-87,927,378.RData"))
names(eqtl.list)
brain.studies <- sprintf("GTEx_V8.Brain_%s", names(eqtl.list))

FIMO.THRESHOLD <- 5e-5
GENE.EXPRESSION.ABSCOR.MINIMUM <- 0.5
BREAK.ABS.PCT.DELTA.THRESHOLD <- 0.15

targetGene <- "CCNH"
tag.snp <- "rs62374257"
tag.snp.hg19 <- 86223195
tag.snp.hg38 <- 86927378
tag.snp.chrom <- "chr5"

if(!exists("tbl.linkage"))
  tbl.linkage <- get(load("../shared/tbl.linkage.hg19.hg38.RData"))

if(!exists("tbl.fimo"))  # these fimo calls may suffice for CCNH
    tbl.fimo <- get(load(sprintf("../shared/tbl.fimo.%s.RData", "COX7C")))

#----------------------------------------------------------------------------------------------------
# cis-tfbs are within 1M of the tss of the target gene.
# this will work for us if it includes the tag.snp and its ld partners.
specify.region <- function()
{
      # based on ampad eqtls for CCNH
   chrom <- "chr5"
   start <- 86608242
   end <- 88413005

   gr.hg38 <- GRanges(seqnames=chrom, IRanges(start=start, end=end))

   if(!file.exists("hg38ToHg19.over.chain")){
      system("curl -O http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz")
      system("gunzip hg38ToHg19.over.chain.gz")
      }
    chain <- import.chain("hg38ToHg19.over.chain")
    x <- liftOver(gr.hg38, chain)
    gr.hg19 <- unlist(x)
    tbl.roi <- data.frame(chrom=tag.snp.chrom,
                          start.hg19=as.data.frame(gr.hg19)$start[1],
                          end.hg19=as.data.frame(gr.hg19)$end[1],
                          start.hg38=start,
                          end.hg38=end,
                          width=1+end-start,
                          stringsAsFactors=FALSE)
    tbl.roi

} # specify.region
#----------------------------------------------------------------------------------------------------
break.motifs <- function(rsids, motif.names)
{
    message(sprintf("--- break.motifs: %d rsids, %d motif.names", length(rsids), length(motif.names)))

    require(motifbreakR)
    require(BiocParallel)
    require(SNPlocs.Hsapiens.dbSNP155.GRCh38)
    require(BSgenome.Hsapiens.UCSC.hg38)

    motifs.selected <- MotifDb[motif.names]
    print(system.time({
       snps.gr.list <- lapply(rsids,
                           function(rsid) snps.from.rsid(rsid = rsid,
                                                         dbSNP=SNPlocs.Hsapiens.dbSNP155.GRCh38,
                                                         search.genome=BSgenome.Hsapiens.UCSC.hg38))
       snps.gr <- unlist(as(snps.gr.list, "GRangesList"))
       }))

    bpparam <- SerialParam()
    print(system.time(results <- motifbreakR(snpList = snps.gr,
                           filterp = TRUE,
                           pwmList = motifs.selected,
                           show.neutral=FALSE,
                           method = c("ic", "log", "notrans")[1],
                           bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                           BPPARAM = bpparam,
                           verbose=TRUE)))
    tbl.out <- as.data.frame(results, row.names=NULL)
    if(nrow(tbl.out) == 0)
        return(data.frame())
    tbl.out <- subset(tbl.out, effect=="strong")
    tbl.out$pctDelta <- with(tbl.out, pctRef - pctAlt)
    tbl.out

} # break.motifs
#----------------------------------------------------------------------------------------------------
test_break.motifs <- function()
{
    motif.names <- c("Hsapiens-HOCOMOCOv11-core-A-RXRA_HUMAN.H11MO.0.A",
                     "Hsapiens-jaspar2018-HNF4G-MA0484.1")

    rsids <- c("rs6979446")
    rsids <- c("rs6979446", "rs73410456")
    tbl.breaks <- break.motifs(rsids, motif.names)
    checkEquals(dim(tbl.breaks), c(2, 26))
    checkEqualsNumeric(tbl.breaks$pctDelta, c(0.1510949, 0.1809402), tol=1e-4)

    big.test <- FALSE
    if(big.test){
        rsids.20 <- c("rs149352678", "rs74504435", "rs77126132", "rs151057105", "rs75061358",
                      "rs76928645", "rs80013346", "rs11765408", "rs142461330", "rs6979446",
                      "rs6951828", "rs115204793", "rs115849102", "rs142950337", "rs58626582",
                      "rs59010780", "rs60596612", "rs61022057", "rs6975304", "rs77148056")
        motif.names.20 <- c("Hsapiens-HOCOMOCOv11-core-A-TYY1_HUMAN.H11MO.0.A",
                            "Hsapiens-jaspar2018-YY1-MA0095.2",
                            "Hsapiens-HOCOMOCOv11-core-A-ZFP42_HUMAN.H11MO.0.A",
                            "Hsapiens-jaspar2018-GSX2-MA0893.1",
                            "Hsapiens-HOCOMOCOv11-core-A-HXB13_HUMAN.H11MO.0.A",
                            "Hsapiens-jaspar2018-MNX1-MA0707.1",
                            "Hsapiens-jaspar2018-HOXA10-MA0899.1",
                            "Hsapiens-HOCOMOCOv11-core-A-MEF2D_HUMAN.H11MO.0.A",
                            "Hsapiens-jaspar2018-CDX1-MA0878.1",
                            "Hsapiens-jaspar2018-HOXA13-MA0650.1",
                            "Hsapiens-jaspar2018-HOXB13-MA0901.1",
                            "Hsapiens-jaspar2018-HOXD13-MA0909.1",
                            "Hsapiens-jaspar2018-MEF2B-MA0660.1",
                            "Hsapiens-HOCOMOCOv11-core-A-TBP_HUMAN.H11MO.0.A",
                            "Hsapiens-jaspar2018-ZNF263-MA0528.1",
                            "Hsapiens-HOCOMOCOv11-core-A-RXRA_HUMAN.H11MO.0.A",
                            "Hsapiens-HOCOMOCOv11-core-A-KLF15_HUMAN.H11MO.0.A",
                            "Hsapiens-jaspar2018-TCF7L1-MA1421.1",
                            "Hsapiens-jaspar2018-RXRG-MA0856.1",
                            "Hsapiens-jaspar2018-HNF4G-MA0484.1")

        t0 <- system.time(
            tbl.breaks.2020 <- break.motifs(rsids.20, motif.names.20)
            )
            #    user  system elapsed
            #  60.409  19.815  80.362
        checkTrue(nrow(tbl.breaks) > 100)
        } # big.test

} # test_break.motifs
#----------------------------------------------------------------------------------------------------
# ampad eqtls apparently have no beta, just pval.
# every trena tf is scored as the sum of the break score, for every variant judged to break
# its motif, * -log10 of the eqtls pvalue
score.breakage <- function(tbl.trena, tbl.ampad.eqtls, tbl.breaks)
{
   tfs.oi <- tbl.trena$gene

   if(!"sig.no.beta" %in% colnames(tbl.ampad.eqtls)){
       sig.no.beta <- with(tbl.ampad.eqtls, -log10(pvalue)) # * abs(beta))
       tbl.ampad.eqtls$sig.no.beta <- sig.no.beta
       }

   breakage.scores <- list()

   for(tf in tfs.oi){
      #tbl.tf.breaks.sub <- subset(tbl.breaks, geneSymbol == tf & SNP_id %in% tbl.eqtl$rsid)
      tbl.tf.breaks.sub <- subset(tbl.breaks, geneSymbol == tf & SNP_id %in% tbl.ampad.eqtls$rsid)
      dups <- which(duplicated(tbl.tf.breaks.sub[, c("SNP_id", "geneSymbol")]))
      if(length(dups) > 0)
          tbl.tf.breaks.sub <- tbl.tf.breaks.sub[-dups,]
      tf.breakage.score <- 0
      for(breaking.snp in tbl.tf.breaks.sub$SNP_id){
          pctDelta <- abs(subset(tbl.tf.breaks.sub, SNP_id==breaking.snp)$pctDelta)
              # multiply this by the tf's rfNorm and the eQTLs sig.no.beta
          sig.no.beta <- max(subset(tbl.ampad.eqtls, rsid==breaking.snp)$sig.no.beta)
          tf.rfNorm <- subset(tbl.trena, gene==tf)$rfNorm
          new.increment <- (pctDelta + sig.no.beta) * tf.rfNorm
             # printf("%20s: %f", breaking.snp, new.increment)
          tf.breakage.score <- tf.breakage.score + new.increment
          } # for breaking.snp
      # printf("--- %s: %5.2f", tf, tf.breakage.score)
      breakage.scores[[tf]] <- tf.breakage.score
      } # for tf

    tbl.trena$breakage.score <- round(as.numeric(breakage.scores), digits=2)
    tbl.trena

} # score.breakage
#----------------------------------------------------------------------------------------------------
write.failure.file <- function(output.filename, short.tissue.name, brain.tissue, tbl.breaks,
                               tbl.tms, tbl.tms.wt, eqtls, status)
{
    new.filename <- sub("gwex", "gwex-FAILED", output.filename)
    tbl.trena <- data.frame()
    save(short.tissue.name, brain.tissue, eqtls, tbl.breaks, tbl.tms, tbl.tms.wt,
         tbl.trena, status, file=new.filename)

} # write.failure.file
#----------------------------------------------------------------------------------------------------
main <- function()
{
    for(brain.tissue in brain.studies[c(5,7)]){
        eqtls <- list()
        status <- "success"
        tbl.breaks <- data.frame()
        tbl.scored <- data.frame()
        tbl.tms <- data.frame()
        tbl.tms.wt <- data.frame()
        output.filename <- sprintf ('gwex.%s.%s.%s.RData', targetGene, brain.tissue,
                                    format (Sys.time(), "%a.%b.%d.%Y-%H:%M:%S"))
        printf("--- starting on %s", output.filename )
        gwex <- gwasExplorer$new(targetGene=targetGene, locusName="COX7C", tagSnp=tag.snp,
                                 shoulder=0, tissueName=brain.tissue,
                                 tbl.prespecifiedRegion=specify.region(),
                                 tbl.linkage=tbl.linkage)

           # previously obtained eqtls
        short.tissue.name <- sub("GTEx_V8.Brain_", "", brain.tissue)
        stopifnot(short.tissue.name %in% names(eqtl.list))
        tbl.eqtl <- eqtl.list[[short.tissue.name]]
        if(nrow(tbl.eqtl) == 0){
            message(sprintf("no gtex eqtls in %s for %s", short.tissue.name, targetGene))
            write.failure.file(output.filename, short.tissue.name, brain.tissue, tbl.breaks,
                               tbl.tms, tbl.tms.wt, eqtls, status="no gtex eqtls of any pval in tissue")
            next;
            }
        deleters <- which(is.na(tbl.eqtl$hg38))
        if(length(deleters) > 0)
            tbl.eqtl <- tbl.eqtl[-deleters,]

        eqtl.pval.max <- 0.05
        tbl.eqtl <- subset(tbl.eqtl, pvalue <= eqtl.pval.max & gene==targetGene)
        if(nrow(tbl.eqtl) == 0){
            message(sprintf("no gtex eqtls in %s for %s with pval <=", short.tissue.name, targetGene, eqtl.pval.max))
            write.failure.file(output.filename, short.tissue.name, brain.tissue, tbl.breaks,
                               tbl.tms, tbl.tms.wt, eqtls, status="no gtex eqtls below threshold")
            next;
            }
        eqtls <- gwex$getEqtlsForGene(eqtl.catalog.studyIDs=brain.tissue,
                                      pval.threshold=eqtl.pval.max,
                                      tbl.eqtls.gtex=tbl.eqtl)

        unique.variants <- unique(c(eqtls$ampad$rsid, eqtls$gtex$rsid))
        length(unique.variants)   # 65
        if(length(unique.variants)==0){
            message(sprintf("no unique variants from gtex or ampad in %s for %s ", short.tissue.name, targetGene))
            write.failure.file(output.filename, short.tissue.name, brain.tissue, tbl.breaks,
                               tbl.tms, tbl.tms.wt, eqtls, status="no variants found")
            next;
            }

           #---------------------------------------------------------
           # now create the tms table. 2 eqtl columns should appear
           #---------------------------------------------------------

        data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
        filename <- "boca-hg38-consensus-ATAC.RData"
        tbl.boca <- get(load(file.path(data.dir, filename)))

        #browser()
        tbl.tms <- gwex$trenaMultiScore(brain.tissue, tbl.fimo=tbl.fimo, tbl.oc=tbl.boca)
        dim(tbl.tms)
        print(colnames(tbl.tms))


           #---------------------------------------------------------
           # build wild-type model
           #---------------------------------------------------------
        build.wt.model <- function(top.tfs, mtx.rna, tbl.tms.wt){
            top.tfs <- sort(unique(tbl.tms.wt$tf))
            length(top.tfs)  # 44
            tbl.trena <- gwex$runTrena(top.tfs, mtx.rna, tbl.tms.wt)
            tbl.trena$rfNorm <- round(tbl.trena$rfNorm, digits=3)
            tbl.trena
            } # build.wt.model


             #--------------------------------------------------------------
             # identify the ampad variants in these regions, and the motifs
             #--------------------------------------------------------------

        tbl.tms.wt <- subset(tbl.tms,
                             (gh > 1 | gtex.eqtl | ampad.eqtl | chip) &
                             abs(cor.all) > GENE.EXPRESSION.ABSCOR.MINIMUM &
                             fimo_pvalue < FIMO.THRESHOLD)
        dim(tbl.tms.wt)

        if(nrow(tbl.tms.wt) == 0){
            message(sprintf("empty tbl.tms.wt in %s for %s ", short.tissue.name, targetGene))
            write.failure.file(output.filename, short.tissue.name, brain.tissue, tbl.breaks,
                               tbl.tms, tbl.tms.wt, eqtls, status="empty tbl.tms.wt")
            next;
            }

        gr.tfbs <- GRanges(tbl.tms.wt[, c("chrom", "start", "end")])
        tbl.ampad.eqtls <- eqtls$ampad[, c("chrom", "hg38", "hg38", "rsid", "pvalue", "genesymbol")]
        colnames(tbl.ampad.eqtls) <- c("chrom", "start", "end", "rsid", "pvalue", "genesymbol")
        gr.eqtls <- GRanges(tbl.ampad.eqtls)
        tbl.ov <- as.data.frame(findOverlaps(gr.eqtls, gr.tfbs))
        message(sprintf("--- overlaps between filtered tbl.tms.wt and ampad eqtls: %d", nrow(tbl.ov)))
        if(nrow(tbl.ov) > 0){
           tbl.ampad.eqtls.tfbs <- tbl.ampad.eqtls[unique(sort(tbl.ov$queryHits)),]
           rsids.in.tfbs <- unique(tbl.ampad.eqtls[tbl.ov$queryHits, "rsid"])
           motif.names <- unique(tbl.tms.wt[tbl.ov[,2], "motif_id"])
           tbl.breaks <- break.motifs(rsids.in.tfbs, motif.names)
           if(nrow(tbl.breaks) == 0){
               write.failure.file(output.filename, short.tissue.name, brain.tissue, tbl.breaks,
                                  tbl.tms, tbl.tms.wt, eqtls, status="no breaks in tfbs")
               next;
               }
           new.order <- order(abs(tbl.breaks$pctDelta), decreasing=TRUE)
           tbl.breaks <- tbl.breaks[new.order,]
           rownames(tbl.breaks) <- NULL
           colnames(tbl.breaks)[1] <- "chrom"
           tbl.breaks$chrom <- as.character(tbl.breaks$chrom)
           broi <- c("chrom", "start", "end", "SNP_id", "geneSymbol", "providerId", "seqMatch", "pctRef", "pctAlt", "pctDelta")
           tbl.breaks.easy <- tbl.breaks[, broi]
           tbl.breaks.easy$start <- tbl.breaks.easy$start - 1
           dups <- which(duplicated(tbl.breaks.easy[, c("chrom", "start", "end", "SNP_id", "geneSymbol", "providerId")]))
           if(length(dups) > 0)
               tbl.breaks.easy <- tbl.breaks.easy[-dups,]
           }

        top.ranked <- table(tbl.tms.wt$tf)
        print(top.ranked)
        top.tfs <- names(top.ranked)
        length(top.tfs)  # 6

        mtx.rna <- gwex$getExpressionMatrix(brain.tissue)
        message(sprintf("--- building wt model for %s in %s", targetGene, short.tissue.name))
        tbl.trena <- build.wt.model(top.tfs, mtx.rna, tbl.tms.wt)
        dim(tbl.trena)
        eqtl.pval.max <- 1e-4
        #browser()
        if(nrow(tbl.breaks) > 0){
           tbl.scored <- score.breakage(tbl.trena, tbl.ampad.eqtls, tbl.breaks)
           message(sprintf("------- results for %s", brain.tissue))
           message(sprintf("------- tbl.scored, sum top 10: %f", sum(head(tbl.scored$breakage.score, n=10))))
           message(sprintf("--- saving to %s", output.filename))
           }
        save(short.tissue.name, brain.tissue, eqtls, tbl.tms, tbl.trena, tbl.breaks, tbl.scored, file=output.filename)
        } # for brain.tissue

} # main
#----------------------------------------------------------------------------------------------------
viz <- function()
{
   if(!exists("igv")){
      igv <- start.igv(targetGene, "hg38")
      shoulder <- 5e5
      showGenomicRegion(igv, sprintf("%s:%d-%d", tag.snp.chrom, tag.snp.hg38-shoulder, tag.snp.hg38+shoulder))
      tbl.track <- data.frame(chrom=tag.snp.chrom, start=tag.snp.hg38-1, end=tag.snp.hg38, stringsAsFactors=FALSE)
      track <- DataFrameAnnotationTrack(tag.snp, tbl.track, trackHeight=25, color="darkred")
      displayTrack(igv, track)
      threshold <- 0.2
      tbl.track <- subset(tbl.linkage, r2 >= threshold)[, c("chrom", "hg38", "hg38", "r2", "rsid")]
      colnames(tbl.track) <- c("chrom", "start", "end", "score", "rsid")
      tbl.track$chrom <- tag.snp.chrom
      tbl.track$start <- tbl.track$start - 1
      track.title <- sprintf("haploreg eur, r^2 > %3.2f", threshold)
      track <- DataFrameQuantitativeTrack(track.title, tbl.track, autoscale=FALSE, min=0, max=1, color="black")
      displayTrack(igv, track)
      }

   etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD", initialize.snpLocs=TRUE)
   tbl.eqtls <- etx$get.ampad.EQTLsForGene()
   #tbl.eqtls <- subset(tbl.eqtls, study=="ampad-rosmap")
   dim(tbl.eqtls)
   head(tbl.eqtls, n=20)
   tail(tbl.eqtls, n=20)
   tbl.track <- subset(tbl.eqtls, pvalue < 1e-2 & study != "GTEx")[, c("chrom", "hg38", "hg38", "rsid", "pvalue")]
   dim(tbl.track)
   colnames(tbl.track)[2:3] <- c("start", "end")
   tbl.track$start <- tbl.track$start - 1
   colnames(tbl.track) <- NULL
   track <- GWASTrack("ampad", tbl.track, trackHeight=100) # chrom.col=0, pos.col=0, pval.col=0)
   displayTrack(igv, track)

       #--------------------------------------
       # eqtls from gtex, one track per gene
       #--------------------------------------

   tbl.eqtl.all <- do.call(rbind, eqtl.list)
   dim(tbl.eqtl.all)
   genes <- sort(unique(tbl.eqtl.all$gene))
   length(genes) # 9
   for(goi in genes) {
      tbl.track <- subset(tbl.eqtl.all, gene==goi)[, c("chrom", "hg38", "hg38", "rsid", "pvalue")]
      rownames(tbl.track) <- NULL
      colnames(tbl.track)[2:3] <- c("start", "end")
      tbl.track$start <- tbl.track$start - 1
      tbl.track$chrom <- paste0("chr", tbl.track$chrom)
      deleters <- which(is.na(tbl.track$start))
      if(length(deleters) > 0)
          tbl.track <- tbl.track[-deleters,]
      colnames(tbl.track) <- NULL
      track <- GWASTrack(goi, tbl.track, trackHeight=100) # chrom.col=0, pos.col=0, pval.col=0)
      displayTrack(igv, track)
      } # for goi

   removeTracksByName(igv,
                      getTrackNames(igv)[5:9]
                      )


      # now broken motifs of tfs in model
   tfs.oi <- tbl.trena$gene[1:20]
   for(TF in tfs.oi){
      tbl.breaks.sub <- subset(tbl.breaks.easy, geneSymbol==TF)[, c("chrom", "start", "end", "pctDelta")]
      tbl.breaks.sub$pctDelta <- -1.0 * tbl.breaks.sub$pctDelta
      if(nrow(tbl.breaks.sub) > 0){
         track <- DataFrameQuantitativeTrack(TF, tbl.breaks.sub, color="random", autoscale=FALSE, min=-0.2, max=0.2)
         displayTrack(igv, track)
         }
      }

   dim(tbl.tms.sub)
   tfs <- sort(unique(tbl.tms.sub$tf))
   length(tfs)


   for(TF in tfs){
      tbl.track <- subset(tbl.tms.sub, tf==TF)[, c("chrom", "start", "end", "cor.all")]
      track <- DataFrameQuantitativeTrack(TF, tbl.track, color="random", autoscale=FALSE, min=-1, max=1)
      displayTrack(igv, track)
      } # for TF

} # viz
#----------------------------------------------------------------------------------------------------
source("../shared/doFile.R")
#----------------------------------------------------------------------------------------------------
# motifBreak.score <- with(tbl.tf, abs(pctDelta) * 100)
#       eqtl.score <- with(tbl.tf, -log10(gtex.eqtl.pval)* abs(gtex.eqtl.beta) * 100)
#     trena.score <- with(tbl.tf, (abs(betaLasso) * 100) + (rfNorm * 10))
#      tfbs.score <- 1/tfbs.count
# score <- trena.score * eqtl.score * motifBreak.score * tfbs.score
review <- function()
{
  files <- list.files(".", pattern=sprintf("gwex.%s.*RData", targetGene))
  files

  tbls.all <- list()
  for(file in files){
      tbl <- do.file(file)
      tbls.all[[file]] <- tbl
      } # for file


  tbl.all <- do.call(rbind, tbls.all)
  dim(tbl.all)
  new.order <- order(tbl.all$total.score, decreasing=TRUE)
  tbl.all <- tbl.all[new.order,]
  rownames(tbl.all) <- NULL

  coi <- c("tissue", "tf", "targetGene", "trena", "eqtl", "motifBreak", "tfbs.count", "rsids", "total.score")
  tbl.all <- tbl.all[, coi]
  tbl.all

     # now inspect
   TF <- "POU3F2"
   tbl.tms.TF <- subset(tbl.tms, tf==TF & ampad.eqtl)
   dim(tbl.tms.TF)
   track <- DataFrameAnnotationTrack(TF, tbl.tms.TF, color="random", trackHeight=25)
   displayTrack(igv, track)

   tbl.breaks.TF <- subset(tbl.breaks, geneSymbol==TF & abs(pctDelta) > 0.05) [, breaks.coi]
   track <- DataFrameQuantitativeTrack(sprintf("%s.rsids", TF), tbl.breaks.TF[, c(1,2,3,7)],
                                       color="brown", autoscale=FALSE, min=-0.2, max=0.2)
   displayTrack(igv, track)

} # review
#----------------------------------------------------------------------------------------------------
