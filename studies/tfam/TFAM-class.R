library(R6)
library(EndophenotypeExplorer)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifbreakR)
library(BiocParallel)
source("~/github/TrenaMultiScore/tools/runner/v2/tmsCore.R")
#----------------------------------------------------------------------------------------------------
#' @title tfam
#' @description A template for building documented, tested R6 classes
#' @name tfam

#' @field id identifier for a class object
#'
#' @examples
#'   rt <- R6Template$new(id="abc")
#' @export

TFAM = R6Class("TFAM",

    #--------------------------------------------------------------------------------
    private = list(targetGene=NULL,
                   trenaProject=NULL,
                   study.region=NULL,
                   tbl.fimo=NULL,
                   tbl.oc=NULL,
                   tbl.ampad.eqtls=NULL,
                   tbl.gtex.eqtls=NULL,
                   tbl.trena=NULL,
                   gtex.eqtl.tissues=NULL,
                   current.tissue=NULL,
                   tbl.rosmap.eqtls=NULL,
                   tms=NULL,
                   tbl.tms=NULL,
                   tbl.tmsFiltered=NULL,   # has just the TFs and regions used by trena
                   etx=NULL,
                   mtx.rna=NULL,
                   tbl.breaks=NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(
         #' @description
         #' Creates a new instance of this [R6][R6::R6Class] class.
         #' @param id character, an indentifier for this object
         #' @return a new instance o tfam
        initialize = function(targetGene, trenaProject, tbl.fimo, tbl.gtex.eqtls, tbl.ampad.eqtls,
                              tbl.oc){
            private$targetGene <- targetGene
            private$trenaProject <- trenaProject
            private$tbl.fimo <- tbl.fimo
            private$tbl.gtex.eqtls <- tbl.gtex.eqtls
            private$tbl.ampad.eqtls <- tbl.ampad.eqtls
            private$gtex.eqtl.tissues <- sort(unique(tbl.gtex.eqtls$id))
            private$current.tissue <- private$gtex.eqtl.tissues[1]
            private$tbl.oc <- tbl.oc
            private$study.region <- self$getFimoGenomicRegion()
            private$etx <- EndophenotypeExplorer$new(targetGene, "hg38", vcf.project="ADNI",
                                                     initialize.snpLocs=FALSE)
            },
        #------------------------------------------------------------
        setStudyRegion = function(chrom, start, end){
            width.kb <- round(((1 + end - start)/1000), digits=2)
            private$study.region <- list(chrom=chrom, start=start, end=end, width.kb=width.kb)
            },
        #------------------------------------------------------------
        getStudyRegion = function(){
            private$study.region
            },
        #------------------------------------------------------------
        getFimoGenomicRegion = function(){
            chrom <- private$tbl.fimo$chrom[1]
            start <- min(private$tbl.fimo$start)
            end   <- max(private$tbl.fimo$end)
            width.kb <- round((1 + end - start)/1000, digits=2)
            return(list(chrom=chrom, start=start, end=end, width.kb=width.kb))
            },
        #------------------------------------------------------------
        get.gtex.eqtls = function(){
            invisible(tbl.gtex.eqtls)
            },
        #------------------------------------------------------------
        get.ampad.eqtls = function(){
            invisible(tbl.ampad.eqtls)
            },
        #------------------------------------------------------------
        getGTEx.eqtl.genomicRegion = function(){
            start <- min(tbl.gtex.eqtls$hg38)
            end   <- max(tbl.gtex.eqtls$hg38)
            width.kb <- round((1 + end - start)/1000, digits=2)
            return(list(start=start, end=end, width.kb=width.kb))
            },
        #------------------------------------------------------------
        getAMPAD.eqtl.genomicRegion = function(){
            start <- min(private$tbl.fimo$start)
            end   <- max(private$tbl.fimo$end)
            width.kb <- round((1 + end - start)/1000, digits=2)
            return(list(start=start, end=end, width.kb=width.kb))
            },
        #------------------------------------------------------------
        getGTEx.eqtl.tissues = function(){
            private$gtex.eqtl.tissues
            },
        #------------------------------------------------------------
        get.current.GTEx.eqtl.tissue = function(){
            private$current.tissue
            },
        #------------------------------------------------------------
        set.current.GTEx.eqtl.tissue = function(new.tissue){
            stopifnot(new.tissue %in% private$gtex.eqtl.tissues)
            private$current.tissue <- new.tissue
            },
        #------------------------------------------------------------
        run.tms = function(){
            current.region <- self$getStudyRegion()
            tbl.fimo.sub <- subset(tbl.fimo, start >= current.region$start & end <= current.region$end)
            private$tms <- TMS$new(private$trenaProject,
                                   private$targetGene,
                                   tbl.fimo.sub,
                                   private$tbl.oc,
                                   quiet=FALSE)

            private$tms$scoreFimoTFBS()   # chip, conservation, genehancer, genic annotations, distance to tss
            current.tissue <- self$get.current.GTEx.eqtl.tissue()
            private$mtx.rna <- private$etx$get.rna.matrix(current.tissue)
            private$tms$add.tf.mrna.correlations(private$mtx.rna, featureName="cor.all")
            private$tbl.tms <- private$tms$getTfTable()
            },
        #------------------------------------------------------------
        add.eqtls.toTmsTable = function(){
            tbl.tms <- private$tbl.tms
            gr.tms <- GRanges(tbl.tms)
            gr.ampad <- with(private$tbl.ampad.eqtls, GRanges(seqnames=chrom[1], IRanges(hg38)))
            tbl.ov <- as.data.frame(findOverlaps(gr.ampad, gr.tms))
            ampad.eqtl <- rep(FALSE, nrow(tbl.tms))
            if(nrow(tbl.ov) > 0)
                ampad.eqtl[unique(tbl.ov[,2])] <- TRUE
            private$tbl.tms$ampad.eqtl <- ampad.eqtl

            tbl.gtex.tissue.eqtls <- subset(private$tbl.gtex.eqtls, id==private$current.tissue)
            dim(tbl.gtex.tissue.eqtls)
            gr.gtex <- with(tbl.gtex.tissue.eqtls, GRanges(seqnames=sprintf("chr%s", chrom[1]), IRanges(hg38)))
            tbl.ov <- as.data.frame(findOverlaps(gr.gtex, gr.tms))
            dim(tbl.ov)
            gtex.eqtl <- rep(FALSE, nrow(tbl.tms))
            if(nrow(tbl.ov) > 0)
                gtex.eqtl[unique(tbl.ov[,2])] <- TRUE
            private$tbl.tms$gtex.eqtl <- gtex.eqtl
            },
        #------------------------------------------------------------
        get.tmsTable = function(){
            private$tbl.tms
            },
        #------------------------------------------------------------
        get.breaksTable = function(){
            private$tbl.breaks
            },
        #------------------------------------------------------------
        run.trena = function(tf.candidates){
            solvers=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest")
            solver <- EnsembleSolver(private$mtx.rna,
                                     targetGene=private$targetGene,
                                     candidateRegulators=tf.candidates,
                                     geneCutoff=1.0,
                                     solverNames=solvers)
            tbl.trena <- trena::run(solver)
            new.order <- order(tbl.trena$rfScore, decreasing=TRUE)
            tbl.trena <- tbl.trena[new.order,]
            rownames(tbl.trena) <- NULL
            tbl.trena$rfNorm <- tbl.trena$rfScore / max(tbl.trena$rfScore)
            tbl.trena <- cbind(tbl.trena[1], as.data.frame(lapply(tbl.trena[-1], function(x) round(x, digits=2))))
            tbl.trena$tfbs <- unlist(lapply(tbl.trena$gene,
                                            function(tf) length(grep(tf, private$tbl.tmsFiltered$tf))))
            private$tbl.trena <- tbl.trena

            },
        #------------------------------------------------------------
        get.trenaTable = function(){
           private$tbl.trena
           },
        #------------------------------------------------------------
        set.tmsFilteredTable = function(tbl.filtered){
           private$tbl.tmsFiltered <- tbl.filtered
           },
        #------------------------------------------------------------
        get.tmsFilteredTable = function(){
           private$tbl.tmsFiltered
           },
        #------------------------------------------------------------
        breakMotifs = function(tbl.trena, tbl.tms){
            tbl.tfbs <- subset(tbl.tms, tf %in% tbl.trena$gene)
            dim(tbl.tfbs)

            browser()
            gr.tfbs <- GRanges(tbl.tfbs[, c("chrom", "start", "end")])
            tbl.gtex.eqtls <- tfam$get.gtex.eqtls()
            tbl.eqtls.tmp <- tbl.gtex.eqtls[, c("chrom", "hg38", "hg38", "rsid")]
            colnames(tbl.eqtls.tmp) <- c("chrom", "start", "end", "rsid")
            tbl.eqtls.tmp$chrom <- sprintf("chr%s", tbl.eqtls.tmp$chrom)
            gr.eqtls <- GRanges(tbl.eqtls.tmp)
            tbl.ov <- as.data.frame(findOverlaps(gr.eqtls, gr.tfbs))
            dim(tbl.ov)
            tbl.eqtls.ov <- tbl.eqtls.tmp[tbl.ov$queryHits,]
            tbl.tfbs.ov <- tbl.tfbs[tbl.ov$subjectHits,]
            tfs.oi <- unique(tbl.tfbs.ov$tf)
            tfs.oi <- tfs.oi[match(tbl.trena$gene, tfs.oi)]
            deleters <- which(is.na(tfs.oi))
            if(length(deleters) > 0)
                tfs.oi <- tfs.oi[-deleters]
            rsids.oi <- unique(tbl.eqtls.ov$rsid)
            length(rsids.oi)

            mdb.human <- query(MotifDb, "sapiens", c("jaspar2022", "hocomoco-core-A"))
            motifs.selected <- query(mdb.human, andStrings="sapiens", orStrings=tfs.oi)
            message(sprintf("--- about to look up %d rsids, may take 5 minutes", length(rsids.oi)))
            x <- system.time(snps.gr <- snps.from.rsid(rsid = rsids.oi,
                                                       dbSNP=SNPlocs.Hsapiens.dbSNP155.GRCh38,
                                                       search.genome=BSgenome.Hsapiens.UCSC.hg38))
            message(sprintf("--- snps.gr for %d snps obtained, elapsed time: %f", length(rsids.oi), x[["elapsed"]]))


            bpparam <- MulticoreParam(workers=3)
            message(sprintf("--- calling motifbreakR"))

            results <- motifbreakR(snpList = snps.gr,
                                   filterp = TRUE,
                                   pwmList = motifs.selected,
                                   show.neutral=FALSE,
                                   method = c("ic", "log", "notrans")[1],
                                   bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                   BPPARAM = bpparam,
                                   verbose=TRUE)
            tbl.breaks <- as.data.frame(results, row.names=NULL)
            colnames(tbl.breaks)[1] <- "chrom"
            tbl.breaks$chrom <- as.character(tbl.breaks$chrom)
            tbl.breaks <- subset(tbl.breaks, effect=="strong")
            tbl.breaks$start <- tbl.breaks$start - 1
            tbl.breaks$pctDelta <- with(tbl.breaks, pctAlt-pctRef)
            new.order <- order(abs(tbl.breaks$pctDelta), decreasing=TRUE)
            tbl.breaks <- tbl.breaks[new.order,]
            motifbreakR.coi <- c("chrom", "start", "end", "SNP_id", "geneSymbol", "providerName", "pctRef", "pctAlt", "pctDelta")
            tbl.breaks <- tbl.breaks[, motifbreakR.coi]

                # motifbreakR will find fresh binding sites, possibly for genes for which we did not
                # so traverse these break regions, make sure they overlap with regions we found for
                # those genes
            gc <- sort(intersect(tbl.tms$tf, tbl.breaks$geneSymbol))
            tbls <- list()
            for(gene in gc){
                tbl.breaks.gene <- subset(tbl.breaks, geneSymbol==gene)
                gr.breaks <- GRanges(tbl.breaks.gene)
                tbl.tms.gene <- subset(tbl.tms, tf==gene)
                gr.tms <- GRanges(tbl.tms.gene)
                tbl.ov <- as.data.frame(findOverlaps(gr.breaks, gr.tms))
                tbl.new <- data.frame()
                if(nrow(tbl.ov) > 0){
                    tbl.new <- tbl.breaks.gene[unique(tbl.ov[,1]),]
                    dups <- which(duplicated(tbl.new[, c("chrom", "start", "SNP_id", "geneSymbol")]))
                    if(length(dups) > 0)
                        tbl.new <- tbl.new[-dups,]
                    }  # if tbl.ov
                tbls[[gene]] <- tbl.new
                } # for gene
            tbl.breaks.tfbs <- do.call(rbind, tbls)
            rownames(tbl.breaks.tfbs) <- NULL
            message(sprintf("--- saving tbl.breaks, %d x %d", nrow(tbl.breaks.tfbs), ncol(tbl.breaks.tfbs)))
            private$tbl.breaks <- tbl.breaks.tfbs
            },  # breakMotifs
        #------------------------------------------------------------
        viz = function(igv){
            shoulder <- 1000
            roi <- self$getStudyRegion()
            roi.string <- with(roi, sprintf("%s:%d-%d", chrom, start-1000, end+1000))
            showGenomicRegion(igv, roi.string)
            tbl.trena <- self$get.trenaTable()
            tbl.tms <- self$get.tmsFilteredTable()
            stopifnot(all(tbl.tms$tf %in% tbl.trena$gene))
                #-----------------------------------
                # draw the tfbs, one track per tf
                #-----------------------------------
            for(TF in tbl.trena$gene){
                tbl.track <- subset(tbl.tms, tf==TF)[, c("chrom", "start", "end")]
                track <- DataFrameAnnotationTrack(TF, tbl.track, color="random", trackHeight=25)
                displayTrack(igv, track)
                }
                #-----------------------------------
                # draw the open chromatin
                #-----------------------------------
            extra.margin <- 2000
            tbl.track <- subset(private$tbl.oc, chrom==roi$chrom &
                                                start >= roi$start-extra.margin &
                                                end <= roi$end+extra.margin)[, c("chrom", "start", "end")]
            track <- DataFrameAnnotationTrack("OC", tbl.track, color="darkgrey")
            displayTrack(igv, track)

                #----------------------------------------------------
                # gtex tissue-specific eqtls, alredy filterd for pval
                #-----------------------------------------------------
            tbl.eqtl <- self$get.gtex.eqtls()
            tbl.eqtl <- subset(tbl.eqtl, hg38 >= roi$start & hg38 <= roi$end)
            tbl.eqtl <- tbl.eqtl[, c("chrom", "hg38", "hg38", "rsid", "pvalue", "beta")]
            colnames(tbl.eqtl)[c(2,3)] <- c("start", "end")
            tbl.eqtl$start <- tbl.eqtl$start - 1
            tbl.eqtl$score <- with(tbl.eqtl, -log10(pvalue) * beta)
            coi <- c("chrom", "start", "end", "score", "rsid")
            tbl.eqtl <- tbl.eqtl[, coi]
            dups <- which(duplicated(tbl.eqtl[, c("chrom", "start", "end", "rsid")]))
            if(length(dups) > 0)
                tbl.eqtl <- tbl.eqtl[-dups,]
            browser()
            title <- self$get.current.GTEx.eqtl.tissue()
            title <- sub("Brain_", "", title)
            track <- DataFrameQuantitativeTrack(title, tbl.eqtl, autoscale=TRUE, color="black")
            displayTrack(igv, track)

                #----------------------------------------------------
                # motifbreakR results
                #-----------------------------------------------------
            tbl.breaks <- ts$get.breaksTable()
            if(!all(is.null(tbl.breaks))){
                tbl.breaks$score <- with(tbl.breaks, pctRef * pctDelta * 100)
                dups <- which(duplicated(tbl.breaks[, c("SNP_id", "geneSymbol")]))
                if(length(dups) > 0)
                    tbl.breaks <- tbl.breaks[-dups,]
                tbl.track <- tbl.breaks[, c("chrom", "start", "end", "score", "SNP_id", "geneSymbol")]
                track <- DataFrameQuantitativeTrack("motif disruption", tbl.track, colo="red", autoscale=TRUE)
                displayTrack(igv, track)
                }  # tbl.breaks

            } # viz
        #------------------------------------------------------------
       ) # public

    ) # class
#--------------------------------------------------------------------------------
