library(R6)
library(EndophenotypeExplorer)
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
                   mtx.rna=NULL
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
        setStudyRegion = function(start, end){
            width.kb <- round(((1 + end - start)/1000), digits=2)
            private$study.region <- list(start=start, end=end, width.kb=width.kb)
            },
        #------------------------------------------------------------
        getStudyRegion = function(){
            private$study.region
            },
        #------------------------------------------------------------
        getFimoGenomicRegion = function(){
            start <- min(private$tbl.fimo$start)
            end   <- max(private$tbl.fimo$end)
            width.kb <- round((1 + end - start)/1000, digits=2)
            return(list(start=start, end=end, width.kb=width.kb))
            },
        #------------------------------------------------------------
        get.gtex.eqtls <- function(){
            invisible(tbl.gtex.eqtls)
            },
        #------------------------------------------------------------
        get.ampad.eqtls <- function(){
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
           }
        #------------------------------------------------------------
       ) # public

    ) # class
#--------------------------------------------------------------------------------
