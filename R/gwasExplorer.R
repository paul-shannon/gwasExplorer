#' @title gwasExplorer
#' @description score non-coding variants in GWAS loci to elucidate function
#' @name gwasExplorer
#'
#' @import R6
#' @import EndophenotypeExplorer
#' @import GwasLocusScores
#' @import ADvariantExplorer
#' @import org.Hs.eg.db
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import TrenaProjectAD
#'
#' @examples
#'   gwex <- gwasExplorer$new(targetGene="NDUFS2", locusName="ADAMTS4", tagSnp="rs4575098")
#' @export

gwasExplorer = R6Class("gwasExplorer",

    #--------------------------------------------------------------------------------
    private = list(targetGene=NULL,
                   entrezID=NULL,
                   locusName=NULL,
                   tagSnp=NULL,
                   tbl.linkage=NULL,
                   tbl.chromLocs=NULL,
                   tbl.eqtls.ampad=NULL,
                   tbl.eqtls.gtex=NULL,
                   tbl.tms=NULL,
                   etx=NULL,
                   avx=NULL,
                   gls=NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(
         #' @description
         #' Creates a new instance of this [R6][R6::R6Class] class.
         #' @param targetGene character
         #' @param locusName haracter, from GWAS studies, often meta-analyses
         #' @param tagSnp character
         #' @return a new instance of gwasExplorer
        initialize = function(targetGene, locusName, tagSnp){
            private$targetGene <- targetGene
            private$locusName <- locusName
            private$tagSnp <- tagSnp
            private$tbl.linkage <- self$loadLinkageTable(tagSnp)
            tbl.locs <- self$getChromLocs(targetGene)
            private$tbl.chromLocs <- tbl.locs
            private$etx <-  EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD",
                                                      initialize.snpLocs=FALSE)
            private$avx <- ADvariantExplorer$new(targetGene,
                                                 tbl.locs$chrom,
                                                 tbl.locs$start,
                                                 tbl.locs$end)

            },
        #------------------------------------------------------------
        #' @description from halporeg, a table with both hg19 and hg38 coordinates
        #' @return data.frame
        getLinkageTable = function(){
            private$tbl.linkage
            },
        #------------------------------------------------------------
        #' @description from halporeg, a table with both hg19 and hg38 coordinates
        #' @param tagSnp character
        #' @return data.frame
        loadLinkageTable = function(tagSnp){
           stopifnot(tagSnp %in% c("rs4575098"))
           filename <- switch(tagSnp,
                              "rs4575098" = system.file(package="gwasExplorer", "extdata",
                                                        "adamts4.locus",
                                                        "haploreg-rs4575098-0.2-hg19-hg38.RData")
                              )
           stopifnot(file.exists(filename))
           tbl <- get(load(filename))
           return(tbl)
           },
        #------------------------------------------------------------
        #' @description hg19 and hg38 chromosomal locations
        #' @param geneSymbol
        #' @return data.frame
        getChromLocs = function(geneSymbol){
           suppressMessages({
             geneID <- AnnotationDbi::select(org.Hs.eg.db, keys=geneSymbol, keytype="SYMBOL", columns=c("ENTREZID"))$ENTREZID
             tbl.tx.hg19 <- AnnotationDbi::select(TxDb.Hsapiens.UCSC.hg19.knownGene, keys=geneID, keytype="GENEID",
                                   columns=c("TXCHROM", "TXSTART", "TXEND"))
             tbl.tx.hg38 <- AnnotationDbi::select(TxDb.Hsapiens.UCSC.hg38.knownGene, keys=geneID, keytype="GENEID",
                                   columns=c("TXCHROM", "TXSTART", "TXEND"))
             })
           data.frame(chrom=tbl.tx.hg19$TXCHROM[1],
                      start.19=min(tbl.tx.hg19$TXSTART),
                      end.19=min(tbl.tx.hg19$TXEND),
                      start.38=min(tbl.tx.hg38$TXSTART),
                      end.38=min(tbl.tx.hg38$TXEND),
                      stringsAsFactors=FALSE)
            },
        #------------------------------------------------------------
        #' @description retrieve a summary of eQTL studies currently available from embl
        #' return data.frame
        getEQTLSummary = function(){
            private$avx$geteQTLSummary()
            },

        #------------------------------------------------------------
        #' @description from all available sources, a table with both hg19 and hg38 coordinates
        #' @param shoulder numeric number of bytes up- and downstream from the largst transcript
        #' @param eqtl.catalog.studyIDs character vector,
        #' @param pval.threshold numeric retain only eQTLs as or more significant than this
        #' @return data.frame
        getEqtlsForGene = function(shoulder, eqtl.catalog.studyIDs, pval.threshold){
           chrom <- private$tbl.chromLocs$chrom
           start <- private$tbl.chromLocs$start.38 - shoulder
           end   <- private$tbl.chromLocs$end.38 + shoulder

           tbl.eqtls.ampad <- private$etx$get.ampad.EQTLsForGene()
           tbl.eqtls.ampad <- subset(tbl.eqtls.ampad, pvalue <= pval.threshold)
           tbl.eqtls.ampad <- subset(tbl.eqtls.ampad, hg38 >= start & hg38 <= end)

           tbl.gtex <- private$avx$geteQTLsByLocationAndStudyID(chrom, start, end,
                                                                studyIDs=eqtl.catalog.studyIDs,
                                                                method="REST", simplify=TRUE)
           tbl.gtex <- subset(tbl.gtex, gene==private$targetGene & pvalue <= pval.threshold)
           tbl.gtex.locs <- private$etx$rsidToLoc(unique(tbl.gtex$rsid))
           tbl.eqtls.gtex <- merge(tbl.gtex, tbl.gtex.locs, by="rsid")
           private$tbl.eqtls.gtex  <- tbl.eqtls.gtex
           private$tbl.eqtls.ampad <- tbl.eqtls.ampad
           list(ampad=tbl.eqtls.ampad, gtex=tbl.eqtls.gtex)
           },

        #' @description build tms model
        #' @param tissueName character as used by GTEx
        #' @param shoulder numeric number of bytes up- and downstream from the largst transcript
        #' @param tbl.fimo data.frame, typically a very large & low threshold set of fimo scores
        #' @param tbl.oc data.frame, open chromatin from, e.g., ATAC-seq, as tissue- and cell-type specific as possible
        #'
        #' @return data.frame
        trenaMultiScore = function(tissueName, shoulder, tbl.fimo, tbl.oc){
           loc.chrom <- private$tbl.chromLocs$chrom
           loc.start <- private$tbl.chromLocs$start.38 - shoulder
           loc.end   <- private$tbl.chromLocs$end.38   + shoulder
           private$gls <- GwasLocusScores$new(private$tagSnp,
                                              loc.chrom, loc.start, loc.end,
                                              tissue.name=tissueName,
                                              targetGene=private$targetGene)
           mtx.rna <- private$etx$get.rna.matrix(tissueName)
           dim(mtx.rna)
           tbl.fimo.sub <- subset(tbl.fimo, start >= loc.start & end <= loc.end)
           tbl.oc.sub <- subset(tbl.oc, start >= loc.start & end <= loc.end)
           tbl.tms <- private$gls$createTrenaMultiScoreTable(TrenaProjectAD(),
                                                             tbl.fimo.sub,
                                                             tbl.oc.sub,
                                                             mtx.rna)
           gr.tms <- GRanges(tbl.tms)
           if(is.data.frame(private$tbl.eqtls.ampad)){
              gr.ampad <- with(private$tbl.eqtls.ampad, GRanges(seqnames=chrom, IRanges(hg38)))
              tbl.ov <- as.data.frame(findOverlaps(gr.ampad, gr.tms))
              ampad.eqtl <- rep(FALSE, nrow(tbl.tms))
              if(nrow(tbl.ov) > 0)
                  ampad.eqtl[unique(tbl.ov[,2])] <- TRUE
              tbl.tms$ampad.eqtl <- ampad.eqtl
              }
           if(is.data.frame(private$tbl.eqtls.gtex)){
              gr.gtex <- with(private$tbl.eqtls.gtex, GRanges(seqnames=paste0("chr", chrom), IRanges(hg38)))
              tbl.ov <- as.data.frame(findOverlaps(gr.gtex, gr.tms))
              gtex.eqtl <- rep(FALSE, nrow(tbl.tms))
              if(nrow(tbl.ov) > 0)
                  gtex.eqtl[unique(tbl.ov[,2])] <- TRUE
              tbl.tms$gtex.eqtl <- gtex.eqtl
              }
           private$tbl.tms <- tbl.tms
           return(invisible(private$tbl.tms))
           } # trenaMultiScore

       ) # public

    ) # class
#--------------------------------------------------------------------------------
