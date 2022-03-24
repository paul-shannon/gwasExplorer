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
                   tbl.eqtls.general=NULL,
                   etx=NULL,
                   avx=NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(
         #' @description
         #' Creates a new instance of this [R6][R6::R6Class] class.
         #' @param targetGene character
         #' @param locusName character, from GWAS studies, often meta-analyses
         #' @param tagSnp character
         #' @return a new instance of gwasExplorer
        initialize = function(targetGene, locusName, tagSnp){
            private$targetGene <- targetGene
            private$locusName <- locusName
            private$tagSnp <- tagSnp
            private$tbl.linkage <- self$loadLinkageTable(tagSnp)
            private$tbl.chromLocs <- self$getChromLocs(targetGene)
            private$etx <-  EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD",
                                                      initialize.snpLocs=FALSE)
            private$avx <- ADvariantExplorer$new(targetGene,
                                                 private$tbl.chromLocs$chrom,
                                                 private$tbl.chromLocs$start,
                                                 private$tbl.chromLocs$end)
            #private$tbl.eqtls.ampad <- etx$get.ampad.EQTLsForGene()

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
           tbl.eqtls.ampad <- private$etx$get.ampad.EQTLsForGene()
           tbl.eqtls.ampad <- subset(tbl.eqtls.ampad, pvalue <= pval.threshold)
           chrom <- private$tbl.locs$chrom
           start <- private$tbl.locs$start.hg38
           end <- private$tbl.locs$end.hg38
           start <- private$tbl.chromLocs$start.38
           end   <- private$tbl.chromLocs$end.38
           chrom <- private$tbl.chromLocs$chrom
           tbl.gtex <- private$avx$geteQTLsByLocationAndStudyID(chrom, start-shoulder, end+shoulder,
                                                                studyIDs=eqtl.catalog.studyIDs,
                                                                method="REST", simplify=TRUE)
           tbl.gtex <- subset(tbl.gtex, gene==private$targetGene & pvalue <= pval.threshold)
           tbl.gtex.locs <- private$etx$rsidToLoc(unique(tbl.gtex$rsid))
           tbl.gtex <- merge(tbl.gtex, tbl.gtex.locs, by="rsid")
           list(ampad=tbl.eqtls.ampad, gtex=tbl.gtex)
           }
       ) # public

    ) # class
#--------------------------------------------------------------------------------
