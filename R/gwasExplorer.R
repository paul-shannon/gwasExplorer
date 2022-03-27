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
                   shoulder=NULL,
                   tissueName=NULL,
                   tbl.linkage=NULL,
                   tbl.chromLocs=NULL,
                   tbl.eqtls.ampad=NULL,
                   tbl.eqtls.gtex=NULL,
                   tbl.tms=NULL,
                   tbl.trena=NULL,
                   tbl.trena.scored=NULL,
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
         #' @param shoulder numeric number of bases up and downstream to consider in lookups
         #' @param tissueName character  selects expression & variants
         #' @parma tbl.prespecifiedRegion data.frame, default NA, use this instead of gene tx +/- shoulder
         #' @return a new instance of gwasExplorer
        initialize = function(targetGene, locusName, tagSnp, shoulder, tissueName, tbl.prespecifiedRegion=NA){
            private$targetGene <- targetGene
            private$locusName <- locusName
            private$tagSnp <- tagSnp
            private$shoulder <- shoulder
            private$tissueName <- tissueName

            private$tbl.linkage <- self$loadLinkageTable(tagSnp)
            tbl.locs <- self$getChromLocs(targetGene)
            if(!is.data.frame(tbl.prespecifiedRegion)){

               }
            if(is.data.frame(tbl.prespecifiedRegion)){
               tbl.locs <- tbl.prespecifiedRegion
               }
            private$tbl.chromLocs <- tbl.locs
            private$etx <-  EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD",
                                                      initialize.snpLocs=FALSE)
            private$avx <- ADvariantExplorer$new(targetGene,
                                                 tbl.locs$chrom,
                                                 tbl.locs$start,
                                                 tbl.locs$end)

            private$tbl.chromLocs$width.hg19 <- with(private$tbl.chromLocs, 1+end.hg19-start.hg19)
            private$tbl.chromLocs$width.hg38 <- with(private$tbl.chromLocs, 1+end.hg38-start.hg38)

            loc.chrom <- private$tbl.chromLocs$chrom
            loc.start <- private$tbl.chromLocs$start.hg38 - private$shoulder
            loc.end   <- private$tbl.chromLocs$end.hg38   + private$shoulder


            private$gls <- GwasLocusScores$new(private$tagSnp,
                                              loc.chrom, loc.start, loc.end,
                                              #tissue.name=private$tissueName,
                                              targetGene=private$targetGene)

            },
        #------------------------------------------------------------
        #' @description the regulatory region being studied
        #' @return data.frame
        getRegion = function(){
            private$tbl.chromLocs
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
        #' @param geneSymbol character, a HUGO human gene symbol
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
                      start.hg19=min(tbl.tx.hg19$TXSTART),
                      end.hg19=min(tbl.tx.hg19$TXEND),
                      start.hg38=min(tbl.tx.hg38$TXSTART),
                      end.hg38=min(tbl.tx.hg38$TXEND),
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
        #' @param eqtl.catalog.studyIDs character vector,
        #' @param pval.threshold numeric retain only eQTLs as or more significant than this
        #' @return data.frame
        getEqtlsForGene = function(eqtl.catalog.studyIDs, pval.threshold){
           chrom <- private$tbl.chromLocs$chrom
           start <- private$tbl.chromLocs$start.hg38 - private$shoulder
           end   <- private$tbl.chromLocs$end.hg38 + private$shoulder

           tbl.eqtls.ampad <- private$etx$get.ampad.EQTLsForGene()
           tbl.eqtls.ampad <- subset(tbl.eqtls.ampad, pvalue <= pval.threshold)
           tbl.eqtls.ampad <- subset(tbl.eqtls.ampad, hg38 >= start & hg38 <= end)

           tbl.gtex <- private$avx$geteQTLsByLocationAndStudyID(chrom, start, end,
                                                                studyIDs=eqtl.catalog.studyIDs,
                                                                simplify=TRUE)
           tbl.gtex <- subset(tbl.gtex, gene==private$targetGene & pvalue <= pval.threshold)
           if(nrow(tbl.gtex) > 0){
              tbl.gtex.locs <- private$etx$rsidToLoc(unique(tbl.gtex$rsid))
              tbl.eqtls.gtex <- merge(tbl.gtex, tbl.gtex.locs, by="rsid")
              }
           if(nrow(tbl.gtex) == 0){
              tbl.eqtls.gtex <- data.frame()
              }
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
        trenaMultiScore = function(tissueName, tbl.fimo, tbl.oc){
           mtx.rna <- private$etx$get.rna.matrix(tissueName)
           dim(mtx.rna)
           loc.chrom <- private$tbl.chromLocs$chrom
           loc.start <- private$tbl.chromLocs$start.hg38 - private$shoulder
           loc.end   <- private$tbl.chromLocs$end.hg38   + private$shoulder
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
           }, # trenaMultiScore

        #' @description getExpressionMatrixCodes: short names (codes) and file paths
        #'
        #' @return list
        getExpressionMatrixCodes = function(){
            private$etx$get.rna.matrix.codes()
            },

        #' @description getExpressionMatrix
        #' @param code character, a short descriptive name, resolves to a file path
        #'
        #' @return data.frame
        getExpressionMatrix = function(code){
            mtx.rna <- private$etx$get.rna.matrix(code)
            invisible(mtx.rna)
            },

        #' @description run trena
        #' @param tfs character vector, one or more transcription factor gene symbols
        #' @param mtx.rna matrix, gene expression of many genes across many conditions or samples
        #' @param tbl.tms data.frame, tf binding sites, annotated and scored in many ways
        #'
        #' @return data.frame
         runTrena = function(tfs, mtx.rna, tbl.tms){
           tbl.trena <- private$gls$runTrena(tfs, mtx.rna, tbl.tms)
           tbl.trena$rfNorm <- tbl.trena$rfScore / max(tbl.trena$rfScore)
           private$tbl.trena <- tbl.trena
           tbl.trena
         },

        #' @description score broken motifs
        #' @param max.pvalue numeric
        #' @param tbl.eqtls data.frame  rsids used as input to motifbreakR
        scoreBrokenMotifs = function(max.pvalue, tbl.eqtls){
            private$gls$set.eqtls(tbl.eqtls)
            x <- system.time(tbl.breaks <- private$gls$breakMotifsAtEQTLs(private$targetGene,
                                                                          pvalue.cutoff=max.pvalue))
            print(x[["elapsed"]])
            tbl.trena.scored <- private$gls$scoreEQTLs(private$tbl.trena, tbl.breaks, private$tbl.eqtls.gtex)
            tbl.trena <- private$tbl.trena
            tbl.eqtls.gtex <- private$tbl.eqtls.gtex
            private$tbl.trena.scored <- tbl.trena.scored
            save(tbl.trena, tbl.breaks, tbl.eqtls.gtex, tbl.trena.scored, file="scoreEQTLs.RData")
            tbl.trena.scored
            }

       ) # public

    ) # class
#--------------------------------------------------------------------------------
