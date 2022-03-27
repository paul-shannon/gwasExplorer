library(RUnit)
library(gwasExplorer)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_setGetRegionToStudy
    test_eqtlsForGene()
    test_trenaMultiScore()
    test_trenaMultiScore_add.eQTLs()
    test_runTrena()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    gwex <- gwasExplorer$new(targetGene="NDUFS2", locusName="ADAMTS4", tagSnp="rs4575098",
                             shoulder=1000, tissueName="GTEx_V8.Brain_Hippocampus")
    checkTrue(all(c("R6", "gwasExplorer") %in% class(gwex)))

    tbl.linkage <- gwex$getLinkageTable()
    checkEquals(dim(tbl.linkage), c(78, 8))
    checkEquals(colnames(tbl.linkage),
                c("chrom", "hg19", "hg38", "rSquared", "dPrime", "rsid", "ref", "alt"))

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_setGetRegionToStudy <- function()
{
   tbl.tagSnp <- data.frame(chrom="chr1",
                            hg19=161155392,
                            hg38=161185602,
                            rSquared=1,
                            dPrime=1,
                            rsid="rs4575098",
                            ref="G",
                            alt="A",
                            stringsAsFactors=FALSE)
   genomic.shoulder <- 10

   tbl.roi <- data.frame(chrom="1",
                         start.hg38=tbl.tagSnp$hg38-genomic.shoulder,
                         end.hg38=tbl.tagSnp$hg38+genomic.shoulder,
                         start.hg19=tbl.tagSnp$hg19-genomic.shoulder,
                         end.hg19=tbl.tagSnp$hg19+genomic.shoulder,
                         stringsAsFactors=FALSE)

   brain.tissue <- "GTEx_V8.Brain_Hippocampus"
   gene.expression.absCor.minimum <- 0.2
   eqtl.pval.max <- 1e-4
   targetGene <- "NDUFS2"
   gwex <- gwasExplorer$new(targetGene=targetGene, locusName="ADAMTS4", tagSnp="rs4575098",
                            shoulder=genomic.shoulder, tissueName=brain.tissue,
                            tbl.prespecifiedRegion=tbl.roi)

   tbl.roi2 <- gwex$getRegion()

   with(tbl.roi2, {
       checkEquals(width.hg19, 21);
       checkEquals(width.hg38, 21);

       checkEquals(start.hg19, 161155382);
       checkEquals(end.hg19,   161155402);

       checkEquals(start.hg38, 161185592);
       checkEquals(end.hg38,   161185612);
       })

   gwex <- gwasExplorer$new(targetGene=targetGene, locusName="ADAMTS4", tagSnp="rs4575098",
                            shoulder=genomic.shoulder, tissueName=brain.tissue,
                            tbl.prespecifiedRegion=NA)
   tbl.roi3 <- gwex$getRegion()

   with(tbl.roi3, {
       checkEquals(width.hg19, 9999);
       checkEquals(width.hg38, 6730);

       checkEquals(start.hg19, 161169105);
       checkEquals(end.hg19,   161179103);

       checkEquals(start.hg38, 161197104);
       checkEquals(end.hg38,   161203833);
       })

} # test_setGetRegionToStudy
#----------------------------------------------------------------------------------------------------
test_eqtlsForGene <- function()
{
    message(sprintf("--- test_eqtlsForGene"))

    gwex <- gwasExplorer$new(targetGene="NDUFS2", locusName="ADAMTS4", tagSnp="rs4575098",
                             shoulder=1000, tissueName="GTEx_V8.Brain_Hippocampus")

    tbl.summary <- gwex$getEQTLSummary()
    checkTrue(nrow(tbl.summary) > 450)
    brain.tissue.study.ids <- grep("GTEx_V8.Brain_[CH]", tbl.summary$unique_id, v=TRUE)
    checkTrue(length(brain.tissue.study.ids) >= 6)

    pval.max <- 1e-5
    x <- gwex$getEqtlsForGene(eqtl.catalog.studyIDs=brain.tissue.study.ids,
                              pval.threshold=pval.max)
    checkEquals(names(x), c("ampad", "gtex"))
    checkTrue(all(as.numeric(lapply(x, nrow)) > 5))   # ampad 7, gtex 8, probably due to sample size
    tbl.ampad <- x$ampad
    tbl.gtex  <- x$gtex

    lapply(x, dim)
    checkEquals(colnames(tbl.ampad),
                c("chrom","hg19","hg38","rsid","pvalue","ensg","genesymbol","study","tissue","assay"))
    checkEquals(colnames(tbl.gtex),
                c("rsid","pvalue","gene","total.alleles","beta","id","chrom","hg19","hg38"))
    shared.rsids <- intersect(tbl.ampad$rsid, tbl.gtex$rsid)
    checkTrue(length(shared.rsids) >= 4)  # 4 on 24 mar 2022

    checkTrue(max(tbl.ampad$pvalue) < pval.max)
    checkTrue(max(tbl.gtex$pvalue)  < pval.max)

} # test_eqtlsForGene
#----------------------------------------------------------------------------------------------------
test_trenaMultiScore <- function()
{
    message(sprintf("--- test_trenaMultiScore"))

    gwex <- gwasExplorer$new(targetGene="NDUFS2", locusName="ADAMTS4", tagSnp="rs4575098",
                             shoulder=5000, tissueName="GTEx_V8.Brain_Hippocampus")
    tbl.fimo.ndufs2 <- get(load("~/github/TrenaProjectAD/prep/bigFimo/from-khaleesi/tbl.fimo.NDUFS2.RData"))
    data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
    filename <- "boca-hg38-consensus-ATAC.RData"
    tbl.boca <- get(load(file.path(data.dir, filename)))

    tbl.tms <- gwex$trenaMultiScore("GTEx_V8.Brain_Hippocampus",
                                    tbl.fimo=tbl.fimo.ndufs2, tbl.oc=tbl.boca)
    top.ranked <- table(subset(tbl.tms, abs(cor.all) > 0.5 & (chip | oc))$tf)

    top.tfs <- unique(subset(tbl.tms, abs(cor.all) > 0.5 & (chip | oc))$tf)
    checkTrue(all(c("ASCL1","EBF1","ZEB1","NFIA","SOX4","SOX21") %in% top.tfs))

    brain.tissue.study.ids <- "GTEx_V8.Brain_Hippocampus"

    pval.max <- 1e-3
    x <- gwex$getEqtlsForGene(eqtl.catalog.studyIDs=brain.tissue.study.ids,
                              pval.threshold=pval.max)
    tbl.ampad <- x$ampad
    tbl.gtex  <- x$gtex

    dim(tbl.ampad)
    dim(tbl.gtex)
    gr.ampad <- with(tbl.ampad, GRanges(seqnames=chrom, IRanges(hg38)))
    gr.gtex <- with(tbl.gtex, GRanges(seqnames=paste0("chr", chrom), IRanges(hg38)))
    gr.tms <- GRanges(tbl.tms[, 1:3])

    tbl.ov.1 <- as.data.frame(findOverlaps(gr.ampad, gr.tms))
    tbl.ov.2 <- as.data.frame(findOverlaps(gr.gtex, gr.tms))
    dim(tbl.ov.1)
    dim(tbl.ov.2)


} # test_trenaMultiScore
#----------------------------------------------------------------------------------------------------
# historically, the tms table does not mark fimo regions for eQTLs.
# but here, now, if eQTLs are called first, the will be added as TRUE/FALSE columns to the tms table
test_trenaMultiScore_add.eQTLs <- function()
{
    message(sprintf("--- test_trenaMultiScore_add.eQTLs"))

    gwex <- gwasExplorer$new(targetGene="NDUFS2", locusName="ADAMTS4", tagSnp="rs4575098",
                             shoulder=5000, tissueName="GTEx_V8.Brain_Hippocampus")

       #-----------------------------------------------
       # find and save eQTLs, ampad and one gtex tissue
       #-----------------------------------------------

    brain.tissue.study.ids <- "GTEx_V8.Brain_Hippocampus"
    pval.max <- 1e-3
    x <- gwex$getEqtlsForGene(eqtl.catalog.studyIDs=brain.tissue.study.ids,
                              pval.threshold=pval.max)

    unique.variants <- unique(c(x$ampad$rsid, x$gtex$rsid))
    length(unique.variants)   # 7

       #---------------------------------------------------------
       # now create the tms table. 2 eqtl columns should appear
       #---------------------------------------------------------

    tbl.fimo.ndufs2 <- get(load("~/github/TrenaProjectAD/prep/bigFimo/from-khaleesi/tbl.fimo.NDUFS2.RData"))
    data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
    filename <- "boca-hg38-consensus-ATAC.RData"
    tbl.boca <- get(load(file.path(data.dir, filename)))

    tbl.tms <- gwex$trenaMultiScore("GTEx_V8.Brain_Hippocampus",
                                    tbl.fimo=tbl.fimo.ndufs2, tbl.oc=tbl.boca)
    checkTrue(all(c("ampad.eqtl", "gtex.eqtl") %in% colnames(tbl.tms)))

    top.ranked <- table(subset(tbl.tms, abs(cor.all) > 0.5 & (chip | oc | ampad.eqtl | gtex.eqtl))$tf)

    top.tfs <- names(top.ranked)
    checkTrue(all(c("ASCL1","EBF1","ZEB1","NFIA","SOX4","SOX21") %in% top.tfs))

} # test_trenaMultiScore_add.eQTLs
#----------------------------------------------------------------------------------------------------
test_runTrena <- function()
{
    message(sprintf("--- test_runTrena"))

    gwex <- gwasExplorer$new(targetGene="NDUFS2", locusName="ADAMTS4", tagSnp="rs4575098",
                             shoulder=1000, tissueName="GTEx_V8.Brain_Hippocampus")
    file <- system.file(package="gwasExplorer", "extdata", "tbl.tms.ndufs2.RData")
    checkTrue(file.exists(file))
    tbl.tms <- get(load(file))

    tbl.sub <-
        subset(tbl.tms,
               abs(cor.all) > 0.5 &
               (chip | oc | ampad.eqtl | gtex.eqtl) &
               gh > 10 &
               fimo_pvalue < 1e-3)
    table(tbl.sub$tf)

    tfs <- unique(tbl.sub$tf)
    length(tfs)  # 6

    # gwex$getExpressionMatrixCodes()
    mtx.rna <- gwex$getExpressionMatrix("GTEx_V8.Brain_Hippocampus")
    tbl.trena <- gwex$runTrena(tfs, mtx.rna, tbl.sub)
    checkEquals(dim(tbl.trena), c(6, 11))

    checkEquals(colnames(tbl.trena),
                c("gene","betaLasso","betaRidge","spearmanCoeff","pearsonCoeff","rfScore",
                  "xgboost","rank","target","tfbs", "rfNorm"))

    brain.tissue.study.ids <- "GTEx_V8.Brain_Hippocampus"
    pval.max <- 1e-3
    x <- gwex$getEqtlsForGene(eqtl.catalog.studyIDs=brain.tissue.study.ids,
                              pval.threshold=pval.max)
    tbl.scored <- gwex$scoreBrokenMotifs(max.pvalue=1e-3, x$gtex)
    checkTrue(sum(tbl.scored$breakage.score) > 10)

} # test_runTrena
#----------------------------------------------------------------------------------------------------
full.run <- function()
{
   f <- system.file(package="gwasExplorer", "extdata", "adamts4.locus",
                    "haploreg-rs4575098-0.2-hg19-hg38.RData")
   file.exists(f)
   tbl.tagHap <- get(load(f))
   tbl.tagSnp <- subset(tbl.tagHap, rSquared == 1)
   genomic.shoulder <- 1000

   tbl.roi <- data.frame(chrom="1",
                         start.hg38=tbl.tagSnp$hg38-genomic.shoulder,
                         end.hg38=tbl.tagSnp$hg38+genomic.shoulder,
                         start.hg19=tbl.tagSnp$hg19-genomic.shoulder,
                         end.hg19=tbl.tagSnp$hg19+genomic.shoulder,
                         stringsAsFactors=FALSE)

   tbl.roi <- data.frame(chrom=sub("chr", "", tbl.tagHap$chrom[1]),
                         start.hg38=min(tbl.tagHap$hg38) - 1000,
                         end.hg38=max(tbl.tagHap$hg38) + 1000,
                         start.hg19=min(tbl.tagHap$hg19) - 1000,
                         end.hg19=max(tbl.tagHap$hg19) + 1000,
                         stringsAsFactors=FALSE)


   brain.studies <- c("GTEx_V8.Brain_Caudate_basal_ganglia",
                      "GTEx_V8.Brain_Cerebellar_Hemisphere",
                      "GTEx_V8.Brain_Cerebellum",
                      "GTEx_V8.Brain_Cortex",
                      "GTEx_V8.Brain_Hippocampus",
                      "GTEx_V8.Brain_Hypothalamus")
   gene.expression.absCor.minimum <- 0.2
   eqtl.pval.max <- 1e-2
   targetGene <- "NDUFS2"
   brain.tissue <- brain.studies[5]
   gwex <- gwasExplorer$new(targetGene=targetGene, locusName="ADAMTS4", tagSnp="rs4575098",
                            shoulder=genomic.shoulder, tissueName=brain.tissue,
                            tbl.prespecifiedRegion=tbl.roi)

   gwex$getRegion()

   x <- gwex$getEqtlsForGene(eqtl.catalog.studyIDs=brain.tissue, pval.threshold=eqtl.pval.max)
   lapply(x, dim)
   unique.variants <- unique(c(x$ampad$rsid, x$gtex$rsid))
   length(unique.variants)   # 7

      #---------------------------------------------------------
      # now create the tms table. 2 eqtl columns should appear
      #---------------------------------------------------------

   tbl.fimo <- get(load(sprintf("~/github/TrenaProjectAD/prep/bigFimo/from-khaleesi/tbl.fimo.%s.RData", targetGene)))
   data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
   filename <- "boca-hg38-consensus-ATAC.RData"
   tbl.boca <- get(load(file.path(data.dir, filename)))

   tbl.tms <- gwex$trenaMultiScore(brain.tissue, tbl.fimo=tbl.fimo, tbl.oc=tbl.boca)
   print(colnames(tbl.tms))

   tbl.tms.sub <- subset(tbl.tms, abs(cor.all) > gene.expression.absCor.minimum & (chip | (oc & (ampad.eqtl | gtex.eqtl))))

   top.ranked <- table(tbl.tms.sub$tf)
   print(top.ranked)
   top.tfs <- names(top.ranked)
   length(top.tfs)  # 6

   mtx.rna <- gwex$getExpressionMatrix(brain.tissue)
   tbl.trena <- gwex$runTrena(top.tfs, mtx.rna, tbl.tms.sub)
   tbl.scored <- gwex$scoreBrokenMotifs(max.pvalue=eqtl.pval.max, x$gtex)
   printf("------- tbl.scored, sum top 10: %f", sum(head(tbl.scored$breakage.score, n=10)))
   filename = sprintf ('gwex.%s.%s.%s.RData', targetGene, brain.tissue, format (Sys.time(), "%a.%b.%d.%Y-%H:%M:%S"))
   printf("--- saving to %s", filename)
   save(mtx.rna, tbl.tms, tbl.trena, tbl.scored, file=filename)

} # full.run
#----------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
