library(RUnit)
library(gwasExplorer)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_getChromLocs()
    test_eqtlsForGene()
    test_trenaMultiScore()
    test_trenaMultiScore_add.eQTLs()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    gwex <- gwasExplorer$new(targetGene="NDUFS2", locusName="ADAMTS4", tagSnp="rs4575098")
    checkTrue(all(c("R6", "gwasExplorer") %in% class(gwex)))

    tbl.linkage <- gwex$getLinkageTable()
    checkEquals(dim(tbl.linkage), c(78, 8))
    checkEquals(colnames(tbl.linkage),
                c("chrom", "hg19", "hg38", "rSquared", "dPrime", "rsid", "ref", "alt"))

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_getChromLocs <- function()
{
    message(sprintf("--- test_getChromLocs"))

    goi <- "NDUFS2"
    gwex <- gwasExplorer$new(targetGene=goi, locusName="ADAMTS4", tagSnp="rs4575098")
    tbl.locs <- gwex$getChromLocs(goi)
    checkEquals(dim(tbl.locs), c(1, 5))
    checkEquals(colnames(tbl.locs), c("chrom","start.19","end.19","start.38","end.38"))
    checkTrue(!is.factor(tbl.locs$chrom))

} # test_getChromLocs
#----------------------------------------------------------------------------------------------------
test_eqtlsForGene <- function()
{
    message(sprintf("--- test_eqtlsForGene"))

    gwex <- gwasExplorer$new(targetGene="NDUFS2", locusName="ADAMTS4", tagSnp="rs4575098")
    tbl.summary <- gwex$getEQTLSummary()
    checkTrue(nrow(tbl.summary) > 450)
    brain.tissue.study.ids <- grep("GTEx_V8.Brain_[CH]", tbl.summary$unique_id, v=TRUE)
    checkTrue(length(brain.tissue.study.ids) >= 6)

    pval.max <- 1e-5
    x <- gwex$getEqtlsForGene(shoulder=0, eqtl.catalog.studyIDs=brain.tissue.study.ids,
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

    shoulder <- 5000
    gwex <- gwasExplorer$new(targetGene="NDUFS2", locusName="ADAMTS4", tagSnp="rs4575098")
    tbl.fimo.ndufs2 <- get(load("~/github/TrenaProjectAD/prep/bigFimo/from-khaleesi/tbl.fimo.NDUFS2.RData"))
    data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
    filename <- "boca-hg38-consensus-ATAC.RData"
    tbl.boca <- get(load(file.path(data.dir, filename)))

    tbl.tms <- gwex$trenaMultiScore("GTEx_V8.Brain_Hippocampus", shoulder=shoulder,
                                    tbl.fimo=tbl.fimo.ndufs2, tbl.oc=tbl.boca)
    top.ranked <- table(subset(tbl.tms, abs(cor.all) > 0.5 & (chip | oc))$tf)

    top.tfs <- unique(subset(tbl.tms, abs(cor.all) > 0.5 & (chip | oc))$tf)
    checkTrue(all(c("ASCL1","EBF1","ZEB1","NFIA","SOX4","SOX21") %in% top.tfs))

    brain.tissue.study.ids <- "GTEx_V8.Brain_Hippocampus"

    pval.max <- 1e-3
    x <- gwex$getEqtlsForGene(shoulder=shoulder, eqtl.catalog.studyIDs=brain.tissue.study.ids,
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


    shoulder <- 5000
    gwex <- gwasExplorer$new(targetGene="NDUFS2", locusName="ADAMTS4", tagSnp="rs4575098")

       #-----------------------------------------------
       # find and save eQTLs, ampad and one gtex tissue
       #-----------------------------------------------

    brain.tissue.study.ids <- "GTEx_V8.Brain_Hippocampus"
    pval.max <- 1e-3
    x <- gwex$getEqtlsForGene(shoulder=shoulder, eqtl.catalog.studyIDs=brain.tissue.study.ids,
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

    tbl.tms <- gwex$trenaMultiScore("GTEx_V8.Brain_Hippocampus", shoulder=shoulder,
                                    tbl.fimo=tbl.fimo.ndufs2, tbl.oc=tbl.boca)
    checkTrue(all(c("ampad.eqtl", "gtex.eqtl") %in% colnames(tbl.tms)))

    top.ranked <- table(subset(tbl.tms, abs(cor.all) > 0.5 & (chip | oc | ampad.eqtl | gtex.eqtl))$tf)

    top.tfs <- names(top.ranked)
    checkTrue(all(c("ASCL1","EBF1","ZEB1","NFIA","SOX4","SOX21") %in% top.tfs))

} # test_trenaMultiScore_add.eQTLs
#----------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
