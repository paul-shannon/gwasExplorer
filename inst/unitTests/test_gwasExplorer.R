library(RUnit)
library(gwasExplorer)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()

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
    checkTrue(all(as.numeric(lapply(x, nrow)) > 5))   # ampad 48, gtex 8, probably due to sample size
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
if(!interactive())
    runTests()
