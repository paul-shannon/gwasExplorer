library(RUnit)
library(TrenaProjectAD)
source("tfam-R6.R")
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()
   test_setStudyRegion()
   test_tissue()
   test_run.tms()
   test_run.trena()
   test_run.small()
   test_run.large()

} # runTests
#----------------------------------------------------------------------------------------------------
targetGene <- "TFAM"
trenaProject <- TrenaProjectAD()
data.dir <- "~/github/gwasExplorer/studies/tfam/data"
if(!exists("tbl.fimo")){
   tbl.fimo <- get(load(file.path(data.dir, "tbl.fimo.TFAM.RData")))
   tbl.ampad.eqtls.raw <- get(load("data/ampad.eqtls.TFAM.RData"))
   tbl.ampad.eqtls <- subset(tbl.ampad.eqtls.raw, study=="ampad-rosmap")
   tbl.gtex.eqtls.raw <- get(load("data/tbl.eqtls.6.tissues-chr10:57385407-59385407.2022-06-08.08:03:52"))
   tbl.gtex.eqtls <- subset(tbl.gtex.eqtls.raw, gene==targetGene)
   data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
   filename <- "boca-hg38-consensus-ATAC.RData"
   tbl.oc <- get(load(file.path(data.dir, filename)))
   }

tfam <- TFAM$new(targetGene, trenaProject, tbl.fimo, tbl.gtex.eqtls, tbl.ampad.eqtls, tbl.oc)
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    checkTrue(all(c("R6", "TFAM") %in% class(tfam)))

    region <- tfam$getFimoGenomicRegion()
    checkEquals(region$start, 57391973)
    checkEquals(region$end,   59392030)
    checkEquals(region$width.kb, 2000.06)

    region <- tfam$getGTEx.eqtl.genomicRegion()
    checkEquals(region$start, 57385692)
    checkEquals(region$end,   59384547)
    checkEquals(region$width.kb, 1998.86)

    region <- tfam$getAMPAD.eqtl.genomicRegion()
    checkEquals(region$start, 57391973)
    checkEquals(region$end,   59392030)
    checkEquals(region$width.kb, 2000.06)

    study.region <- tfam$getStudyRegion()
    checkEquals(study.region$start, 57391973)
    checkEquals(study.region$end, 59392030)

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_setStudyRegion <- function()
{
    message(sprintf("--- test_setStudyRegion"))

    study.region <- tfam$getStudyRegion()
    checkEquals(study.region$start, 57391973)
    checkEquals(study.region$end,   59392030)

    new.start <- 58384216
    new.end   <- 58390229
    tfam$setStudyRegion(start=new.start, end=new.end)
    new.region <- tfam$getStudyRegion()

    checkEquals(new.region$start, new.start)
    checkEquals(new.region$end, new.end)
    checkEquals(new.region$width.kb, 6.01)

} # test_setStudRegion
#----------------------------------------------------------------------------------------------------
test_tissue <- function()
{
    message(sprintf("--- test_tissue"))
    tissues <- tfam$getGTEx.eqtl.tissues()
    checkEquals(length(tissues), 6)
    for(tissue in tissues){
       tfam$set.current.GTEx.eqtl.tissue(tissue)
       checkEquals(tfam$get.current.GTEx.eqtl.tissue(), tissue)
       } # for tissue

    TRUE

} # test_tissue
#----------------------------------------------------------------------------------------------------
test_run.tms <- function()
{
    message(sprintf("--- test_run.tms"))

    tissues <- tfam$getGTEx.eqtl.tissues()
    tfam$set.current.GTEx.eqtl.tissue(tissues[1])

    new.start <- 58385244
    new.end   <- 58385416
    tfam$setStudyRegion(start=new.start, end=new.end)

    tfam$run.tms()
    tbl.tms <- tfam$get.tmsTable()
    checkEquals(ncol(tbl.tms), 16)
    checkTrue(nrow(tbl.tms) > 500 & nrow(tbl.tms) < 650)
    new.cols <- setdiff(colnames(tbl.tms), colnames(tbl.fimo))
    checkTrue(all(c("chip","phast7","phast100","gh","oc","tss","cor.all") %in% new.cols))

    tfam$add.eqtls.toTmsTable()
    tbl.tms <- tfam$get.tmsTable()
    checkEquals(ncol(tbl.tms), 18)
    checkEquals(colnames(tbl.tms)[17:18], c("ampad.eqtl", "gtex.eqtl"))
    checkTrue(as.numeric(table(tbl.tms$ampad.eqtl))[2] > 60)
    checkTrue(as.numeric(table(tbl.tms$gtex.eqtl))[2] > 100)

} # test_run.tms
#----------------------------------------------------------------------------------------------------
test_run.trena <- function()
{
    message(sprintf("--- test_run.trena"))
    tbl.tms <- tfam$get.tmsTable()
    dim(tbl.tms)
    tbl.tms.filtered <- subset(tbl.tms, ampad.eqtl & gtex.eqtl & abs(cor.all) > 0.4)
    tfam$set.tmsFilteredTable(tbl.tms.filtered)
    tf.candidates <- unique(tbl.tms.filtered)$tf
    if(length(tf.candidates) > 0)
        tfam$run.trena(tf.candidates)
    tbl.trena <- tfam$get.trenaTable()
    checkEquals(dim(tbl.trena), c(3,8))

} # test_run.trena
#----------------------------------------------------------------------------------------------------
# the above tests were incremental.  this one runs all steps together
test_run.small <- function()
{
    message(sprintf("--- test_run.small"))

    tfam <- TFAM$new(targetGene, trenaProject, tbl.fimo, tbl.gtex.eqtls, tbl.ampad.eqtls, tbl.oc)
    new.start <- 58385244
    new.end   <- 58385416
    tfam$setStudyRegion(start=new.start, end=new.end)

    study.region <- tfam$getStudyRegion()
    checkTrue(study.region$width.kb > 0.16)  # 2M
    tissues <- tfam$getGTEx.eqtl.tissues()
    tfam$set.current.GTEx.eqtl.tissue(tissues[1])

    tfam$run.tms()
    tfam$add.eqtls.toTmsTable()

    tbl.tms <- tfam$get.tmsTable()
    dim(tbl.tms)
    tbl.tms.filtered <- subset(tbl.tms, ampad.eqtl & gtex.eqtl & abs(cor.all) > 0.4)
    dim(tbl.tms.filtered)
    tfam$set.tmsFilteredTable(tbl.tms.filtered)
    tfam$get.tmsFilteredTable()
    tf.candidates <- unique(tbl.tms.filtered$tf)
    checkTrue(length(tf.candidates) == 3)
    tfam$run.trena(tf.candidates)
    tbl.trena <- tfam$get.trenaTable()
    checkEquals(dim(tbl.trena), c(3, 8))
    checkEquals(tbl.trena$gene, c("SP4", "GABPA", "PLAG1"))
    checkEquals(tbl.trena$tfbs, c(1,1,1))

} # test_run.small
#----------------------------------------------------------------------------------------------------
test_run.large <- function()
{
    message(sprintf("--- test_run.large"))

    tfam <- TFAM$new(targetGene, trenaProject, tbl.fimo, tbl.gtex.eqtls, tbl.ampad.eqtls, tbl.oc)
    study.region <- tfam$getStudyRegion()
    checkTrue(study.region$width.kb > 2000)  # 2M
    tissues <- tfam$getGTEx.eqtl.tissues()
    tfam$set.current.GTEx.eqtl.tissue(tissues[1])

    tfam$run.tms()
    tfam$add.eqtls.toTmsTable()

    tbl.tms <- tfam$get.tmsTable()
    dim(tbl.tms)
    tbl.tms.filtered <- subset(tbl.tms, ampad.eqtl & gtex.eqtl & abs(cor.all) > 0.4)
    dim(tbl.tms.filtered)
    tfam$set.tmsFilteredTable(tbl.tms.filtered)
    tf.candidates <- unique(tbl.tms.filtered$tf)
    checkTrue(length(tf.candidates) > 30)
    tfam$run.trena(tf.candidates)
    tbl.trena <- tfam$get.trenaTable()
    checkEquals(dim(tbl.trena), c(47, 8))
    checkTrue(all(tbl.trena$tfbs > 0))

} # test_run.large
#----------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
