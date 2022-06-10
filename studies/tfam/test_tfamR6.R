library(RUnit)
library(TrenaProjectAD)
source("TFAM-class.R")
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
   tbl.ampad.eqtls <- subset(tbl.ampad.eqtls.raw, study=="ampad-rosmap" & pvalue < 0.05)
   tbl.gtex.eqtls.raw <- get(load("data/tbl.eqtls.6.tissues-chr10:57385407-59385407.2022-06-08.08:03:52"))
   tbl.gtex.eqtls <- subset(tbl.gtex.eqtls.raw, gene==targetGene & pvalue < 0.05)
   data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
   filename <- "mayoAllPeaks.merged.96064x4.RData"
   tbl.mayoAtac <- get(load(file.path(data.dir, filename)))
    #   data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
    # filename <- "boca-hg38-consensus-ATAC.RData"
   tbl.oc <- tbl.mayoAtac
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
    checkEquals(study.region$chrom, "chr10")
    checkEquals(study.region$start, 57391973)
    checkEquals(study.region$end,   59392030)

    new.start <- 58384216
    new.end   <- 58390229
    tfam$setStudyRegion(chrom="chr10", start=new.start, end=new.end)
    new.region <- tfam$getStudyRegion()

    checkEquals(new.region$chrom, "chr10")
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
run.small <- function(shoulder=86)
{
    tfam <- TFAM$new(targetGene, trenaProject, tbl.fimo, tbl.gtex.eqtls, tbl.ampad.eqtls, tbl.oc)
    chrom <- "chr10"
    center <- 58385330
    new.start <- center - shoulder
    new.end   <- center + shoulder
    tfam$setStudyRegion(chrom=chrom, start=new.start, end=new.end)

    study.region <- tfam$getStudyRegion()
    tissues <- tfam$getGTEx.eqtl.tissues()
    tfam$set.current.GTEx.eqtl.tissue(tissues[1])

    tfam$run.tms()
    tfam$add.eqtls.toTmsTable()

    tbl.tms <- tfam$get.tmsTable()
    tbl.tms.filtered <- subset(tbl.tms, ampad.eqtl & gtex.eqtl & abs(cor.all) > 0.4)
    tfam$set.tmsFilteredTable(tbl.tms.filtered)
    tfam$get.tmsFilteredTable()
    tf.candidates <- unique(tbl.tms.filtered$tf)
    tfam$run.trena(tf.candidates)
    tbl.trena <- tfam$get.trenaTable()

    return(tfam)

} # run.small
#----------------------------------------------------------------------------------------------------
# the above tests were incremental.  this one runs all steps together
test_run.small <- function()
{
    message(sprintf("--- test_run.small"))

    tfam.small <- run.small()
    study.region <- tfam.small$getStudyRegion()
    checkTrue(study.region$width.kb < 0.20)

    tbl.tms <- tfam.small$get.tmsFilteredTable()
    checkEquals(dim(tbl.tms), c(3, 18))

    tbl.trena <- tfam.small$get.trenaTable()
    checkEquals(dim(tbl.trena), c(3, 8))
    checkEquals(tbl.trena$gene, c("SP4", "GABPA", "PLAG1"))
    checkEquals(tbl.trena$tfbs, c(1,1,1))

} # test_run.small
#----------------------------------------------------------------------------------------------------
# the above tests were incremental.  this one runs all steps together
test_run.2k.break.motifs <- function()
{
    message(sprintf("--- test_run.2k.break.motifs"))

    tfam <- TFAM$new(targetGene, trenaProject, tbl.fimo, tbl.gtex.eqtls, tbl.ampad.eqtls, tbl.oc)
    new.start <- 58384447
    new.end   <- 58386447
    tfam$setStudyRegion(start=new.start, end=new.end)

    study.region <- tfam$getStudyRegion()
    checkTrue(study.region$width.kb >= 2)
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
    checkTrue(length(tf.candidates) == 5)
    tfam$run.trena(tf.candidates)
    tbl.trena <- tfam$get.trenaTable()
    checkEquals(dim(tbl.trena), c(5, 8))
    checkEquals(tbl.trena$gene, c("SP4", "GABPA", "HES7", "PLAG1", "NFYA"))
    checkEquals(tbl.trena$tfbs, c(2,1,1,1,1))
    tfam$breakMotifs(tbl.trena, tbl.tms)
    tbl.breaks <- tfam$getBreaksTable()

} # test_run.2k.break.motifs
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
    length(tf.candidates)
    checkTrue(length(tf.candidates) > 30)
    tfam$run.trena(tf.candidates)
    tbl.trena <- tfam$get.trenaTable()
    checkEquals(dim(tbl.trena), c(42, 8))
    checkTrue(all(tbl.trena$tfbs > 0))

    tbl.trena.strong <- subset(tbl.trena, abs(betaLasso) > .05 | rfNorm > 0.3)
    rownames(tbl.trena.strong) <- NULL
    tbl.tfbs <- subset(tbl.tms, tf %in% tbl.trena.strong$gene)
    save(tbl.trena.strong, tbl.tfbs, file="full-run-large.RData")

} # test_run.large
#----------------------------------------------------------------------------------------------------
test_viz <- function()
{
   message(sprintf("--- test_viz"))

   ts <- run.small(shoulder=300)

   tbl.trena <- ts$get.trenaTable()
   tbl.tms   <- ts$get.tmsFilteredTable()
   ts$breakMotifs(tbl.trena, tbl.tms)
   tbl.breaks <- ts$get.breaksTable()

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
       browser(); xyz <-99
       }
   tbl.breaks.tfbs <- do.call(rbind, tbls)
   rownames(tbl.breaks.tfbs) <- NULL


   igv <- start.igv(targetGene, "hg38")
   ts$viz(igv)


} # test_viz
#----------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
