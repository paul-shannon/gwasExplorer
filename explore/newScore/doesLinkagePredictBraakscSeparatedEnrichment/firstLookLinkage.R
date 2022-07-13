library(GenomicRanges)
library(RUnit)
library(EndophenotypeExplorer)
library(gwasExplorer)
library(plyr)
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
test.variant.for.braaksc.separation <- function()
{
   targetGene <- "NDUFS2"

      # first, rna-seq separated on variant
   etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD")
   tbl.eqtls.ndufs2 <- etx$get.ampad.EQTLsForGene()
   dim(tbl.eqtls.ndufs2) # 11310    10
   tbl.eqtls.ndufs2 <- subset(tbl.eqtls.ndufs2, study=="ampad-rosmap")
   dim(tbl.eqtls.ndufs2) #  3774    10
   tbl.track <- tbl.eqtls.ndufs2[, c("chrom", "hg19", "hg19", "pvalue", "rsid")]
   colnames(tbl.track) <- c("chrom", "start", "end", "pvalue", "rsid")
   tbl.track$start <- tbl.track$start - 1
   tbl.track$score <- round(-log10(tbl.track$pvalue), digits=3)
   track <- DataFrameQuantitativeTrack("NDUFS2 eQTL", tbl.track[, c("chrom", "start", "end", "score", "rsid")],
                                       color="brown", autoscale=FALSE,
                                       min=0, max=20)
   displayTrack(igv, track)

   rsid.oi <- subset(tbl.track, score == max(tbl.track$score))$rsid
   mtx.geno.1 <- etx$getGenoMatrixByRSID(rsid.oi)
   mtx.geno.pt.rosmap <- etx$subsetAndRelabelGenoMatrixByPatientIDs(mtx.geno.1, "rosmap")

   tbl.pt <- etx$get.rosmap.patient.table(NA)
   rosmap.patients <- intersect(tbl.pt$individualID, colnames(mtx.geno.pt.rosmap))
   length(rosmap.patients)   # 1143

   tbl.pt.rosmap <- subset(tbl.pt, individualID %in% rosmap.patients)
   dim(tbl.pt.rosmap)  # 1143 18

   dim(mtx.geno.pt.rosmap)
   table(mtx.geno.pt.rosmap)  # mtx.geno.pt.rosmap
                              # 0/0 0/1 1/1
                              # 827 296  28
   pt.ad <-  subset(tbl.pt.rosmap, braaksc >=5)$individualID   # 303
   pt.ctl <- subset(tbl.pt.rosmap, braaksc <=1)$individualID  # 85


   table(mtx.geno.pt.rosmap)  # mtx.geno.pt.rosmap
                              #  0/0 0/1 1/1
                              #  827 296  28

   pt.ad <-  subset(tbl.pt.rosmap, braaksc >=5)$individualID   # 303
   pt.ctl <- subset(tbl.pt.rosmap, braaksc <=1)$individualID  # 85

   table(mtx.geno.pt.rosmap[,pt.ad])    # 0/0 0/1 1/1
                                        # 215  83   5
   table(mtx.geno.pt.rosmap[,pt.ctl])   # 0/0 0/1 1/1
                                        #  60  24   1
       #----------------------------------------
       #  the deprecated t-test
       #----------------------------------------

   ad.dist <- c(rep(0, 215), rep(1, 83), rep(2, 5))
   ctl.dist <- c(rep(0,60), rep(1, 24), rep(2, 1))
   t.test(ad.dist, ctl.dist)$p.value    # 0.986

       #----------------------------------------
       #  preferred # 1: chi-square
       #----------------------------------------

   tbl.summary <- data.frame(wt=c(215, 60), het=c(83, 24), hom=c(5, 1),
                             stringsAsFactors=FALES,
                             row.names=c("ad", "ctl"))
   mosaicplot(tbl.summary, color=TRUE, main=rsid.oi)
   chisq.test(tbl.summary)$p.value   # 0.944
   fisher.test(tbl.summary)$p.value  # 0.968

} # test.variant.for.braaksc.separation <- function()
#----------------------------------------------------------------------------------------------------
tag.snp.haplotype.relations <- function()
{
  f <- system.file(package="gwasExplorer", "extdata", "adamts4.locus",
                   "haploreg-rs4575098-0.2-hg19-hg38.RData")
  tbl.haploReg <- get(load(f))
  tbl.track <- tbl.haploReg[, c("chrom", "hg19", "hg19", "rSquared")]
  colnames(tbl.track)[2:3] <- c("start", "end")
  tbl.track$start <- tbl.track$start - 1
  track <- DataFrameQuantitativeTrack("haploreg", tbl.track, color="black", autoscale=FALSE,
                                      min=0, max=1)
  displayTrack(igv, track)

} # tag.snp.haplotype.relations
#----------------------------------------------------------------------------------------------------
viz <- function()
{
   igv <- start.igv("NDUFS2", "hg19")
   showGenomicRegion(igv, "chr1:161,039,989-161,313,300")

} # viz
#----------------------------------------------------------------------------------------------------
rosmap.ad.ctl.separation <- function(rsid)
{
   mtx.geno.1 <- etx$getGenoMatrixByRSID(rsid)
   mtx.geno.pt.rosmap <- etx$subsetAndRelabelGenoMatrixByPatientIDs(mtx.geno.1, "rosmap")

   tbl.pt <- etx$get.rosmap.patient.table(NA)
   rosmap.patients <- intersect(tbl.pt$individualID, colnames(mtx.geno.pt.rosmap))
   length(rosmap.patients)   # 1143

   tbl.pt.rosmap <- subset(tbl.pt, individualID %in% rosmap.patients)
   dim(tbl.pt.rosmap)  # 1143 18

   dim(mtx.geno.pt.rosmap)
   table(mtx.geno.pt.rosmap)  # mtx.geno.pt.rosmap
                              # 0/0 0/1 1/1
                              # 827 296  28
   pt.ad <-  subset(tbl.pt.rosmap, braaksc >=5)$individualID   # 303
   pt.ctl <- subset(tbl.pt.rosmap, braaksc <=1)$individualID  # 85


   table(mtx.geno.pt.rosmap)  # mtx.geno.pt.rosmap
                              #  0/0 0/1 1/1
                              #  827 296  28

   pt.ad <-  subset(tbl.pt.rosmap, braaksc >=5)$individualID   # 303
   pt.ctl <- subset(tbl.pt.rosmap, braaksc <=1)$individualID  # 85

   ad.wt  <- length(grep("0/0", (mtx.geno.pt.rosmap[, pt.ad])))
   ad.het <- length(grep("0/1", (mtx.geno.pt.rosmap[, pt.ad])))
   ad.hom <- length(grep("1/1", (mtx.geno.pt.rosmap[, pt.ad])))

   ctl.wt  <- length(grep("0/0", (mtx.geno.pt.rosmap[, pt.ctl])))
   ctl.het <- length(grep("0/1", (mtx.geno.pt.rosmap[, pt.ctl])))
   ctl.hom <- length(grep("1/1", (mtx.geno.pt.rosmap[, pt.ctl])))

   tbl.summary <- data.frame(wt=c(ad.wt, ctl.wt), het=c(ad.het, ctl.het), hom=c(ad.hom, ctl.hom),
                             row.names=c("ad", "ctl"))
   print(tbl.summary)

   fisher.test(tbl.summary)$p.value

} # rosmap.ad.ctl.separation
#----------------------------------------------------------------------------------------------------
test_rosmap.ad.ctl.separation <- function()
{
   rsid <- "rs1136207"
   x <- rosmap.ad.ctl.separation(rsid)
   checkEqualsNumeric(x, 0.0400, tol=1e-2)

   rsid <- "rs528321"
   x <- rosmap.ad.ctl.separation(rsid)
   checkEqualsNumeric(x, 0.0820, tol=1e-2)

   rsid <- "rs352680"
   x <- rosmap.ad.ctl.separation(rsid)
   checkEqualsNumeric(x, 0.666, tol=1e-2)

} # test_rosap.ad.ctl.separation
#----------------------------------------------------------------------------------------------------
is.ad.ctl.separation.by.extreme.braaksc.predicted.by.haploreg.and.rosmap.eqtl <- function()
{
   etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD")
   tbl.eqtls.ndufs2 <- etx$get.ampad.EQTLsForGene()
   tbl.eqtls.ndufs2 <- subset(tbl.eqtls.ndufs2, study=="ampad-rosmap")
   dim(tbl.eqtls.ndufs2) #  3774    10

   rsids.shared <- intersect(tbl.eqtls.ndufs2$rsid, tbl.haploReg$rsid)
   length(rsids.shared)  # 72
   rsids.rosmap.only <- setdiff(subset(tbl.eqtls.ndufs2, pvalue <= 1e-2)$rsid, tbl.haploReg$rsid)
   length(rsids.rosmap.only) # 61

   tbl.track <- subset(tbl.eqtls.ndufs2, rsid %in% rsids.shared)[, c("chrom", "hg19", "hg19", "rsid")]
   track <- DataFrameAnnotationTrack("shared", tbl.track, color="brown", trackHeight=25)
   displayTrack(igv, track)

   shared.sep.sig.list <- lapply(rsids.shared, rosmap.ad.ctl.separation)
   names(shared.sep.sig.list) <- rsids.shared

   rosmap.only.sep.sig.list <- lapply(rsids.rosmap.only, rosmap.ad.ctl.separation)
   names(rosmap.only.sep.sig.list) <- rsids.rosmap.only

   if(FALSE){
       boxplot(-log10(as.numeric(rosmap.only.sep.sig.list)),
               -log10(as.numeric(shared.sep.sig.list)), main="-log10 braaksc separation",
               names=c("rosmap, no haploreg", "rosmap & haploreg"))
       hist(-log10(as.numeric(rosmap.only.sep.sig.list)), main="rosmap only", xlim=c(0,3))
       hist(-log10(as.numeric(shared.sep.sig.list)), main="rosmap & haploreg")
       } # FALSE: boxplot


   # which(rosmap.only.sep.sig.list == min(as.numeric(rosmap.only.sep.sig.list)))
   # rs6696411  46
   if(FALSE){
      indices <- match(names(shared.sep.sig.list), tbl.eqtls.ndufs2$rsid)
      tbl.track <- tbl.eqtls.ndufs2[indices, c("chrom", "hg19", "hg19")]
      tbl.track$score <- -log10(as.numeric(shared.sep.sig.list))
      colnames(tbl.track)[2:3] <- c("start", "end")
      tbl.track$start <- tbl.track$start - 1
      track <- DataFrameQuantitativeTrack("shared braaksc separation", tbl.track, color="blue",
                                          autoscale=FALSE, min=0, max=3)
      displayTrack(igv, track)

      indices <- match(names(rosmap.only.sep.sig.list), tbl.eqtls.ndufs2$rsid)
      tbl.track <- tbl.eqtls.ndufs2[indices, c("chrom", "hg19", "hg19")]
      tbl.track$score <- -log10(as.numeric(rosmap.only.sep.sig.list))
      colnames(tbl.track)[2:3] <- c("start", "end")
      tbl.track$start <- tbl.track$start - 1
      track <- DataFrameQuantitativeTrack("rosmap only braaksc separation", tbl.track, color="green",
                                          autoscale=FALSE, min=0, max=3)
      displayTrack(igv, track)

      } # FALSE: sep


} # is.ad.ctl.separation.by.extreme.braaksc.predicted.by.haploreg.and.rosmap.eqtl
#----------------------------------------------------------------------------------------------------
# rs429358 rs7412  Name
#          C	T   ε1
#          T	T   ε2
#          T	C   ε3
#          C	C   ε4
apoe.and.braaksc <- function()
{
   rsid <- "rs429358"
   etx <- EndophenotypeExplorer$new("APOE", "hg19", vcf.project="AMPAD")
   x <- rosmap.ad.ctl.separation(rsid) # 0.00298

   rsid <- "rs7412"
   x <- rosmap.ad.ctl.separation(rsid) #  0.000374


} # apoe.and.braaksc
#----------------------------------------------------------------------------------------------------
#   T             C        rs429358 & rs7412
# rs429358      rs7412        Genotype
#
# (C;C)         (T;T)           ε1/ε1
# (C;T)         (T;T)           ε1/ε2
# (C;T)         (C;T)           ε2/ε4
# (C;C)         (C;T)           ε1/ε4
# (T;T)         (T;T)           ε2/ε2
# (T;T)         (C;T)           ε2/ε3
# (T;T)         (C;C)           ε3/ε3
# (C;T)         (C;C)           ε3/ε4
# (C;C)         (C;C)           ε4/ε4
apoe.and.braaksc.allCombinations <- function()
{
   targetGene <- "APOE"
   etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD", initialize.snpLocs=TRUE)
   snp.1 <- "rs429358"  # 19:45411941_T/C
   snp.2 <- "rs7412"    # 19:45412079_C/T
   mtx.geno.1 <- etx$getGenoMatrixByRSID(snp.1)
   mtx.geno.2 <- etx$getGenoMatrixByRSID(snp.2)
   dim(mtx.geno.1)
   dim(mtx.geno.2)
   mtx.geno <- rbind(mtx.geno.1, mtx.geno.2)
   mtx.geno <- etx$subsetAndRelabelGenoMatrixByPatientIDs(mtx.geno, "rosmap")

                                                                                                                  # patients
   e1.e1 <- intersect(names(which(mtx.geno[1,] == "1/1")), names(which(mtx.geno[2,] == "1/1"))); length(e1.e1)    #   0
   e1.e2 <- intersect(names(which(mtx.geno[1,] == "0/1")), names(which(mtx.geno[2,] == "1/1"))); length(e1.e2)    #   0
   e2.e4 <- intersect(names(which(mtx.geno[1,] == "0/1")), names(which(mtx.geno[2,] == "0/1"))); length(e2.e4)    #  21
   e1.e4 <- intersect(names(which(mtx.geno[1,] == "1/1")), names(which(mtx.geno[2,] == "0/1"))); length(e1.e4)    #   0
   e2.e2 <- intersect(names(which(mtx.geno[1,] == "0/0")), names(which(mtx.geno[2,] == "1/1"))); length(e2.e2)    #   7
   e2.e3 <- intersect(names(which(mtx.geno[1,] == "0/0")), names(which(mtx.geno[2,] == "0/1"))); length(e2.e3)    # 142
   e3.e3 <- intersect(names(which(mtx.geno[1,] == "0/0")), names(which(mtx.geno[2,] == "0/0"))); length(e3.e3)    # 699
   e3.e4 <- intersect(names(which(mtx.geno[1,] == "0/1")), names(which(mtx.geno[2,] == "0/0"))); length(e3.e4)    # 256
   e4.e4 <- intersect(names(which(mtx.geno[1,] == "1/1")), names(which(mtx.geno[2,] == "0/0"))); length(e4.e4)    #  20

   tbl.pt.rosmap <- subset(tbl.pt, individualID %in% rosmap.patients)
   dim(tbl.pt.rosmap)  # 1143 18
                                                                                               # fraction with braaksc >= 5
   x <- subset(tbl.pt.rosmap, individualID %in% e2.e4)$braaksc; length(which(x >=5))/length(x) #  0.286
   x <- subset(tbl.pt.rosmap, individualID %in% e2.e2)$braaksc; length(which(x >=5))/length(x) #  0
   x <- subset(tbl.pt.rosmap, individualID %in% e2.e3)$braaksc; length(which(x >=5))/length(x) #  0.063
   x <- subset(tbl.pt.rosmap, individualID %in% e3.e3)$braaksc; length(which(x >=5))/length(x) #  0.233
   x <- subset(tbl.pt.rosmap, individualID %in% e3.e4)$braaksc; length(which(x >=5))/length(x) #  0.443
   x <- subset(tbl.pt.rosmap, individualID %in% e4.e4)$braaksc; length(which(x >=5))/length(x) #  0.600

   x <- subset(tbl.pt.rosmap, individualID %in% e2.e2)$braaksc; printf("hi: %d  lo: %d mid: %d", length(which(x >=5)), length(which(x <= 1)), length(intersect(which(x > 1), which(x <5))))
   x <- subset(tbl.pt.rosmap, individualID %in% e2.e3)$braaksc; printf("hi: %d  lo: %d mid: %d", length(which(x >=5)), length(which(x <= 1)), length(intersect(which(x > 1), which(x <5))))
   x <- subset(tbl.pt.rosmap, individualID %in% e2.e4)$braaksc; printf("hi: %d  lo: %d mid:  %d", length(which(x >=5)), length(which(x <= 1)), length(intersect(which(x > 1), which(x <5))))
   x <- subset(tbl.pt.rosmap, individualID %in% e3.e3)$braaksc; printf("hi: %d  lo: %d mid: %d", length(which(x >=5)), length(which(x <= 1)), length(intersect(which(x > 1), which(x <5))))
   x <- subset(tbl.pt.rosmap, individualID %in% e3.e4)$braaksc; printf("hi: %d  lo: %d mid: %d", length(which(x >=5)), length(which(x <= 1)), length(intersect(which(x > 1), which(x <5))))
   x <- subset(tbl.pt.rosmap, individualID %in% e4.e4)$braaksc; printf("hi: %d  lo: %d mid: %d", length(which(x >=5)), length(which(x <= 1)), length(intersect(which(x > 1), which(x <5))))

   length(e4.e4)  # 20

   tbl.pt <- etx$get.rosmap.patient.table(NA)
   # mtx.geno.pt.rosmap <- etx$subsetAndRelabelGenoMatrixByPatientIDs(mtx.geno, "rosmap")

   rosmap.patients <- intersect(tbl.pt$individualID, colnames(mtx.geno.pt.rosmap))
   length(rosmap.patients)   # 1143

   dim(mtx.geno.pt.rosmap)
   table(mtx.geno.pt.rosmap[1,]) # 0/0 0/1 1/1
                                 # 854  277 20
   table(mtx.geno.pt.rosmap[2,]) # 0/0 0/1 1/1
                                 # 981 163   7

   table(tbl.pt.rosmap$braaksc)  #   0   1   2   3   4   5   6
                                 #  13  72 110 285 360 292  11

   sum(as.integer(table(tbl.pt.rosmap$braaksc))) # [1] 1143

   length(intersect(colnames(mtx.geno.pt.rosmap),  rosmap.patients)) # [1] 1143

   tbl.apoe <- data.frame(braak.lo=c(0,   12,   3,  56,  13,  1),
                          braak.mid=c(7, 121,  12, 480, 129,  7),
                          braak.hi=c(0,    9,   6, 163, 113, 12),
                          row.names=c("e2.e2", "e2.e3", "e2.e4", "e3.e3", "e3.e4", "e4.e4"))
   tbl.apoe$not.hi <- with(tbl.apoe, braak.lo + braak.mid)
   tbl.apoe <- tbl.apoe[, c("braak.lo", "braak.mid", "not.hi", "braak.hi")]
   fisher.test(tbl.apoe[, c(1,2)])$p.value   # 7.49 e-6
   fisher.test(tbl.apoe[, c(2,3)], simulate.p.value=TRUE)$p.value   # 0.999
   fisher.test(tbl.apoe[, c(3,4)], simulate.p.value=TRUE)$p.value   # 4.997e-4
   fisher.test(tbl.apoe[, c(1,4)], simulate.p.value=TRUE)$p.value   # 4.997e-4
   fisher.test(tbl.apoe, simulate.p.value=TRUE)$p.value   # 4.997e-4

} # apoe.and.braaksc.allCombinations
#----------------------------------------------------------------------------------------------------
genotype.row.to.mutation.table.row <- function(named.xtab.list)
{
   name <- names(named.xtab.list[1])
   xtab <- named.xtab.list[[1]]
   tbl.raw <- as.data.frame(xtab, stringsAsFactors=FALSE)
   tbl <- as.data.frame(matrix(tbl.raw$Freq, nrow=1, dimnames=list(name, c(tbl.raw$Var1))))

   if(colnames(tbl)[1] == "./.")
       tbl <- tbl[,-1, drop=FALSE]
   if(! "0/0" %in% colnames(tbl)){
       tbl[1, "0/0"] <- 0
       tbl <- tbl[, sort(colnames(tbl))]
       }

   #if(ncol(tbl) == 1) # just wildtype
   #    pc <- 0
   #if(ncol(tbl) > 1)
   #    pc <- round(100 * sum(tbl[1, -1])/sum(tbl[1,]), digits=0)
   #tbl$pc <- pc

   return(tbl)

} # genotype.row.to.mutation.table.row
#----------------------------------------------------------------------------------------------------
test_genotype.row.to.mutation.table.row <- function()
{
    message(sprintf("--- test_genotype.row.to.mutation.percentage"))
    load("xtabs.RData")
    tbl.ctl <- genotype.row.to.mutation.table.row(list(ctl=xtab.ctl))
    tbl.ad  <- genotype.row.to.mutation.table.row(list(ad=xtab.ad))
    mtx.both <- as.matrix(rbind.fill(tbl.ctl, tbl.ad))
    rownames(mtx.both) <- c("ctl", "ad")
    mtx.both[is.na(mtx.both)] <- 0
    tbl.both <- as.data.frame(mtx.both)
    checkEquals(dim(tbl.both), c(2, 3))
    checkEquals(rownames(tbl.both), c("ctl", "ad"))
    checkEquals(colnames(tbl.both), c("0/0", "0/1", "1/1"))
    checkEquals(as.numeric(tbl.both[1,,drop=TRUE]), c(63, 22, 0))
    checkEquals(as.numeric(tbl.both[2,,drop=TRUE]), c(179, 106, 18))
    pval <- fisher.test(tbl.both, simulate.p.value=TRUE)$p.value   # 0.00799
    checkEqualsNumeric(pval, 0.005, tol=1e-4)

} # test_genotype.row.to.mutation.table.row
#----------------------------------------------------------------------------------------------------
rsid.to.rosmap.braaksc.association <- function(rsid, chromosome)
{
     # targetGene APOE will be ignored
   etx <- EndophenotypeExplorer$new("APOE", "hg19", vcf.project="AMPAD", initialize.snpLocs=TRUE)
   etx$setTargetChromosome(chromosome, "hg19")
   mtx.geno <- etx$getGenoMatrixByRSID(rsid)
   mtx.geno <- etx$subsetAndRelabelGenoMatrixByPatientIDs(mtx.geno, "rosmap")
   dim(mtx.geno)
   tbl.pt <- etx$get.rosmap.patient.table(NA)
   #browser()

   rosmap.patients <- intersect(tbl.pt$individualID, colnames(mtx.geno))
   length(rosmap.patients)   # 1143

   pt.ad <-  subset(tbl.pt, braaksc >=5)$individualID   # 431
   pt.ctl <- subset(tbl.pt, braaksc <=1)$individualID   # 121
   pt.ad <- intersect(pt.ad, colnames(mtx.geno))        # 303
   pt.ctl <- intersect(pt.ctl, colnames(mtx.geno))      # 85

   xtab.ctl <- list(table(mtx.geno[, pt.ctl]))
   names(xtab.ctl) <- sprintf("%s-ctl", rsid)
   xtab.ad  <- list(table(mtx.geno[, pt.ad]))
   names(xtab.ad) <- sprintf("%s-ad", rsid)

   tbl.ctl <- genotype.row.to.mutation.table.row(xtab.ctl)
   tbl.ad  <- genotype.row.to.mutation.table.row(xtab.ad)

   mtx.both <- as.matrix(rbind.fill(tbl.ctl, tbl.ad))
   rownames(mtx.both) <- c(names(xtab.ctl), names(xtab.ad))

   mtx.both[is.na(mtx.both)] <- 0
   tbl.both <- as.data.frame(mtx.both)
   tbl.both

} # rsid.to.rosmap.braaksc.association
#----------------------------------------------------------------------------------------------------
test_rsid.to.rosmap.braaksc.association <- function()
{
   message(sprintf("--- test_rsid.to.rosmap.braaksc.association"))
   rsid.to.rosmap.braaksc.association("rs4575098", "ADAMTS4")



} # test_rsid.to.rosmap.braaksc.association
#----------------------------------------------------------------------------------------------------
run_rsid.to.rosmap.braaksc.association <- function()
{
   data.dir <- "~/github/TrenaProjectAD/inst/extdata/gwasLoci"
   file.1 <- "schwartzentruber-2021-with-hg38-locs.RData"
   full.path.1 <- file.path(data.dir, file.1)
   checkTrue(file.exists(full.path.1))
   tbl.snps.1 <- get(load(full.path.1))
   file.2 <- "posthuma-2019-with-hg38-locs.RData"
   full.path.2 <- file.path(data.dir, file.2)
   checkTrue(file.exists(full.path.2))
   tbl.snps.2 <- get(load(full.path.2))
   tbl.snps.2$source <- "posthuma"
   tbl.snps.1$source <- "schwartzentruber"

   tbl.snps <- rbind(tbl.snps.1, tbl.snps.2)
   dups <- which(duplicated(tbl.snps$rsid))
   length(dups)
   if(length(dups) > 0)
       tbl.snps <- tbl.snps[-dups,]

   tbls <- list()

   for(r in 46:147){
       tryCatch({
          rsid <- tbl.snps$rsid[r]
          chrom <- sprintf("chr%s", tbl.snps$chrom[r])
          source <- tbl.snps$source[r]
          printf("%3d) %s", r, rsid)
          suppressMessages({
             tbl.both <- rsid.to.rosmap.braaksc.association(rsid, chrom)
             sig.pval <-  fisher.test(tbl.both, simulate.p.value=TRUE)$p.value
             })
          tbls[[r]] <- data.frame(rsid=rsid, sig=sig.pval, source=source, stringsAsFactors=FALSE)
          },
       error = function(e){
           printf("failed on rsid %s", rsid)
           })
    } # for r

   tbl.out <- do.call(rbind, tbls)
   new.order <- order(tbl.out$sig, decreasing=FALSE)
   tbl.out <- tbl.out[new.order,]
   rownames(tbl.out) <- NULL
   data.dir <- "~/github/TrenaProjectAD/inst/extdata/gwasLoci"
   filename <- "schwartzentruber.posthuma.braaksc.sepration.RData"
   save(tbl.out, file=file.path(data.dir, filename))


} # run_rsid.to.rosmap.braaksc.association
#----------------------------------------------------------------------------------------------------
picalm.and.apoe4.createFeatureTable <- function()
{
   igv <- start.igv("PICALM", "hg19")
   igv.hg38 <- start.igv("PICALM", "hg38")
   tbl.schwartz <- get(load("~/github/TrenaProjectAD/inst/extdata/gwasLoci/schwartzentruber-2021-with-hg38-locs.RData"))
   tbl.pos <- get(load("~/github/TrenaProjectAD/inst/extdata/gwasLoci/tbl.posthuma-38-geneAssociations-curated-3828x12.RData"))
   locs.hg38 <- getGenomicRegion(igv.hg38)
   locs.hg19 <- getGenomicRegion(igv)
   # dim(with(locs.hg38, subset(tbl.schwartz, hg38 > start & hg38 < end & chrom=="11")))

   etx <- EndophenotypeExplorer$new("PICALM", "hg19", vcf.project="AMPAD", initialize.snpLocs=TRUE)
   mtx.geno.picalm <- with(locs.hg19,  etx$getGenoMatrix("11", start, end))

   apoe.variants <- c("rs429358", "rs7412")
   etx$setTargetChromosome("chr19", "hg19")
   mtx.geno.apoe <- etx$getGenoMatrixByRSID(apoe.variants)

   mtx.geno <- rbind(mtx.geno.picalm, mtx.geno.apoe)
   mtx.geno <- etx$subsetAndRelabelGenoMatrixByPatientIDs(mtx.geno, "rosmap")

   tbl.pt <- etx$get.rosmap.patient.table(NA)
   dim(tbl.pt)
   rosmap.patients <- intersect(tbl.pt$individualID, colnames(mtx.geno))
   length(rosmap.patients)   # 1143
   tbl.pt.rosmap <- subset(tbl.pt, individualID %in% rosmap.patients)
   dim(tbl.pt.rosmap)
   table(tbl.pt.rosmap$braaksc)   #   0   1   2   3   4   5   6
                                  #  13  72 110 285 360 292  11
   table(tbl.pt.rosmap$ceradsc)   #   1   2   3   4
                                  # 373 397 106 267
   names(table(mtx.geno))
   #  "./." "0/0" "0/1" "0/2" "0/3" "0/4" "0/5" "0/6" "1/1" "1/2" "1/3" "1/4" "1/5" "1/6"
   #  "2/2" "2/3" "2/4" "2/5" "2/6" "3/3" "3/4" "3/5" "3/6" "4/4" "4/5" "4/6" "5/5" "5/6" "6/6"
   mtx.geno.fixed <- mtx.geno
   mtx.geno.fixed[mtx.geno.fixed == "./."] <- "0"
   mtx.geno.fixed[mtx.geno.fixed == "0/0"] <- "0"
   mtx.geno.fixed[mtx.geno.fixed == "0/1"] <- "1"
   mtx.geno.fixed[mtx.geno.fixed == "0/2"] <- "1"
   mtx.geno.fixed[mtx.geno.fixed == "0/3"] <- "1"
   mtx.geno.fixed[mtx.geno.fixed == "0/4"] <- "1"
   mtx.geno.fixed[mtx.geno.fixed == "0/5"] <- "1"
   mtx.geno.fixed[mtx.geno.fixed == "0/6"] <- "1"
   mtx.geno.fixed[mtx.geno.fixed == "1/1"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "1/2"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "1/3"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "1/4"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "1/5"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "1/6"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "2/2"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "2/3"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "2/4"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "2/5"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "2/6"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "3/3"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "3/4"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "3/5"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "3/6"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "4/4"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "4/5"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "4/6"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "5/5"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "5/6"] <- "2"
   mtx.geno.fixed[mtx.geno.fixed == "6/6"] <- "2"

   table(mtx.geno.fixed)  #       0       1       2
                          # 3340634  183989  126349
   mtx.geno.fixed.2 <- matrix(data=as.numeric(mtx.geno.fixed), nrow=nrow(mtx.geno.fixed))

   colnames(mtx.geno.fixed.2) <- colnames(mtx.geno)
   rownames(mtx.geno.fixed.2) <- rownames(mtx.geno)
   all(is.numeric(colSums(mtx.geno.fixed.2)))

   library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
   locs.raw <- sub("_.*$", "", rownames(mtx.geno.fixed.2))
   gr.locs <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP144.GRCh37, GRanges(locs.raw))
   length(locs.raw)   # 3172
   length(gr.locs)    # 1794
     # todo: create a complete map of locs.raw to rsids, accounting for 1378 failures
   rsids <- gr.locs$RefSNP_id

   tokens <- strsplit(locs.raw, ":")
   chroms <- unlist(lapply(tokens, "[", 1))
   loc <- as.numeric(unlist(lapply(tokens, "[", 2)))

   tbl.features <- as.data.frame(t(mtx.geno.fixed.2))
   dups <- setdiff(rownames(tbl.features), rosmap.patients)
   length(dups)
   if(length(dups) > 0){
      dup.indices <- match(dups, rownames(tbl.features))
      tbl.features <- tbl.features[-dup.indices,]
      }
   dim(tbl.features)
   stopifnot(all(unlist(lapply(tbl.features, class) == "numeric")))

   pt.indices <- match(rownames(tbl.features), tbl.pt$individualID)
   tbl.features$braak <- as.numeric(tbl.pt$braaksc[pt.indices])
   tbl.features$cerad <- as.numeric(tbl.pt$ceradsc[pt.indices])
   tbl.features$cogdx <- as.numeric(tbl.pt$cogdx[pt.indices])
   stopifnot(all(unlist(lapply(tbl.features, class) == "numeric")))

   dim(tbl.features)
   tbl.features[1:5, c(1:5, 3170:3175)]
   save(tbl.features, file="tbl.features.1143ptx3175variants.RData")

} # picalm.and.apoe4.createFeatureTable
#----------------------------------------------------------------------------------------------------
picalm.and.apoe4.analyze <- function()
{
   dim(tbl.features)

  tbl.features.sub <- subset(tbl.features, braak == 1 | braak >=5)
  dim(tbl.features.sub)
  tbl.features.sub$braak2 <- 0
  tbl.features.sub$braak2[tbl.features.sub$braak >=5] <- 1
  candidates <- c("19:45411941_T/C", "11:85693353_C/T", "11:85663499_T/G", "11:85749737_G/A")
  tbl.features.sub <- tbl.features.sub[, c(candidates, "braak2")]
  table(tbl.features.sub$braak2)
  model <- lm(braak2 ~ 0 + ., data=tbl.features.sub)
  mtx.coef <- summary(model)$coefficients
  colnames(mtx.coef)[4] <- "pval"

  new.order <- order(mtx.coef[, "pval"], decreasing=FALSE)
  tbl.coef <- as.data.frame(mtx.coef[new.order,])
  tbl.coef <- subset(tbl.coef, pval < 0.01)


  model <- lm(braak ~ ., data=tbl.features[, c(1:3173)])
  mtx.coef <- summary(model)$coefficients
  colnames(mtx.coef)[4] <- "pval"
  new.order <- order(mtx.coef[, "pval"], decreasing=FALSE)
  tbl.coef <- as.data.frame(mtx.coef[new.order,])
  tbl.coef <- subset(tbl.coef, pval < 0.01)

  candidates <- rownames(subset(tbl.coef, pval < 0.01))
  candidates <- gsub("`", "", candidates)
  lapply(candidates, function(candidate) table(tbl.features[,candidate]))
  table(tbl.features[,candidates[1]])

  features.hclust <- hclust(dist(tbl.features[, 3173]))
  plot(features.hclust, main="ad")

  cuts <- cutree(features.hclust, 3)
  tree.count <- length(table(cuts))
  print(tree.count)

  for(cluster in 1:tree.count){
     cluster.patients <- names(cuts[cuts==cluster])
     printf("cluster %2d: %d patients", cluster, length(cluster.patients))
     model <- lm(braak2 ~ ., data=tbl.features.sub[cluster.patients, c(1:3173, 3176)])
     model.summary <- summary(model)
     names(model.summary)
     mtx.coef <- model.summary$coefficients
     skip <- any(!is.nan(mtx.coef[, "Pr(>|t|)"]))
     if(!skip){
        sig.snps <- as.integer(which(mtx.coef[, "Pr(>|t|)"] < 0.05))
        printf("sig.snps count: %d", length(sig.snps))
        sig.count <- length(sig.snps)
        if(sig.count > 0){
           mtx.coef.strong <- mtx.coef[sig.snps,]
           print(mtx.coef.strong)
           } # > 0
      } # any
  } # for cluster



} # picalm.and.apoe4.analyze
#----------------------------------------------------------------------------------------------------
