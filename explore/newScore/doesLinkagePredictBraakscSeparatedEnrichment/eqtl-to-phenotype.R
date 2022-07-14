library(EndophenotypeExplorer)
library(RUnit)
targetGene <- "NDUFS2"
if(!exists("etx")){
   etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD")
   tbl.eqtls <- etx$get.ampad.EQTLsForGene()
   dim(tbl.eqtls)
   head(tbl.eqtls)
   table(tbl.eqtls$study)   #   ampad-mayo ampad-rosmap         GTEx
                            #        7532         3774            4
   }
if(!exists("tbl.gwascat.ad")){
   data.dir <- data.dir <- "~/github/TrenaProjectAD/inst/extdata/gwasLoci"
   file.exists(data.dir)
   file <- "alzheimerSubsetOfGWASCatatalog-29apr2022.RData"
   full.path <- file.path(data.dir, file)
   file.exists(full.path)
   tbl.gwascat.ad <- get(load(full.path))
   stopifnot(colnames(tbl.gwascat.ad)[c(1,2,33)] == c("seqnames", "start", "P.VALUE"))
   }
#----------------------------------------------------------------------------------------------------
rosmap.ad.ctl.separation <- function(rsid)
{
   printf("----- rosmap.ad.ctl.separation: %s", rsid)

   result <- list(pval.t=NA,
                  pval.fisher=NA,
                  tbl.geno=data.frame(),
                  tbl.pt=data.frame(),
                  mtx.geno=matrix(),
                  mtx.geno.study=matrix(),
                  pt.ad=c(),
                  pt.ctl=c())

   mtx.geno.1 <- etx$getGenoMatrixByRSID(rsid)
   if(all(is.na(mtx.geno.1)))
       return(result)

   mtx.geno.pt.rosmap <- etx$subsetAndRelabelGenoMatrixByPatientIDs(mtx.geno.1, "rosmap")
   if(all(is.na(mtx.geno.pt.rosmap)))
       return(result)

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

   ad.wt <- 0
   at.het <- 0
   ad.hom <- 0

   ad.wt  <- length(grep("0/0", (mtx.geno.pt.rosmap[, pt.ad])))
   ad.het <- length(grep("0/1", (mtx.geno.pt.rosmap[, pt.ad])))
   ad.hom <- length(grep("1/1", (mtx.geno.pt.rosmap[, pt.ad])))

   ctl.wt  <- length(grep("0/0", (mtx.geno.pt.rosmap[, pt.ctl])))
   ctl.het <- length(grep("0/1", (mtx.geno.pt.rosmap[, pt.ctl])))
   ctl.hom <- length(grep("1/1", (mtx.geno.pt.rosmap[, pt.ctl])))

   tbl.summary <- data.frame(wt=c(ad.wt, ctl.wt), het=c(ad.het, ctl.het), hom=c(ad.hom, ctl.hom),
                             row.names=c("ad", "ctl"))
   print(tbl.summary)

   ad.vector <- with(tbl.summary["ad",],
                     c(rep(0, wt),
                       rep(1, het),
                       rep(2, hom)))
   ctl.vector <- with(tbl.summary["ctl",],
                     c(rep(0, wt),
                       rep(1, het),
                       rep(2, hom)))

   #pval.t <- t.test(ad.vector, ctl.vector)$p.value
   #pval.fisher <- fisher.test(tbl.summary)$p.value

   pval.t <- tryCatch({
       t.test(ad.vector, ctl.vector)$p.value
       }, error=function(e){return(1)})

   pval.fisher <- tryCatch({
      fisher.test(tbl.summary)$p.value
      }, error=function(e){return(1)})



   return(list(pval.t=pval.t,
               pval.fisher=pval.fisher,
               tbl.geno=tbl.summary,
               tbl.pt=tbl.pt.rosmap,
               mtx.geno=mtx.geno.1,
               mtx.geno.study=mtx.geno.pt.rosmap,
               pt.ad=pt.ad,
               pt.ctl=pt.ctl))


} # rosmap.ad.ctl.separation
#----------------------------------------------------------------------------------------------------
test_rosmap.ad.ctl.separation <- function()
{
   targetGene <- "MS4A4A"
   if(!exists("etx"))
       etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD")

   RSID <- "rs9916042"
   x <- rosmap.ad.ctl.separation(RSID)
   checkEquals((names(x)), c("pval.t", "pval.fisher", "tbl.geno", "tbl.pt", "mtx.geno", "mtx.geno.study", "pt.ad", "pt.ctl"))

   RSID <- "rs1582763"  # very strong in mayo, almost nothing in rosmap
   checkEquals(sort(names(x)), c("mtx.geno", "mtx.geno.study", "pt.ad", "pt.ctl", "pval",
                                 "tbl.geno", "tbl.pt"))
   #checkEqualsNumeric(x$pval, 0.82851, tol=1e-5)
   checkEqualsNumeric(x$pval, 0.54275, tol=1e-5)
   checkEquals(as.integer(rowSums(x$tbl.geno)), c(303, 85))

   rsid <- "rs199636781"
   x <- rosmap.ad.ctl.separation(rsid)
   checkTrue(is.na(x))

   rsid <- "rs4575098"
   x <- rosmap.ad.ctl.separation(rsid)
   checkEqualsNumeric(x$pval, 0.0058, tol=2e5)

   rsid <- "rs1136207"
   x <- rosmap.ad.ctl.separation(rsid)
   checkEqualsNumeric(x$pval, 0.0400, tol=1e-2)

   rsid <- "rs528321"
   x <- rosmap.ad.ctl.separation(rsid)
   checkEqualsNumeric(x$pval, 0.0820, tol=1e-2)

   rsid <- "rs352680"
   x <- rosmap.ad.ctl.separation(rsid)
   checkEqualsNumeric(x$pval, 0.666, tol=1e-2)

} # test_rosmap.ad.ctl.separation
#----------------------------------------------------------------------------------------------------
mayo.ad.ctl.separation <- function(rsid)
{
   printf("----- mayo.ad.ctl.separation: %s", rsid)
   mtx.geno.1 <- etx$getGenoMatrixByRSID(rsid)
   result <- list(pval.t=NA,
                  pval.fisher=NA,
                  tbl.geno=data.frame(),
                  tbl.pt=data.frame(),
                  mtx.geno=matrix(),
                  mtx.geno.study=matrix(),
                  pt.ad=c(),
                  pt.ctl=c())
   if(all(is.na(mtx.geno.1)))
       return(result)
   mtx.geno.pt.mayo <- etx$subsetAndRelabelGenoMatrixByPatientIDs(mtx.geno.1, "mayo")
   if(all(is.na(mtx.geno.pt.mayo)))
       return(result)
   dim(mtx.geno.pt.mayo)   # 1 349

   tbl.pt <- etx$get.mayo.patient.table(NA)
   mayo.patients <- intersect(tbl.pt$individualID, colnames(mtx.geno.pt.mayo))
   printf("mayo patients with genomes: %d", length(mayo.patients))   # 349

   tbl.pt.mayo <- subset(tbl.pt, individualID %in% mayo.patients)
   dim(tbl.pt.mayo)  # 349 19

   table(mtx.geno.pt.mayo)  # mtx.geno.pt.mayo
                            #  0/0 0/1 1/1
                            #  187 152  10
   pt.ad <-  subset(tbl.pt.mayo, Braak >=5)$individualID   # 86
   pt.ctl <- subset(tbl.pt.mayo, Braak <=1)$individualID   # 32

   pt.ad <- intersect(pt.ad, colnames(mtx.geno.pt.mayo))   # 86
   pt.ctl <- intersect(pt.ctl, colnames(mtx.geno.pt.mayo)) # 32

   ad.wt  <- length(grep("0/0", (mtx.geno.pt.mayo[, pt.ad])))
   ad.het <- length(grep("0/1", (mtx.geno.pt.mayo[, pt.ad])))
   ad.hom <- length(grep("1/1", (mtx.geno.pt.mayo[, pt.ad])))

   ctl.wt  <- length(grep("0/0", (mtx.geno.pt.mayo[, pt.ctl])))
   ctl.het <- length(grep("0/1", (mtx.geno.pt.mayo[, pt.ctl])))
   ctl.hom <- length(grep("1/1", (mtx.geno.pt.mayo[, pt.ctl])))

   tbl.summary <- data.frame(wt=c(ad.wt, ctl.wt), het=c(ad.het, ctl.het), hom=c(ad.hom, ctl.hom),
                             row.names=c("ad", "ctl"))
   print(tbl.summary)

   ad.vector <- with(tbl.summary["ad",],
                     c(rep(0, wt),
                       rep(1, het),
                       rep(2, hom)))
   ctl.vector <- with(tbl.summary["ctl",],
                     c(rep(0, wt),
                       rep(1, het),
                       rep(2, hom)))

   pval.t <- tryCatch({
       t.test(ad.vector, ctl.vector)$p.value
       }, error=function(e){return(1)})

   pval.fisher <- tryCatch({
      fisher.test(tbl.summary)$p.value
      }, error=function(e){return(1)})

   return(list(pval.t=pval.t,
               pval.fisher=pval.fisher,
               tbl.geno=tbl.summary,
               tbl.pt=tbl.pt.mayo,
               mtx.geno=mtx.geno.1,
               mtx.geno.study=mtx.geno.pt.mayo,
               pt.ad=pt.ad,
               pt.ctl=pt.ctl))

} # mayo.ad.ctl.separation
#----------------------------------------------------------------------------------------------------
test_mayo.ad.ctl.separation <- function()
{
   message(sprintf("--- test_mayo.ad.ctl.separation"))

   #targetGene <- "COX7C"
   #etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD")
   #x <- mayo.ad.ctl.separation("rs199636781")
   #tbl.ld <- get(load("~/github/gwasExplorer/studies/locus.cox7c/shared/tbl.linkage.hg19.hg38.RData"))
   #rsids <- tbl.ld$rsid

   targetGene <- "MS4A4A"
   RSID <- "rs1582763"  # very strong in mayo, almost nothing in rosmap
   etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD")
   x <- mayo.ad.ctl.separation(RSID)
   checkEquals(sort(names(x)), c("mtx.geno", "mtx.geno.study", "pt.ad", "pt.ctl", "pval",
                                 "tbl.geno", "tbl.pt"))
   #checkEqualsNumeric(x$pval, 0.000627002, tol=1e-5)
   checkEqualsNumeric(x$pval, 0.69156, tol=1e-5)
     # contrast with rosmap at this locus:  c(303, 85))

   checkEquals(as.integer(rowSums(x$tbl.geno)), c(86, 32))


} # test_mayo.ad.ctl.separation
#----------------------------------------------------------------------------------------------------
test_apoe_variants <- function()
{
   message(sprintf("--- test_apoe_variants"))
   rsids <- c("rs429358", "rs7412")
   targetGene <- "APOE"
   etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD")

   lapply(rsids, mayo.ad.ctl.separation)
     # rs429358  0.0003012019   braak > 5 <= 1
     # rs7412    0.02480738
     # rs429358  0.004436996    braak >= 5 <= 1
     # rs7412    0.08458145


} # test_apoe_variants
#----------------------------------------------------------------------------------------------------
run.gwas.loci <- function()
{
   coi <- c("chrom", "hg38", "rsid", "source", "pvalue")
   data.dir <- "~/github/TrenaProjectAD/inst/extdata/gwasLoci"
   file <- "bellenguez-2022-83variants-hg38.RData"
   full.path <- file.path(data.dir, file)

   tbl.gwas.1 <- get(load(full.path))
   tbl.gwas.1$source <- "bellenguez"
   colnames(tbl.gwas.1)[grep("stage1.pval", colnames(tbl.gwas.1))] <- "pvalue"
   stopifnot(all(coi %in% colnames(tbl.gwas.1)))
   dim(tbl.gwas.1)
   head(tbl.gwas.1)


   file2 <- "schwartzentruber-2021-with-hg38-locs.RData"
   file3 <- "posthuma-2019-with-hg38-locs.RData"

   full.path <- file.path(data.dir, file2)
   tbl.gwas.2 <- get(load(full.path))
   tbl.gwas.2$source <- "schwartzentruber"
   dim(tbl.gwas.2)
   head(tbl.gwas.2)
   stopifnot(all(coi %in% colnames(tbl.gwas.2)))

   full.path <- file.path(data.dir, file3)
   tbl.gwas.3 <- get(load(full.path))
   dim(tbl.gwas.3)
   head(tbl.gwas.3)
   tbl.gwas.3$source <- "posthuma"
   stopifnot(all(coi %in% colnames(tbl.gwas.3)))

   tbl.gwas <- rbind(tbl.gwas.1[, coi],
                     tbl.gwas.2[, coi],
                     tbl.gwas.3[, coi])

   tbl.gwas$chrom <- sub("chr", "", tbl.gwas$chrom)
   tbl.gwas$chrom <- paste0("chr", tbl.gwas$chrom)
   tbl.gwas <- unique(tbl.gwas)
   dim(tbl.gwas)  # 235 5
   save(tbl.gwas, file="~/github/TrenaProjectAD/inst/extdata/gwasLoci/tbl.gwas.3studies.235x5.RData")

   etx <- EndophenotypeExplorer$new("NDUFS2", "hg19", vcf.project="AMPAD")

   tbls <- list()

   for(r in seq_len(nrow(tbl.gwas))){
       chrom <- tbl.gwas$chrom[r]
       rsid <- tbl.gwas$rsid[r]
       etx$setTargetChromosome(chrom, "hg19")
       printf("------ %s (%s)", rsid, targetGene)
       mayo.score <- mayo.ad.ctl.separation(rsid)
       rosmap.score <- rosmap.ad.ctl.separation(rsid)
       printf("        mayo: %10.5f", mayo.score$pval.t)
       printf("      rosmap: %10.5f", rosmap.score$pval.t)
       tbls[[rsid]] <- list(mayo=mayo.score, rosmap=rosmap.score)
                         #data.frame(rsid=rsid, mayo=mayo.score$pval, rosmap=rosmap.score$pval, stringsAsFactors=FALSE)
       save(tbls, file="gwas.bellenquez.braak-truber-posthuma-in-process-3.RData")
       } # for r


   tbl.all <- do.call(rbind, tbls)
   dim(tbl.all)
   save(tbl.all, file="gwas.bellenquez-truber-posthuma.braak-3.RData")

   if(FALSE){  # visualize the outliers, mayo, rosmap, then both
       if(!exists("igv")) igv <- start.igv("all")
       showGenomicRegion(igv, "all")
       tbl.track <- GWASTrack("gwascat", tbl.gwascat.ad, chrom.col=1, pos.col=2, pval.col=33)
       displayTrack(igv, tbl.track)
       rosmap.outliers <- rownames(subset(tbl.rsids, rosmap.f > 2))
       tbl.track <- subset(tbl.gwas, rsid %in% rosmap.outliers)
       track <- GWASTrack("rosmap outliers", tbl.track,  chrom.col=1, pos.col=2, pval.col=5)
       displayTrack(igv, track)
       }

} # run.gwas.loci
#----------------------------------------------------------------------------------------------------
tbls.to.summary <- function()
{
    if(!exists("tbls"))
        load("gwas.bellenquez.braak-truber-posthuma-in-process-3.RData")

    rsids <- names(tbls)
    tbl.list <- list()
    for(rsid in rsids){
       x.mayo <- tbls[[rsid]]$mayo
       x.rosmap <- tbls[[rsid]]$rosmap
       mayo.geno.ad <- as.integer(x.mayo$tbl.geno[1,])
       mayo.geno.ctl <- as.integer(x.mayo$tbl.geno[2,])
       rosmap.geno.ad <- as.integer(x.rosmap$tbl.geno[1,])
       rosmap.geno.ctl <- as.integer(x.rosmap$tbl.geno[2,])
       tbl.rsid <- data.frame(rsid=rsid,
                              mayo.t=-log10(x.mayo$pval.t),
                              mayo.f=-log10(x.mayo$pval.fisher),
                              rosmap.t=-log10(x.rosmap$pval.t),
                              rosmap.f=-log10(x.rosmap$pval.fisher),
                              mayo.geno.ctl=sprintf("%d,%d,%d",
                                                    mayo.geno.ctl[1],
                                                    mayo.geno.ctl[2],
                                                    mayo.geno.ctl[3]),
                              mayo.geno.ad=sprintf("%d,%d,%d",
                                                    mayo.geno.ad[1],
                                                    mayo.geno.ad[2],
                                                    mayo.geno.ad[3]),
                              rosmap.geno.ctl=sprintf("%d,%d,%d",
                                                    rosmap.geno.ctl[1],
                                                    rosmap.geno.ctl[2],
                                                    rosmap.geno.ctl[3]),
                              rosmap.geno.ad=sprintf("%d,%d,%d",
                                                    rosmap.geno.ad[1],
                                                    rosmap.geno.ad[2],
                                                    rosmap.geno.ad[3]),

                              #rosmap.geno=sprintf("ctl: %d,%d,%d  ad: %d,%d,%d",
                              #                  rosmap.geno.ctl[1],
                              #                  rosmap.geno.ctl[2],
                              #                  rosmap.geno.ctl[3],
                              #                  rosmap.geno.ad[1],
                              #                  rosmap.geno.ad[2],
                              #                  rosmap.geno.ad[3]),

                              #mayo.ctl.wt=mayo.geno.ctl[1],
                              #mayo.ctl.het=mayo.geno.ctl[2],
                              #mayo.ctl.hom=mayo.geno.ctl[3],
                              #mayo.ad.wt=mayo.geno.ad[1],
                              #mayo.ad.het=mayo.geno.ad[2],
                              #mayo.ad.hom=mayo.geno.ad[3],

                              #rosmap.ctl.wt=rosmap.geno.ctl[1],
                              #rosmap.ctl.het=rosmap.geno.ctl[2],
                              #rosmap.ctl.hom=rosmap.geno.ctl[3],
                              #rosmap.ad.wt=rosmap.geno.ad[1],
                              #rosmap.ad.het=rosmap.geno.ad[2],
                              #rosmap.ad.hom=rosmap.geno.ad[3],

                              stringsAsFactors=FALSE)
       tbl.list[[rsid]] <- tbl.rsid
       } # for rsid

    tbl.rsids <- do.call(rbind, tbl.list)
    tbl.rsids <- tbl.rsids[, -1]
    filename <- sprintf("tbl.rsids.summary.%s-hetAndHomSeparate.RData",
                        format(Sys.time(), "%a.%b.%d.%Y-%H:%M:%S"))
    save(tbl.rsids, file=filename)

    with(tbl.rsids, plot(y=rosmap.f, x=mayo.f, col="blue", pch=16, cex=1.2,
                         main="-log10(fisher's pval), association\n of genotype and extreme braak score")
         )


    max.values <- as.numeric(apply(as.matrix(tbl.rsids[, c(2,4)]), 1, max))
    tbl.rsids$max <- max.values
    coi <- c("mayo.t", "mayo.f", "rosmap.t", "rosmap.f", "max", "mayo.geno.ctl",
             "mayo.geno.ad", "rosmap.geno.ctl", "rosmap.geno.ad")

    tbl.rsids <- tbl.rsids[, coi]
    tbl.strong <- subset(tbl.rsids, max > 2)
    new.order <- order(tbl.strong$max, decreasing=TRUE)
    tbl.strong <- tbl.strong[new.order,]

    tbl.strong$mayo.t <- round(tbl.strong$mayo.t, digits=2)
    tbl.strong$mayo.f <- round(tbl.strong$mayo.f, digits=2)
    tbl.strong$rosmap.t <- round(tbl.strong$rosmap.t, digits=2)
    tbl.strong$rosmap.f <- round(tbl.strong$rosmap.f, digits=2)
    tbl.strong$max <- round(tbl.strong$max, digits=2)

    tbl.rsids

} # tbls.to.summary
#----------------------------------------------------------------------------------------------------
test_rs34173062 <- function()
{
   rsid.oi <- "rs34173062"

   tbl.rsids <- tbls.to.summary()
   tbl.rsids[rsid.oi,]
      #    mayo.t   mayo.f  rosmap.t rosmap.f    max
      #      3.37     1.89      0.66     0.29   3.37

  tbl.ids <- etx$getIdMap()

  xm <- tbls[[rsid.oi]]$mayo
  xr <- tbls[[rsid.oi]]$rosmap

  tbl.pt.mayo <- xm$tbl.pt #subset(xm$tbl.pt, individualIdSource=="MayoBrainBank")
  dim(tbl.pt.mayo) # 349 19
  pt.mayo <- tbl.pt.mayo$individualID
  pt.mayo.geno <- intersect(pt.mayo, colnames(xm$mtx.geno))
  length(pt.mayo.geno)
  mtx.geno.m <- xm$mtx.geno[, pt.mayo.geno]
  table(mtx.geno.m)

  tbl.pt.rosmap <- xr$tbl.pt
  dim(tbl.pt.rosmap)  # 1143 18
  #pt.rosmap <-
  mtx.geno.r <- xr$mtx.geno
  length(subset(tbl.ids, study=="rosmap")$sample)

  pt.rosmap.geno <- intersect(subset(tbl.ids, study=="rosmap")$sample, colnames(xr$mtx.geno))
  length(pt.rosmap.geno)  # 1151
  mtx.geno.r <- xr$mtx.geno[, pt.rosmap.geno]
  table(mtx.geno.r)


} # test_rs34173062
#----------------------------------------------------------------------------------------------------
test_rs6448451 <- function()
{
  rsid.oi <- "rs6448451"

  tbl.rsids[rsid.oi,]
      #              mayo.t    mayo.f rosmap.t rosmap.f      max
      # rs6448451 0.2445721 0.1163162 3.816741 3.097074 3.816741

  tbl.ids <- etx$getIdMap()

  xm <- tbls[[rsid.oi]]$mayo
  xr <- tbls[[rsid.oi]]$rosmap

  tbl.pt.mayo <- xm$tbl.pt #subset(xm$tbl.pt, individualIdSource=="MayoBrainBank")
  dim(tbl.pt.mayo) # 349 19
  pt.mayo <- tbl.pt.mayo$individualID
  pt.mayo.geno <- intersect(pt.mayo, colnames(xm$mtx.geno))
  length(pt.mayo.geno)
  mtx.geno.m <- xm$mtx.geno[, pt.mayo.geno]
  xm$tbl.geno
  table(mtx.geno.m)


  tbl.pt.rosmap <- xr$tbl.pt
  dim(tbl.pt.rosmap)  # 1143 18
  mtx.geno.r <- xr$mtx.geno
  length(subset(tbl.ids, study=="rosmap")$sample)

  pt.rosmap.geno <- intersect(subset(tbl.ids, study=="rosmap")$sample, colnames(xr$mtx.geno))
  length(pt.rosmap.geno)  # 1151
  mtx.geno.r <- xr$mtx.geno[, pt.rosmap.geno]
  xm$tbl.geno
  table(mtx.geno.r)

} # test_rs6448451
#----------------------------------------------------------------------------------------------------
test_one <- function()
{
  subset(tbl.rsids, rosmap.f > 2 & mayo.f < 0.3)
       #              mayo.t    mayo.f   rosmap.t rosmap.f      max
       # rs429358  3.1721201 2.3529109 4.81230347 3.525583 4.812303
       # rs7412    0.9000917 1.0727249 2.36826934 3.427160 3.427160
       # rs6448451 0.6061332 0.1163162 0.85922285 3.097074 3.097074
       # rs6448453 0.6061332 0.1163162 0.85922285 3.097074 3.097074
       # rs9381563 0.1464635 0.1797320 0.06475329 3.642303 3.642303

  rsid.oi <- "rs9381563"
  rsid.oi <- "rs6448453"
  rsid.oi <- "rs7767350"

  tbl.rsids[rsid.oi,]
      #              mayo.t    mayo.f rosmap.t rosmap.f      max
      # rs6448451 0.2445721 0.1163162 3.816741 3.097074 3.816741

  tbl.ids <- etx$getIdMap()

  xm <- tbls[[rsid.oi]]$mayo
  xr <- tbls[[rsid.oi]]$rosmap

  tbl.pt.mayo <- xm$tbl.pt #subset(xm$tbl.pt, individualIdSource=="MayoBrainBank")
  dim(tbl.pt.mayo) # 349 19
  pt.mayo <- tbl.pt.mayo$individualID
  pt.mayo.geno <- intersect(pt.mayo, colnames(xm$mtx.geno))
  length(pt.mayo.geno)
  mtx.geno.m <- xm$mtx.geno[, pt.mayo.geno]
  xm$tbl.geno
  table(mtx.geno.m)


  tbl.pt.rosmap <- xr$tbl.pt
  dim(tbl.pt.rosmap)  # 1143 18
  mtx.geno.r <- xr$mtx.geno
  length(subset(tbl.ids, study=="rosmap")$sample)

  pt.rosmap.geno <- intersect(subset(tbl.ids, study=="rosmap")$sample, colnames(xr$mtx.geno))
  length(pt.rosmap.geno)  # 1151
  mtx.geno.r <- xr$mtx.geno[, pt.rosmap.geno]
  xr$tbl.geno
  table(mtx.geno.r)
  xm$tbl.geno
  xr$tbl.geno

} # test_rs6448451
#----------------------------------------------------------------------------------------------------
explore.L.shaped.distribution <- function()
{
   RSID <- "rs1582763"  # very strong in mayo, almost nothing in rosmap
   subset(tbl.all, rsid==RSID)
      #                rsid        mayo    rosmap
      # rs1582763 rs1582763 0.000627002 0.8285095

   lapply(list(tbl.gwas.1, tbl.gwas.2, tbl.gwas.3), function(tbl) subset(tbl, rsid==RSID))
      #      rsid chrom     hg38 nearest.gene locus  alleles   maf stage1.or stage1.pval stage2.or stage2.pval stage12.or stage12.pval stage12.hetero stage12.hetero.pval
      #  rs1582763 chr11 60254475       MS4A4A  MS4A     A/G 0.371      0.92    1.65e-24      0.89    2.87e-20       0.91     3.74e-42           33.5               0.101
      #      rsid chrom     hg38 allele pvalue
      # rs1582763    11 60254475      R  5e-09
      # dbSNP says G>A



    # variant <- "19:11178759_GTATATATATATATATATA/G" # or in this case, the row index 81 would also work
    # vcf <- readVcf(url, "hg19", roi)
    # alt(vcf)[variant]

   etx$setTargetChromosome("chr11", "hg19")
   url <- etx$getVcfUrl()
   mtx.geno <- etx$getGenoMatrixByRSID(RSID)
   dim(mtx.geno)
   rownames(mtx.geno)

   tbl.alfa <- etx$getAggregatedAlleleFrequencies(RSID, quiet=FALSE)
   mayo.score <- mayo.ad.ctl.separation(RSID)
   rosmap.score <- rosmap.ad.ctl.separation(RSID)

     # now look closely at the data used in the calculation
   targetGene <- "MS4A4A"
   RSID <- "rs1582763"  # very strong in mayo, almost nothing in rosmap
   etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD")
   x.may <- mayo.ad.ctl.separation(RSID)
   x.ros <- rosmap.ad.ctl.separation(RSID)

   tbl.pt <- etx$get.rosmap.patient.table(NA)
   dim(tbl.pt)  # 3583 18
   dir <- system.file(package="EndophenotypeExplorer", "extdata", "clinical")
   tbl.clinical.mayo <- read.table(file.path(dir, "MayoRNAseq_individual_metadata.csv"),
                                   sep=",", header=TRUE, as.is=TRUE)

   tbl.std.mayo.pt <- etx$standardizeMayoPatientTable(tbl.clinical.mayo)
     # 370 mayo samples
   dim(tbl.std.mayo.pt) # 370 10
   etx$get.rna.matrix.codes()
   mtx.mayo.tcx <- etx$get.rna.matrix("sage-eqtl-tcx")
   dim(mtx.mayo.tcx)  # 17009 257.  also rna-seq mayo matrices with 264, 263, 262 columns
   length(intersect(tbl.std.mayo.pt$patientID, colnames(mtx.mayo.tcx)))  # 257
   length(intersect(tbl.std.mayo.pt$patientID, colnames(mtx.geno)))  # 349

   mayo.pt <- tbl

   dir <- system.file(package="EndophenotypeExplorer", "extdata", "clinical")
   full.path <- file.path(dir, "ROSMAP_clinical.csv")
   tbl.clinical.rosmap <- read.table(full.path, sep=",", header=TRUE, as.is=TRUE)
   checkEquals(dim(tbl.clinical.rosmap), c(3583, 18))

   tbl.std.rosmap.pt <- etx$standardizeRosmapPatientTable(tbl.clinical.rosmap)
   dim(tbl.std.rosmap.pt)  # 3583 10
   length(intersect(tbl.std.rosmap.pt$patientID, colnames(mtx.geno)))  # 349

   tbl.map <- etx$getIdMap()
   table(tbl.map$study)        #   mayo rosmap  sinai
                               #    872   1783   1371

   length(intersect(colnames(mtx.geno), subset(tbl.map, study=="rosmap")$sample))  #  1151
   length(intersect(colnames(mtx.geno), subset(tbl.map, study=="mayo")$sample))    #   349
   length(intersect(colnames(mtx.geno), subset(tbl.map, study=="sinai")$sample))   #   345
   sum(1151, 349, 345) # 1845
   dim(mtx.geno)       # 1 1894    49 samples unaccounted for
   setdiff(colnames(mtx.geno), c(subset(tbl.map, study=="rosmap")$sample,
                                 subset(tbl.map, study=="mayo")$sample,
                                 subset(tbl.map, study=="sinai")$sample))

     #  [1] "71729"       "71823"       "76354"       "76655"       "MAP15387421" "MAP22868024"
     #  [7] "MAP26637867" "MAP29629849" "MAP33332646" "MAP34726040" "MAP46246604" "MAP46251007"
     # [13] "MAP50103679" "MAP50104134" "MAP50104846" "MAP50106442" "MAP50106992" "MAP50108462"
     # [19] "MAP50301099" "MAP50302680" "MAP50402431" "MAP50409406" "MAP61344957" "MAP85980779"
     # [25] "MAP87264456" "MAP89164957" "MAP93787649" "MAP95330358" "MAP95453354" "ROS10524640"
     # [31] "ROS11430815" "ROS11697592" "ROS15114174" "ROS15738428" "ROS20225051" "ROS20251553"
     # [37] "ROS20254452" "ROS20275399" "ROS20327084" "ROS20376029" "ROS20945666" "ROS20946257"
     # [43] "ROS20990085" "ROS20998065" "ROS21001807" "ROS21112011" "ROS21113864" "ROS21274866"
     # [49] "ROS79590778"

    # percentage genotypes at any variant, for each study, ought to be roughly the same

   pt.rosmap <- intersect(colnames(mtx.geno), subset(tbl.map, study=="rosmap")$sample)
   pt.mayo <- intersect(colnames(mtx.geno), subset(tbl.map, study=="mayo")$sample)
   pt.sinai <- intersect(colnames(mtx.geno), subset(tbl.map, study=="sinai")$sample)

   percent.dists <- function(table){
      tbl.freq <- as.data.frame(table)
      tbl.freq$percent <- round(100 *tbl.freq$Freq/sum(tbl.freq$Freq), digits=2)
      tbl.freq
      } # percent.dists

} # explore.L.shaped.distribution
#----------------------------------------------------------------------------------------------------
explore.braak.distribution <- function()
{
  tbl.r.pt <- etx$get.rosmap.patient.table(NA)
  dim(tbl.r.pt)
  tbl.m.pt <- etx$get.mayo.patient.table(NA)
  dim(tbl.m.pt)
  par(mfrow=c(1,2))
  hist(tbl.m.pt$Braak, breaks=0:6, ylim=c(0,600), main="Mayo Braak Scores n=370")
  hist(tbl.r.pt$braaksc, main="ROSMAP Braak Scores n=3583", breaks=0:6, ylim=c(0,600))


} # explore.braak.distribution
#----------------------------------------------------------------------------------------------------

# shared.sep.sig.list <- lapply(rsids.shared, rosmap.ad.ctl.separation)


