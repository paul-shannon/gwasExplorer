if(!exists("igv"))
   igv <- start.igv("all")

if(!exists("tbl.gwas3")){ # this 235 row, 215 rsid table, combines posthuma, schartzentruber, bellenquez
    data.dir <- "~/github/TrenaProjectAD/inst/extdata/gwasLoci"
    file.exists(data.dir)
    filename <- "tbl.gwas.3studies.235x5.RData"
    full.path <- file.path(data.dir, filename)
    tbl.gwas3 <- get(load(full.path))
    stopifnot(colnames(tbl.gwas3)[c(1,2,5)] == c("chrom", "hg38", "pvalue"))
    head(tbl.gwas3)
    }

if(!exists("tbl.gwascat.ad")){
   data.dir <- system.file(package="igvR", "extdata", "gwas")
   file.exists(data.dir)
   file <- "alzheimerSubsetOfGWASCatatalog-12jun2022.RData"
   full.path <- file.path(data.dir, file)
   file.exists(full.path)
   tbl.gwascat.ad <- get(load(full.path))
   stopifnot(colnames(tbl.gwascat.ad)[c(1,2,33)] == c("seqnames", "start", "P.VALUE"))
   }

if(!exists("tbl.rsids")){  # with mayo & rosmap scores for association of each variant to extreme braak score
   tbl.rsids <- get(load("tbl.rsids.summary.Tue.May.31.2022-13:26:55-hetAndHomSeparate.RData"))
   }

#----------------------------------------------------------------------------------------------------
showTracks <- function()
{
    showGenomicRegion(igv, "all")
    tbl.track <- GWASTrack("AD gwascat", tbl.gwascat.ad, chrom.col=1, pos.col=2, pval.col=33)
    displayTrack(igv, tbl.track)

    rosmap.outliers <- rownames(subset(tbl.rsids, rosmap.f > 2))
    tbl.track <- subset(tbl.gwas3, rsid %in% rosmap.outliers)
    track <- GWASTrack("rosmap outliers", tbl.track,  chrom.col=1, pos.col=2, pval.col=5)
    displayTrack(igv, track)

    mayo.outliers <- rownames(subset(tbl.rsids, mayo.f > 2))
    tbl.track <- subset(tbl.gwas3, rsid %in% mayo.outliers)
    track <- GWASTrack("mayo outliers", tbl.track,  chrom.col=1, pos.col=2, pval.col=5)
    displayTrack(igv, track)



} # showTracks
#----------------------------------------------------------------------------------------------------



