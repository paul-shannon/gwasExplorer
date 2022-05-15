f <- "gwex.CCNH.GTEx_V8.Brain_Hippocampus.Sun.May.15.2022-07:46:37.RData"
print(load(f))
viz <- function()
{
   if(!exists("igv")){
      igv <- start.igv(targetGene, "hg38")
      shoulder <- 5e5
      showGenomicRegion(igv, sprintf("%s:%d-%d", tag.snp.chrom, tag.snp.hg38-shoulder, tag.snp.hg38+shoulder))
      tbl.track <- data.frame(chrom=tag.snp.chrom, start=tag.snp.hg38-1, end=tag.snp.hg38, stringsAsFactors=FALSE)
      track <- DataFrameAnnotationTrack(tag.snp, tbl.track, trackHeight=25, color="darkred")
      displayTrack(igv, track)
      threshold <- 0.2
      tbl.track <- subset(tbl.linkage, r2 >= threshold)[, c("chrom", "hg38", "hg38", "r2", "rsid")]
      colnames(tbl.track) <- c("chrom", "start", "end", "score", "rsid")
      tbl.track$chrom <- tag.snp.chrom
      tbl.track$start <- tbl.track$start - 1
      track.title <- sprintf("haploreg eur, r^2 > %3.2f", threshold)
      track <- DataFrameQuantitativeTrack(track.title, tbl.track, autoscale=FALSE, min=0, max=1, color="black")
      displayTrack(igv, track)
      }

   etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD", initialize.snpLocs=TRUE)
   tbl.eqtls <- etx$get.ampad.EQTLsForGene()
   #tbl.eqtls <- subset(tbl.eqtls, study=="ampad-rosmap")
   dim(tbl.eqtls)
   head(tbl.eqtls, n=20)
   tail(tbl.eqtls, n=20)
   tbl.track <- subset(tbl.eqtls, pvalue < 1e-2 & study != "GTEx")[, c("chrom", "hg38", "hg38", "rsid", "pvalue")]
   dim(tbl.track)
   colnames(tbl.track)[2:3] <- c("start", "end")
   tbl.track$start <- tbl.track$start - 1
   colnames(tbl.track) <- NULL
   track <- GWASTrack("ampad", tbl.track, trackHeight=100) # chrom.col=0, pos.col=0, pval.col=0)
   displayTrack(igv, track)

       #--------------------------------------
       # eqtls from gtex, one track per gene
       #--------------------------------------

   tbl.eqtl.all <- do.call(rbind, eqtl.list)
   dim(tbl.eqtl.all)
   genes <- sort(unique(tbl.eqtl.all$gene))
   length(genes) # 9
   for(goi in genes) {
      tbl.track <- subset(tbl.eqtl.all, gene==goi)[, c("chrom", "hg38", "hg38", "rsid", "pvalue")]
      rownames(tbl.track) <- NULL
      colnames(tbl.track)[2:3] <- c("start", "end")
      tbl.track$start <- tbl.track$start - 1
      tbl.track$chrom <- paste0("chr", tbl.track$chrom)
      deleters <- which(is.na(tbl.track$start))
      if(length(deleters) > 0)
          tbl.track <- tbl.track[-deleters,]
      colnames(tbl.track) <- NULL
      track <- GWASTrack(goi, tbl.track, trackHeight=100) # chrom.col=0, pos.col=0, pval.col=0)
      displayTrack(igv, track)
      } # for goi

   removeTracksByName(igv,
                      getTrackNames(igv)[5:9]
                      )


      # now broken motifs of tfs in model
   tfs.oi <- tbl.trena$gene[1:20]
   for(TF in tfs.oi){
      tbl.breaks.sub <- subset(tbl.breaks.easy, geneSymbol==TF)[, c("chrom", "start", "end", "pctDelta")]
      tbl.breaks.sub$pctDelta <- -1.0 * tbl.breaks.sub$pctDelta
      if(nrow(tbl.breaks.sub) > 0){
         track <- DataFrameQuantitativeTrack(TF, tbl.breaks.sub, color="random", autoscale=FALSE, min=-0.2, max=0.2)
         displayTrack(igv, track)
         }
      }

   dim(tbl.tms.sub)
   tfs <- sort(unique(tbl.tms.sub$tf))
   length(tfs)


   for(TF in tfs){
      tbl.track <- subset(tbl.tms.sub, tf==TF)[, c("chrom", "start", "end", "cor.all")]
      track <- DataFrameQuantitativeTrack(TF, tbl.track, color="random", autoscale=FALSE, min=-1, max=1)
      displayTrack(igv, track)
      } # for TF

} # viz
