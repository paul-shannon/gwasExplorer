library(EndophenotypeExplorer)
targetGene <- "TFAM"
#----------------------------------------------------------------------------------------------------
ampad.eqtls <- function()
{
   etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD", initialize.snpLocs=FALSE)
   tbl.eQTL <- etx$get.ampad.EQTLsForGene()
   filename <- sprintf("tbl.eqtls.ampad.%s.%s.RData", targetGene,  gsub(" ", ".", Sys.time()))
   save(tbl.eQTL, file=filename)


                                        #         8740         4429           28
   viz <- FALSE
   if(viz){
      table(tbl.eQTL$study)  #   ampad-mayo ampad-rosmap         GTEx
      tbl.track <- subset(tbl.eQTL, study=="ampad-rosmap")[, c("chrom", "hg38", "hg38", "pvalue", "rsid")]
      colnames(tbl.track) <- c("chrom", "start", "end", "score", "rsid")
      tbl.track$start <- tbl.track$start - 1
      tbl.track$score <- -log10(tbl.track$score)
      track <- DataFrameQuantitativeTrack("rosmap dlpfc", tbl.track, color="random", autoscale=TRUE)
      displayTrack(igv, track)
      }

} # ampad.eqtls
#----------------------------------------------------------------------------------------------------
library(ebi.eqtls)
selected.studies <- c("GTEx_V8.Brain_Amygdala",
                      "GTEx_V8.Brain_Anterior_cingulate_cortex_BA24",
                      "GTEx_V8.Brain_Caudate_basal_ganglia",
                      "GTEx_V8.Brain_Cerebellar_Hemisphere",
                      "GTEx_V8.Brain_Cerebellum",
                      "GTEx_V8.Brain_Cortex",
                      "GTEx_V8.Brain_Frontal_Cortex_BA9",
                      "GTEx_V8.Brain_Hippocampus",
                      "GTEx_V8.Brain_Hypothalamus",
                      "GTEx_V8.Brain_Nucleus_accumbens_basal_ganglia",
                      "GTEx_V8.Brain_Putamen_basal_ganglia",
                      "GTEx_V8.Brain_Spinal_cord_cervical_c-1",
                      "GTEx_V8.Brain_Substantia_nigra")[c(4,5,6,7,8,9)]

chrom <- "chr10"
tss <- 58385407
#shoulder <- 100
shoulder <- 1000000
start <- tss - shoulder
end <- tss + shoulder
f <- function(study){
    ee <- ebi.eqtls$new()
    ee$fetch.eqtls.in.chunks(chrom=chrom,
                             start=start,
                             end=end,
                             study=study,
                             simplify=TRUE,
                             chunk.size=5000)
    }
x <- lapply(selected.studies, function(study) f(study))  # 7 jun 3:27pm
length(x)
names(x) <- selected.studies
tbl.all <- do.call(rbind, x)
rownames(tbl.all) <- NULL
dim(tbl.all)
filename <- sprintf("tbl.eqtls.%d.tissues-%s:%d-%d.%s",
                    length(x), chrom, start, end, gsub(" ", ".", Sys.time()))

save(tbl.all, file=filename)

further.work <- FALSE
if(further.work){
    tbl.all <- subset(tbl.all, gene==targetGene)
    print(dim(tbl.all))
    as.data.frame(sort(table(tbl.all$id)))


    if(!exists("igv"))
        igv <- start.igv(targetGene, "hg38")

    for(tissue in selected.studies){
        tbl.track <- subset(tbl.all, id==tissue)[, c("chrom", "hg38", "hg38", "pvalue", "beta", "rsid")]
        colnames(tbl.track)[1:3] <- c("chrom", "start", "end")
        tbl.track$start <- tbl.track$start - 1
        tbl.track$score <- -log10(tbl.track$pvalue)
        tissue.name.short <- sub("GTEx_V8.Brain_", "", tissue, fixed=TRUE)
        title <- sprintf("%s -log10(pval)", tissue.name.short)
        track <- DataFrameQuantitativeTrack(title,
                                            tbl.track[, c("chrom", "start", "end", "score", "rsid")],
                                            autoscale=TRUE, color="random")
        displayTrack(igv, track)
        tbl.track$score <- tbl.track$score * tbl.track$beta
        title <- sprintf("%s -log10(pval) * beta", tissue.name.short)
        track <- DataFrameQuantitativeTrack(title,
                                            tbl.track[, c("chrom", "start", "end", "score", "rsid")],
                                            autoscale=TRUE, color="random")
        displayTrack(igv, track)
        } # for tissue

    } # further.work


