library(EndophenotypeExplorer)
library(ADvariantExplorer)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(RUnit)
#----------------------------------------------------------------------------------------------------
targetGene <- "ADAMTS4"
tag.snp <- "rs4575098"
tag.snp.chrom <- "chr1"
tag.snp.hg19 <- 161155392
tag.snp.hg38 <- 161185602
if(!exists("tbl.snpLocs")){
    if(!file.exists("snpLocs.RData")){
        msg <- "need to create tbl.snpLocs for region, save as .RData"
        stop(msg)
        }
    load("snpLocs.RData")
    }
#----------------------------------------------------------------------------------------------------
build.snpLocs.fromScratch <- function(shoulder=1e6)
{
   message("")
   message(sprintf("--- about to build tbl.snpLocs from scratch, width=%d", shoulder))
   message("")


           #--------------
           # hg38 first
           #--------------

      gr <- GRanges(seqnames=sub("chr", "", tag.snp.chrom),
                    IRanges(start=tag.snp.hg38-shoulder, end=tag.snp.hg38+shoulder))
      gr.snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, gr)
      tbl.snpLocs.hg38 <- as.data.frame(gr.snps)[, c(1,2,4)]
      colnames(tbl.snpLocs.hg38) <- c("chrom", "hg38", "rsid")
      tbl.snpLocs.hg38$chrom <- as.character(tbl.snpLocs.hg38$chrom)
      dim(tbl.snpLocs.hg38)

           #--------------
           # hg19
           #--------------

      gr <- GRanges(seqnames=sub("chr", "", tag.snp.chrom),
                    IRanges(start=tag.snp.hg19-shoulder, end=tag.snp.hg19+shoulder))
      gr.snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP144.GRCh37, gr)
      tbl.snpLocs.hg19 <- as.data.frame(gr.snps)[, c(1,2,4)]
      colnames(tbl.snpLocs.hg19) <- c("chrom", "hg19", "rsid")
      tbl.snpLocs.hg19$chrom <- as.character(tbl.snpLocs.hg19$chrom)
      dim(tbl.snpLocs.hg19)
      tbl.snpLocs <- merge(tbl.snpLocs.hg38, tbl.snpLocs.hg19, by=c("rsid", "chrom"), all=TRUE)
      rownames(tbl.snpLocs) <- tbl.snpLocs$rsid
      coi <- c("chrom", "hg19", "hg38", "rsid")
      tbl.snpLocs <- tbl.snpLocs[, coi]
      message(sprintf("saving tbl.snpLocs with %d rows", nrow(tbl.snpLocs)))
      save(tbl.snpLocs, file="snpLocs.RData")
      } # function

#----------------------------------------------------------------------------------------------------
if(!exists("etx")){
   etx <- EndophenotypeExplorer$new(targetGene, "hg19", vcf.project="AMPAD", initialize.snpLocs=TRUE)
   tbl.alfa <- etx$getAggregatedAlleleFrequencies(tag.snp)
     #          rsid ref       population    T     C total    T.freq    C.freq   min.freq
     # 1  rs76928645   C         European 1755 12531 14286 12.284754  87.71525  12.284754
     # 2  rs76928645   C   African Others    1   113   114  0.877193  99.12281   0.877193
     # 3  rs76928645   C       East Asian    0    86    86  0.000000 100.00000 100.000000
     # 4  rs76928645   C African American   52  2780  2832  1.836158  98.16384   1.836158
     # 5  rs76928645   C Latin American 1   11   135   146  7.534247  92.46575   7.534247
     # 6  rs76928645   C Latin American 2   45   565   610  7.377049  92.62295   7.377049
     # 7  rs76928645   C      Other Asian    0    26    26  0.000000 100.00000 100.000000
     # 8  rs76928645   C      South Asian    0    98    98  0.000000 100.00000 100.000000
     # 9  rs76928645   C          African   53  2893  2946  1.799050  98.20095   1.799050
     # 10 rs76928645   C            Asian    0   112   112  0.000000 100.00000 100.000000
     # 11 rs76928645   C            Total 1922 16968 18890 10.174696  89.82530  10.174696
     # 12 rs76928645   C            Other   58   634   692  8.381503  91.61850   8.381503
   }
if(!exists("avx")){
   avx <- ADvariantExplorer$new(targetGene, "chr4", 23564566, 24179593)
   }
#----------------------------------------------------------------------------------------------------
ampad.eqtls <- function()
{
    tbl.eqtls <- etx$get.ampad.EQTLsForGene()
    dim(tbl.eqtls)
    head(tbl.eqtls, n=20)
    tail(tbl.eqtls, n=20)
    tbl.track <- subset(tbl.eqtls, pvalue < 1e-2 & study != "GTEx")[, c("chrom", "hg19", "hg19", "rsid", "pvalue")]
    dim(tbl.track)
    colnames(tbl.track)[2:3] <- c("start", "end")
    tbl.track$start <- tbl.track$start - 1
    colnames(tbl.track) <- NULL
    track <- GWASTrack("ampad", tbl.track, trackHeight=100) # chrom.col=0, pos.col=0, pval.col=0)
    displayTrack(igv, track)

} # ampad.eqtls
#----------------------------------------------------------------------------------------------------
fetch.eqtls <- function(chrom, start, end, study, simplify, chunk.size)
{
  roi.width <- 1 + end - start

  if(roi.width <= chunk.size){
     message(sprintf("--- fetch.eqtls, just one chunk"))
     tbl <- avx$geteQTLsByLocationAndStudyID(chrom, start, end, study, simplify=simplify)

  } else {
     boundaries.needed <- 1 + (roi.width %/% chunk.size)
     starts <- as.integer(seq(from=start, to=end, length.out=boundaries.needed))
     ends <- starts[-1]
     starts <- starts[-(length(starts))]
     tbls <- list()
     intervals <- length(starts)
     message(sprintf("==== fetch.eqtls, %d chunks", intervals))
     for(i in seq_len(intervals)){
        message(sprintf("--- fetching chunk %2d/%d for %s", i, intervals, study))
        tbl.chunk <- avx$geteQTLsByLocationAndStudyID(chrom,
                                                      as.integer(starts[i]),
                                                      as.integer(ends[i]),
                                                      study,
                                                      targetGene.only=FALSE,
                                                      simplify=simplify)
        tbls[[i]] <- tbl.chunk
        } # for i
     tbl <- do.call(rbind, tbls)
     } # else

  invisible(tbl)

} # fetch.eqtls
#----------------------------------------------------------------------------------------------------
test_fetch.eqtls <- function()
{
   tbl.small.no.chunks <- fetch.eqtls("chr4", start=23738124, end=23786367,
                                       study="GTEx_V8.Brain_Cerebellar_Hemisphere",
                                       simplify=TRUE, chunk.size=100000)
   checkTrue(nrow(tbl.small.no.chunks) > 200)

   tbl.small.20k.chunks <- fetch.eqtls("chr4", start=23738124, end=23786367,
                                       study="GTEx_V8.Brain_Cerebellar_Hemisphere",
                                       simplify=TRUE, chunk.size=20000)
   checkEquals(nrow(tbl.small.no.chunks), nrow(tbl.small.20k.chunks))

   tbl.small.10k.chunks <- fetch.eqtls("chr4", start=23738124, end=23786367,
                                       study="GTEx_V8.Brain_Cerebellar_Hemisphere",
                                       simplify=TRUE, chunk.size=10000)
   checkEquals(nrow(tbl.small.no.chunks), nrow(tbl.small.10k.chunks))

   tbl.big <- fetch.eqtls("chr4", start=22412644, end=25189900,
                          study="GTEx_V8.Brain_Cerebellar_Hemisphere",
                          chunk.size=100000,
                          simplify=TRUE)
   viz.big <- function(){
      tbl.sub <- subset(tbl.big, pvalue < 1.5)
      dim(tbl.sub)
      deleters <- which(is.na(tbl.sub$rsid))
      length(deleters)
      if(length(deleters) > 0)
          tbl.sub <- tbl.sub[-deleters,]
      tbl.snpLocs <- as.data.frame(snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, tbl.sub$rsid, ifnotfound="drop"))
      tbl.merged <- merge(tbl.sub, tbl.snpLocs, by.x="rsid", by.y="RefSNP_id")
      tbl.track <- tbl.merged[, c("seqnames", "pos", "pos", "rsid", "pvalue")]
      colnames(tbl.track) <- c("chrom", "start", "end", "rsid", "pvalue")
      fivenum(tbl.track$pvalue)
      tbl.track$chrom <- sprintf("chr%s", tbl.track$chrom)
      tbl.track$start <- tbl.track$start - 1
      dim(tbl.track)

      colnames(tbl.track) <- NULL
      track <- GWASTrack("big test 1.0", tbl.track, trackHeight=100)
      #track <- DataFrameQuantitative("tbl.
      displayTrack(igv, track)
      }

} # test_fetch.eqtls
#----------------------------------------------------------------------------------------------------
fetch.gtex.eqtls <- function(shoulder=1e6)
{
   loc <- list(chrom=tag.snp.chrom,
               start=tag.snp.hg38 - shoulder,
               end=tag.snp.hg38 + shoulder,
               string = sprintf("%s:%d-%d", tag.snp.chrom, tag.snp.hg38 - shoulder, tag.snp.hg38 + shoulder))

   message(sprintf("--- fetching gtex brain eqtls, region size: %dk", as.integer(loc$width/1000)))

   avx <- with(loc, ADvariantExplorer$new(targetGene, chrom, start, end))
   tbl.eCat <- avx$geteQTLSummary()
   dim(tbl.eCat)   # 498 12

   brain.geneExpression.studies <- subset(tbl.eCat, grepl("brain", tissue_label, ignore.case=TRUE) &
                                                    quant_method=="ge")$unique_id
   gtex.v8.brain.studies <- grep("GTEx_V8", brain.geneExpression.studies, value=TRUE)
   length(gtex.v8.brain.studies) # 13
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
                         "GTEx_V8.Brain_Substantia_nigra")[c(3,4,5,6,8,9)]

   chunk.size <- 10000

   for(study in selected.studies){
      tbl.eqtl <- with(loc, fetch.eqtls(chrom, start, end,
                                        study=study,
                                        chunk.size=chunk.size,
                                        simplify=TRUE))
      study.title <- sub("GTEx_V8.Brain_", "", study)
      study.filename <- sprintf("tbl.eqtls.gtex.%s.%s.RData", study.title, loc$string)
      dim(tbl.eqtl)
      deleters <- which(is.na(tbl.eqtl$rsid))
      length(deleters)
      if(length(deleters) > 0)
          tbl.eqtl <- tbl.eqtl[-deleters,]
      tbl.eqtl <- subset(tbl.eqtl, pvalue < 0.05)
      tbl.eqtl <- merge(tbl.eqtl, tbl.snpLocs, by="rsid", all.x=TRUE)
      tbl.eqtl$score <- -log10(tbl.eqtl$pvalue) * abs(tbl.eqtl$beta)
      dim(tbl.eqtl)
      new.order <- order(tbl.eqtl$score, decreasing=TRUE)
      tbl.eqtl <- tbl.eqtl[new.order,]
      tbl.eqtl$name <- sprintf("%s-%s", tbl.eqtl$rsid, tbl.eqtl$gene)
      rownames(tbl.eqtl) <- NULL
      deleters <- grep("RP11", tbl.eqtl$gene)
      if(length(deleters) > 0)
          tbl.eqtl <- tbl.eqtl[-deleters,]
      tbls.eqtl[[study.title]] <- tbl.eqtl
      message("")
      message(sprintf("---- saving %d %s eqtls to %s", nrow(tbl.eqtl), study.title, filename))
      message("")
      save(tbl.eqtl, file=study.filename)
      examine <- FALSE
      if(examine){
         tbl.eqtl.strong <- subset(tbl.eqtl, score > 0.2)
         gois <- sort(unique(tbl.eqtl.strong$gene))
         for(goi in gois){
             tbl.track <- subset(tbl.eqtl.strong, gene==goi)[, c("chrom", "hg38", "hg38", "score", "rsid")]
             deleters <- which(is.na(tbl.track$hg38))
             if(length(deleters) > 0)
             tbl.track <- tbl.track[-deleters,]
             colnames(tbl.track) <- c("chrom", "start", "end", "score", "name")
             tbl.track$chrom <- sprintf("chr%s", tbl.track$chrom)
             tbl.track$start <- tbl.track$start - 1
             dim(tbl.track)
             track.title <- sprintf("%s %s", study.title, goi)
             track <- DataFrameQuantitativeTrack(track.title, tbl.track, autoscale=FALSE,
                                                 min=0, max=5, color="random",
                                                 trackHeight=25)
             displayTrack(igv, track)
             } # for goi
          } # if examine
      } # for each brain study

   fs <- list.files(".", "tbl.*RData")
   tbls.eqtl <- list()
   filename <- sprintf("tbl.eqtls.gtex.%d.brain.tissues.%s.RData", length(selected.studies), loc$string)

   for(f in fs){
      tissue <- sub("tbl.eqtls.gtex\\.", "", f)
      tissue <- sub("\\.chr.*RData", "", tissue)
      tbl.eqtl <- get(load(f))
      tbl.eqtl$tissue <- tissue
      tbls.eqtl[[tissue]] <- tbl.eqtl
      } # for f


   tbl.all <- do.call(rbind, tbls.eqtl)
   save(tbl.all, file=filename)

} # gtex.eqtls
#----------------------------------------------------------------------------------------------------
