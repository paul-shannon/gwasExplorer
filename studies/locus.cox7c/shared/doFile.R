#----------------------------------------------------------------------------------------------------
do.file <- function(file)
{
   print(file)
   print(load(file))
   short.tissue.name <- sub("GTEx_V8.Brain_", "", brain.tissue)
   message(sprintf("=== survey %s, %d tfs in trena model, %d motif breaks, %d tbl.scored rows",
                   short.tissue.name, nrow(tbl.trena), nrow(tbl.breaks), nrow(tbl.scored)))
   if(nrow(tbl.scored) == 0){
       return(data.frame())
       }

   tbls.tfs <- list()

   breaks.coi <- c("chrom", "start", "end", "SNP_id", "pctRef", "pctAlt", "pctDelta", "geneSymbol")
   nrow(tbl.breaks)
   tfs <- head(tbl.scored$gene, n=10)
   for(TF in tfs){
          # just the tfbs for current TF
      tbl.tms.TF <- subset(tbl.tms, tf==TF & ampad.eqtl)
      dups <- which(duplicated(tbl.tms.TF[, c(1:3)]))
      if(length(dups) > 0)
          tbl.tms.TF <- tbl.tms.TF[-dups,]
      if(nrow(tbl.breaks) == 0) next;
          # just the breaks for this TF and a meaninful change in binding
          # get the whole table width in case we want to look at the values
      tbl.breaks.TF <- subset(tbl.breaks, geneSymbol==TF & abs(pctDelta) > BREAK.ABS.PCT.DELTA.THRESHOLD) [, breaks.coi]
          # just the breaks which fall within previously selected tfbs
      if(nrow(tbl.breaks.TF) == 0)
          next;

      gr.tms <- GRanges(tbl.tms.TF)
      gr.breaks <- GRanges(tbl.breaks.TF)
      tbl.ov <- as.data.frame(findOverlaps(gr.breaks, gr.tms))
      if(nrow(tbl.ov) == 0)
          next;
      breaking.snps <- unique(tbl.breaks.TF$SNP_id[tbl.ov[, "queryHits"]])
          # how many TFBS for this TF?
      tfbs.count <- subset(tbl.trena, gene==TF)$tfbs
      if(nrow(tbl.tms.TF) == 0) next;
      if(nrow(tbl.breaks.TF) == 0) next;
      tbl.breaks.TF.eQTL <- unique(subset(tbl.breaks.TF, SNP_id %in% breaking.snps))
      dups <- which(duplicated(tbl.breaks.TF.eQTL[, 1:4]))
      if(length(dups) > 0)
        tbl.breaks.TF.eQTL <- tbl.breaks.TF.eQTL[-dups,]
      motifBreak.score <- with(tbl.breaks.TF.eQTL, sum(abs(pctDelta)) * 100)
      browser()
      tbl.eqtl.TF <- subset(eqtls$ampad, rsid %in% breaking.snps) # & study=="ampad-rosmap")
      eqtl.score <-  sum(-log10(tbl.eqtl.TF$pvalue))
      trena.score <- with(subset(tbl.trena, gene==TF), (abs(betaLasso) * 100) + (rfNorm * 10))
      tfbs.count <- subset(tbl.trena, gene==TF)$tfbs
      score.mult <- round(max(trena.score * eqtl.score * motifBreak.score), digits=0)
      score.sum <- max(trena.score + eqtl.score + motifBreak.score) # * tfbs.score)
      #printf("%10s: %6.1f %6.1f %6.1f %d  %d", TF,
      #       trena.score, eqtl.score, motifBreak.score,
      #       as.integer(score.sum), as.integer(score.mult))
      #browser()
      tbl <- data.frame(targetGene=targetGene,
                        tf=TF,
                        tissue=short.tissue.name,
                        trena=trena.score,
                        tfbs.count=tfbs.count,
                        eqtl=eqtl.score,
                        motifBreak=motifBreak.score,
                        total.score=as.integer(score.mult),
                        rsids=paste(breaking.snps, collapse=" "),
                        stringsAsFactors=FALSE)
       name <- sprintf("%s-%s", short.tissue.name, TF)
       tbls.tfs[[name]] <- unique(tbl)
      } # for TF

   if(length(tbls.tfs) == 0)
       return(data.frame())
   tbl <- do.call(rbind, tbls.tfs)
   rownames(tbl) <- NULL
   new.order <- order(tbl$total, decreasing=TRUE)
   tbl[new.order,]

} # do.file
#----------------------------------------------------------------------------------------------------
