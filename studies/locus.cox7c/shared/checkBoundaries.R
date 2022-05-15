# looking for 2M region, approximately chr5:85,893,798-87,960,958
tbl.snpLocs <- get(load("snpLocs.RData"))
dim(tbl.snpLocs)  # 647319      4
min <- min(tbl.snpLocs$hg38, na.rm=TRUE)  # 85927379
max <- max(tbl.snpLocs$hg38, na.rm=TRUE)  # 85927379
(max-min)/1000   # 1999.993

eqtl.list <- get(load("tbl.eqtls.gtex.13brain.tissues.chr5:85,927,378-87,927,378.RData"))
length(eqtl.list) # 13
names(eqtl.list)
    #  [1] "Amygdala"
    #  [2] "Anterior_cingulate_cortex_BA24"
    #  [3] "Caudate_basal_ganglia"
    #  [4] "Cerebellar_Hemisphere"
    #  [5] "Cerebellum"
    #  [6] "Cortex"
    #  [7] "Frontal_Cortex_BA9"
    #  [8] "Hippocampus"
    #  [9] "Hypothalamus"
    # [10] "Nucleus_accumbens_basal_ganglia"
    # [11] "Putamen_basal_ganglia"
    # [12] "Spinal_cord_cervical_c-1"
    # [13] "Substantia_nigra"
lapply(tbl.eqtls, nrow)
   # eqtl span is not the full 2M, since cis-eQTLs tend to cluster near the target gene
min(eqtl.list[[1]]$hg38, na.rm=TRUE)   # [1] 86427516
max(eqtl.list[[1]]$hg38, na.rm=TRUE)   # [1] 87425406

tbl.fimo  <- get(load("tbl.fimo.COX7C.RData"))
dim(tbl.fimo) # 3511559       9
min(tbl.fimo$start) # 85927348
max(tbl.fimo$start) # 87927400


