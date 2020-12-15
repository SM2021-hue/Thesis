rm(list=ls())
setwd("~/Documents/Thesis_Suneeta_Modekurty/ad_psp_try")
## -----------------------------------------------------------------------------
library(knitr) 


## -----------------------------------------------------------------------------
library(pacman)

p_load(plyr,
       tidyverse,
       VennDiagram,
       pheatmap,
       gplots,
       gridExtra,
       EnhancedVolcano,
       venn,
       openxlsx,
       gprofiler2,
       corrplot,
       PerformanceAnalytics,
       Hmisc,
       edgeR,
       limma,
       Glimma,
       gplots,
       org.Mm.eg.db,
       RColorBrewer)

p_load_gh("kassambara/ggpubr")
source("./scripts/TAMPOR.R") 


## -----------------------------------------------------------------------------
Mayo_Proteomics_TC_traits <- read.csv("./data_prot/Mayo_Proteomics_TC_traits.csv",
                                      sep = ",",
                                      header = TRUE,
                                      stringsAsFactors = FALSE) 


## -----------------------------------------------------------------------------
Mayo_Proteomics_TC_proteinoutput <- read.table("./data_prot/Mayo_Proteomics_TC_proteinoutput.txt",
                                               sep = "\t",
                                               header = TRUE,
                                               stringsAsFactors = FALSE) 


## -----------------------------------------------------------------------------
## Fitler out Potential.contaminant and Reverse "+"
dat <- Mayo_Proteomics_TC_proteinoutput %>%
  filter(Potential.contaminant != "+") %>%
  filter(Reverse != "+") 


## -----------------------------------------------------------------------------
rownames(dat) <- dat$Protein.IDs 


## -----------------------------------------------------------------------------
## Select LFQ + egis columns only
dat.m <- dat %>%
  dplyr::select(starts_with("LFQ")) %>%
  dplyr::select(!contains("egis"))

dim(dat.m) 


## -----------------------------------------------------------------------------
traits <- Mayo_Proteomics_TC_traits[1:4] 


## -----------------------------------------------------------------------------
rownames(traits) <- traits$Samples_Simple %>%
  gsub("_", "\\.", .)

traits.e <- traits[1:215, ] # egis
traits.m <- traits[c(1:200, 216:230), ] # mgis 


## -----------------------------------------------------------------------------
dat.mt <- dat.m

## Renaming scheme:
## egis samples: "LFQ.intensity.mayo_b5_egis_24" -> "b5.24"
## non-egis:     "LFQ.intensity.mayo_b5_197_41" -> "b5.197"
colnames(dat.mt) <- dat.m %>%
  names %>%
  gsub("LFQ.intensity.mayo_", "", .) %>%
  gsub("mgis_", "", .) %>%
  gsub("([0-9]{3})_[0-9]{2}$", "\\1", .) %>% # to strip 2 last digits from non egis samples only
  gsub("_", "\\.", .) 


## -----------------------------------------------------------------------------
dat.mt <- dat.mt[, rownames(traits.m)] 


## -----------------------------------------------------------------------------
## see: https://github.com/edammer/TAMPOR/issues/2
dat.mt[dat.mt<=0]<- NA 


## -----------------------------------------------------------------------------
## ## ## Samples designated as GIS were removed before visualization of variance and MDS
## TAMPORlist.MGIS.all.true <- TAMPOR(dat.mt, traits.m,
##                                  noGIS=FALSE,
##                                  useAllNonGIS=TRUE, # This assumes
##                                                     # batches are
##                                                     # generally
##                                                     # trait-balanced
##                                                     # and randomized,
##                                                     # but GIS and
##                                                     # non-GIS samples
##                                                     # can be grossly
##                                                     # different.
##                                  batchPrefixInSampleNames=TRUE,
##                                  GISchannels=c("01","23", "44", "45"),
##                                  samplesToIgnore = NULL,
##                                  removeGISafter = TRUE,
##                                  parallelThreads=2,
##                                  outputSuffix="MGIS")

## saveRDS(TAMPORlist.MGIS.all.true, "./results_prot/TAMPORlist.MGIS.all.true.rds")


TAMPORlist.MGIS.all.true <- readRDS("./results_prot/TAMPORlist.MGIS.all.true.rds")


## plotMDS(TAMPORlist.GIS$cleanRelAbun %>% log2, col = labels2colors(colnames(TAMPORlist.GIS$cleanRelAbun)))

cra <- TAMPORlist.MGIS.all.true$cleanRelAbun

cra.log <- log2(cra)
 


## -----------------------------------------------------------------------------
conditions.col <- traits$Diagnosis[1:200] %>%
  sub("Control", "green", .) %>%
  sub("AD", "red", .) %>%
  sub("PSP", "blue", .)

traw.log <- log2(dat.mt[1:200])
colnames(traw.log) <- traits[1:200,4] %>% make.unique

tcra.log <- cra.log
colnames(tcra.log) <- traits[1:200,4] %>% make.unique

require(limma)

## opar <- par(no.readonly=TRUE)
png("./figs_prot/mds_tampor_conditions2.png", height = 2000, width = 6000, units = "px")
## tiff("./figs_prot/mds_tampor_conditions2.tiff", height = 2000, width = 6000, units = "px")
par(mfrow = c(1, 2), cex = 5)

plotMDS(traw.log, col = conditions.col)
title("log2 (raw)")

plotMDS(tcra.log, col = conditions.col)
title("After TAMPOR (mgis)")

dev.off()
## par(opar)
## dev.off()

rm(traw.log, tcra.log)
 


## ----echo=FALSE, label="", fig.show='asis', fig.cap="Number of samples: AD = 84, PSP = 85, Control = 31"----
include_graphics("./figs_prot/mds_tampor_conditions2.png") 


## -----------------------------------------------------------------------------
pv_cut <- 5e-2
fc_cut <- 0.3 


## -----------------------------------------------------------------------------
traits <- Mayo_Proteomics_TC_traits %>%
  filter(Diagnosis != "EmoryGlobalInternalStandard" &
           Diagnosis != "MayoGlobalInternalStandard") %>%
  dplyr::select(1:4)


control_cols <- traits %>%
  mutate(n = rownames(.) %>% as.numeric) %>%
  filter(Diagnosis == "Control") %>%
  .$n

ad_cols <- traits %>%
  mutate(n = rownames(.) %>% as.numeric) %>%
  filter(Diagnosis == "AD") %>%
  .$n

psp_cols <- traits %>%
  mutate(n = rownames(.) %>% as.numeric) %>%
  filter(Diagnosis == "PSP") %>%
  .$n

## mean normalisation
## replace_with_mean <- function(x) replace(x, x==0, mean(x))
## median normalisation
replace_with_median <- function(x) replace(x, is.na(x), median(x, na.rm = TRUE))

##
crad <- as.data.frame(cra)
## crad <- as.data.frame(cra[-c(23, 53, 84)])

control <- crad %>%
  dplyr::select(all_of(control_cols)) %>%
  apply(., 1, replace_with_median) %>%
  t %>%
  data.frame

ad <- crad %>%
  dplyr::select(all_of(ad_cols)) %>%
  apply(., 1, replace_with_median) %>%
  t %>%
  data.frame

psp <- crad %>%
  dplyr::select(all_of(psp_cols)) %>%
  apply(., 1, replace_with_median) %>%
  t %>%
  data.frame

data_imp <- cbind.data.frame(control, ad, psp) 


## -----------------------------------------------------------------------------
# Quantile normalisation : the aim is to give different distributions the
# same statistical properties
quantile_normalisation <- function(df){

  # Find rank of values in each column
  df_rank <- map_df(df,rank,ties.method="average")
  # Sort observations in each column from lowest to highest
  df_sorted <- map_df(df,sort)
  # Find row mean on sorted columns
  df_mean <- rowMeans(df_sorted)

  # Function for substiting mean values according to rank
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }

  # Replace value in each column with mean according to rank
  df_final <- map_df(df_rank,index_to_mean, my_mean=df_mean)

  return(df_final)
} 


## -----------------------------------------------------------------------------
data_norm <- data_imp %>%
  quantile_normalisation() 


## -----------------------------------------------------------------------------
data_imp %>% ggplot() +
  geom_density(aes(b1.001 %>% log2), col = "red") +
  geom_density(aes(b2.002 %>% log2), col = "green") +
  geom_density(aes(b3.011 %>% log2), col = "blue") +
  geom_density(aes(b4.005 %>% log2), col = "magenta") +
  geom_density(aes(b5.175 %>% log2)) +
  ggtitle("Dstribution before normalization (5 samples from different batches)") +
  xlab("log2 Intensity") +
  ylab("Density")

data_norm %>% ggplot() +
  geom_density(aes(b1.001 %>% log2), col = "red", lwd=2.5) +
  geom_density(aes(b2.002 %>% log2), col = "green", lwd=2) +
  geom_density(aes(b3.011 %>% log2), col = "blue", lwd = 1.5) +
  geom_density(aes(b4.005 %>% log2), col = "magenta", lwd = 1) +
  geom_density(aes(b5.175 %>% log2), lwd=0.5) +
  ggtitle("Dstribution after normalization (The same 5 samples)") +
  xlab("log2 Intensity") +
  ylab("Density") 


## -----------------------------------------------------------------------------
traits <- Mayo_Proteomics_TC_traits %>%
  filter(Diagnosis != "EmoryGlobalInternalStandard" &
           Diagnosis != "MayoGlobalInternalStandard") %>%
  dplyr::select(1:4)


control_cols <- traits %>%
  mutate(n = rownames(.) %>% as.numeric) %>%
  filter(Diagnosis == "Control") %>%
  .$n

ad_cols <- traits %>%
  mutate(n = rownames(.) %>% as.numeric) %>%
  filter(Diagnosis == "AD") %>%
  .$n

psp_cols <- traits %>%
  mutate(n = rownames(.) %>% as.numeric) %>%
  filter(Diagnosis == "PSP") %>%
  .$n

## rearrange columns in cra (optional)
cra <- cra %>%
  data.frame %>%
  dplyr::select(all_of(c(control_cols, ad_cols, psp_cols))) 


## -----------------------------------------------------------------------------
diagnosis <- c(rep("Control", length(control_cols)), rep("AD", length(ad_cols)), rep("PSP", length(psp_cols)))
design <- model.matrix( ~ diagnosis + 0) 


## -----------------------------------------------------------------------------
data_norm <- data.frame(data_norm)
rownames(data_norm) <- rownames(cra) %>% sub("^.+\\|.+\\|(.+)_HUMAN.?", "\\1", .) %>% make.unique

require(limma)
fit <- lmFit(data_norm %>% log2, design) 


## -----------------------------------------------------------------------------
AvsC <- makeContrasts(diagnosisAD - diagnosisControl, levels=design) 


## -----------------------------------------------------------------------------
fit_AvsC <-contrasts.fit(fit, AvsC)

bf_AvsC <- eBayes(fit_AvsC)

bf_AvsC %>% decideTests %>% summary 


## -----------------------------------------------------------------------------
AvsC_top <- topTable(bf_AvsC, adjust="BH", sort.by = "logFC", number = 1000, p.value = pv_cut, lfc = fc_cut)

write.csv(AvsC_top, "./results_prot/prot_ad_vs_control.csv")

knitr::kable(AvsC_top,
             caption = "Top DE proteins: AD vs Control")

 


## -----------------------------------------------------------------------------
PvsC <- makeContrasts(diagnosisPSP - diagnosisControl, levels=design) 


## -----------------------------------------------------------------------------
fit_PvsC <-contrasts.fit(fit, PvsC)

bf_PvsC <- eBayes(fit_PvsC)

bf_PvsC %>% decideTests %>% summary 


## -----------------------------------------------------------------------------
PvsC_top <- topTable(bf_PvsC, adjust="BH", sort.by = "logFC", number = 1000, p.value = pv_cut, lfc = fc_cut)

write.csv(PvsC_top, "./results_prot/prot_psp_vs_control.csv")

knitr::kable(PvsC_top,
             caption = "Top DE proteins: PSP vs Control") 


## -----------------------------------------------------------------------------
AvsP <- makeContrasts(diagnosisAD - diagnosisPSP, levels=design) 


## -----------------------------------------------------------------------------
fit_AvsP <-contrasts.fit(fit, AvsP)

bf_AvsP <- eBayes(fit_AvsP)

bf_AvsP %>% decideTests %>% summary 


## -----------------------------------------------------------------------------
AvsP_top <- topTable(bf_AvsP, adjust="BH", sort.by = "logFC", number = 1000, p.value = pv_cut, lfc = fc_cut)

write.csv(AvsP_top, "./results_prot/prot_ad_vs_psp.csv")

knitr::kable(AvsP_top,
             caption = "Top DE proteins: AD vs PSP") 


## -----------------------------------------------------------------------------
dat_fc_pv <- bind_cols(ProteinIDs = rownames(data_imp),
                       fc_control_ad = topTable(bf_AvsC, sort = "none", n = Inf)$logFC,
                       fc_control_psp = topTable(bf_PvsC, sort = "none", n = Inf)$logFC,
                       fc_ad_psp = topTable(bf_AvsP, sort = "none", n = Inf)$logFC,
                       ad_control_pv = topTable(bf_AvsC, sort = "none", n = Inf)$adj.P.Val,
                       psp_control_pv = topTable(bf_PvsC, sort = "none", n = Inf)$adj.P.Val,
                       ad_psp_pv = topTable(bf_AvsP, sort = "none", n = Inf)$adj.P.Val,
                       PIDshort = rownames(data_norm)) 


## ----fig.asp = 1--------------------------------------------------------------

AD_C_UP <- dat_fc_pv %>%
  filter(ad_control_pv <= pv_cut) %>%
  filter(fc_control_ad >= fc_cut) %>%
  arrange(-abs(fc_control_ad)) %>%
  dplyr::select(PIDshort, ad_control_pv, fc_control_ad)

AD_C_DOWN <- dat_fc_pv %>%
  filter(ad_control_pv <= pv_cut) %>%
  filter(fc_control_ad <= -fc_cut) %>%
  arrange(-abs(fc_control_ad)) %>%
  dplyr::select(PIDshort, ad_control_pv, fc_control_ad)

PSP_C_UP <- dat_fc_pv %>%
  filter(psp_control_pv <= pv_cut) %>%
  filter(fc_control_psp >= fc_cut) %>%
  arrange(-abs(fc_control_psp)) %>%
  dplyr::select(PIDshort, psp_control_pv, fc_control_psp)

PSP_C_DOWN <- dat_fc_pv %>%
  filter(psp_control_pv <= pv_cut) %>%
  filter(fc_control_psp <= -fc_cut) %>%
  arrange(-abs(fc_control_psp)) %>%
  dplyr::select(PIDshort, psp_control_pv, fc_control_psp)

AD_PSP_UP <- dat_fc_pv %>%
  filter(ad_psp_pv <= pv_cut) %>%
  filter(fc_ad_psp >= fc_cut) %>%
  arrange(-abs(fc_ad_psp)) %>%
  dplyr::select(PIDshort, ad_psp_pv, fc_ad_psp)

AD_PSP_DOWN <- dat_fc_pv %>%
  filter(ad_psp_pv <= pv_cut) %>%
  filter(fc_ad_psp <= -fc_cut) %>%
  arrange(-abs(fc_ad_psp)) %>%
  dplyr::select(PIDshort, ad_psp_pv, fc_ad_psp)

require(openxlsx)

write.xlsx(AD_C_UP,
           file = paste0("./results_prot/DE_genes_p-value_", pv_cut, "_log2FC_", fc_cut, ".xlsx"),
           append = TRUE,
           sheetName = "AD vs Control UP")

write.xlsx(AD_C_DOWN,
           file = paste0("./results_prot/DE_genes_p-value_", pv_cut, "_log2FC_", fc_cut, ".xlsx"),
           append = TRUE,
           sheetName = "AD vs Control DOWN")

write.xlsx(PSP_C_UP,
           file = paste0("./results_prot/DE_genes_p-value_", pv_cut, "_log2FC_", fc_cut, ".xlsx"),
           append = TRUE,
           sheetName = "PSP vs Control UP")

write.xlsx(PSP_C_DOWN,
           file = paste0("./results_prot/DE_genes_p-value_", pv_cut, "_log2FC_", fc_cut, ".xlsx"),
           append = TRUE,
           sheetName = "PSP vs Control DOWN")

write.xlsx(AD_PSP_UP,
           file = paste0("./results_prot/DE_genes_p-value_", pv_cut, "_log2FC_", fc_cut, ".xlsx"),
           append = TRUE,
           sheetName = "AD vs PSP UP")

write.xlsx(AD_PSP_DOWN,
           file = paste0("./results_prot/DE_genes_p-value_", pv_cut, "_log2FC_", fc_cut, ".xlsx"),
           append = TRUE,
           sheetName = "AD vs PSP DOWN") 


## -----------------------------------------------------------------------------
data_log2 <- data_norm %>% log2

d <- data_log2

## Set clear column names
colnames(d) <- c(rep("Control", 31), rep("AD", 84), rep("PSP", 85)) %>% make.unique

## add columns with short protein IDs, log2 folds,  adjusted p-values
data_all <- data.frame(d, dat_fc_pv)

rm(d) 


## -----------------------------------------------------------------------------
dat_filt_ad <- data_all %>% filter(abs(fc_control_ad) >= fc_cut &
                                     ad_control_pv <= pv_cut) %>%
  arrange(-abs(fc_control_ad)) %>%
  ## slice_head(n=25) %>%
  dplyr::select(PIDshort, (starts_with("AD.") | starts_with("Control")) ) 


## -----------------------------------------------------------------------------
## data.matrix() leaves numeric data as is.
ad_matrix <- data.matrix(dat_filt_ad)

rownames(ad_matrix) <- dat_filt_ad$PIDshort

ad_matrix <- ad_matrix[,-1] 


## -----------------------------------------------------------------------------
ad_scaled <- ad_matrix %>% t %>% scale %>% t 


## -----------------------------------------------------------------------------
ad_dist <- ad_scaled %>% t %>%
  dist(., method = "euclidean", diag = FALSE, upper = FALSE) 


## -----------------------------------------------------------------------------
ad_dist_prot <- ad_scaled %>% 
  dist(., method = "euclidean", diag = FALSE, upper = FALSE) 


## ----fig.asp = 1.5------------------------------------------------------------
ad_hclust <- hclust(ad_dist, method = "ward.D2", members = NULL)
ad_hclust_prot <- hclust(ad_dist_prot, method = "ward.D2", members = NULL)


png("./figs_prot/hclust_ad.png", width = 2500, height = 3000, res = 300)
par(mfrow = c(2,1))
ad_hclust %>% plot(cex = 0.4, main = "Conditions (Number of samples: AD = 84, Control = 31)")
ad_hclust_prot %>% plot(cex = 0.6, main = "Proteins (Number of samples: AD = 84, Control = 31)")
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_prot/hclust_ad.png") 


## ----fig.asp = 1--------------------------------------------------------------
# Set colours for heatmap, 25 increments
my_palette <- colorRampPalette(c("blue","white","red"))(n = 25)

png("./figs_prot/heatmap_ad.png", width = 2000, height = 2000, res = 300)
# Plot heatmap with heatmap.2
par(cex.main=0.75) # Shrink title fonts on plot
ad_scaled %>%
  # Plot heatmap
  gplots::heatmap.2(.,                     # Tidy, normalised data
                    Colv=as.dendrogram(ad_hclust),     # Experiments clusters in cols
                    Rowv=as.dendrogram(ad_hclust_prot),     # Protein clusters in rows
                    revC=TRUE,                  # Flip plot to match pheatmap
                    density.info="histogram",   # Plot histogram of data and colour key
                    trace="none",               # Turn of trace lines from heat map
                    col = my_palette,           # Use my colour scheme
                    cexRow=0.3,cexCol=0.2,      # Amend row and column label fonts
                    xlab = paste0("Number of samples: AD = ", length(ad_cols), ", Control = ", length(control_cols))
                    )
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_prot/heatmap_ad.png") 


## -----------------------------------------------------------------------------
dat_filt_psp <- data_all %>% filter(abs(fc_control_psp) >= fc_cut &
                                     psp_control_pv <= pv_cut) %>%
  arrange(-abs(fc_control_psp)) %>%
  ## slice_head(n=25) %>%
  dplyr::select(PIDshort, (starts_with("PSP.") | starts_with("Control")) ) 


## -----------------------------------------------------------------------------
## data.matrix() leaves numeric data as is.
psp_matrix <- data.matrix(dat_filt_psp)

rownames(psp_matrix) <- dat_filt_psp$PIDshort

psp_matrix <- psp_matrix[,-1] 


## -----------------------------------------------------------------------------
psp_scaled <- psp_matrix %>% t %>% scale %>% t 


## -----------------------------------------------------------------------------
psp_dist <- psp_scaled %>% t %>%
  dist(., method = "euclidean", diag = FALSE, upper = FALSE) 


## -----------------------------------------------------------------------------
psp_dist_prot <- psp_scaled %>%
  dist(., method = "euclidean", diag = FALSE, upper = FALSE) 


## ----fig.asp = 1.5------------------------------------------------------------
psp_hclust <- hclust(psp_dist, method = "ward.D2", members = NULL)
psp_hclust_prot <- hclust(psp_dist_prot, method = "ward.D2", members = NULL)


png("./figs_prot/hclust_psp.png", width = 2500, height = 3000, res = 300)
par(mfrow = c(2,1))
psp_hclust %>% plot(cex = 0.4, main = "Conditions (Number of samples: PSP = 85, Control = 31)")
psp_hclust_prot %>% plot(cex = 0.6, main = "Proteins (Number of samples: PSP = 85, Control = 31)")
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_prot/hclust_psp.png") 


## ----fig.asp = 1--------------------------------------------------------------
# Set colours for heatmap, 25 increments
my_palette <- colorRampPalette(c("blue","white","red"))(n = 25)

# Plot heatmap with heatmap.2
png("./figs_prot/heatmap_psp.png", width = 2000, height = 2000, res = 300)
par(cex.main=0.75) # Shrink title fonts on plot
psp_scaled %>%
  # Plot heatmap
  gplots::heatmap.2(.,                     # Tidy, normalised data
          Colv=as.dendrogram(psp_hclust),     # Experiments clusters in cols
          Rowv=as.dendrogram(psp_hclust_prot),     # Protein clusters in rows
          revC=TRUE,                  # Flip plot to match pheatmap
          density.info="histogram",   # Plot histogram of data and colour key
          trace="none",               # Turn of trace lines from heat map
          col = my_palette,           # Use my colour scheme
          cexRow=0.3, cexCol=0.3,     # Amend row and column label fonts
          xlab = paste0("Number of samples: PSP = ", length(psp_cols), ", Control = ", length(control_cols))
          )
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_prot/heatmap_psp.png") 


## -----------------------------------------------------------------------------
dat_filt_adpsp <- data_all %>% filter(abs(fc_ad_psp) >= fc_cut &
                                     ad_psp_pv <= pv_cut) %>%
  arrange(-abs(fc_ad_psp)) %>%
  ## slice_head(n=25) %>%
  dplyr::select(PIDshort, (starts_with("AD.") | starts_with("PSP.")) ) 


## -----------------------------------------------------------------------------
## data.matrix() leaves numeric data as is.
adpsp_matrix <- data.matrix(dat_filt_adpsp)

rownames(adpsp_matrix) <- dat_filt_adpsp$PIDshort

adpsp_matrix <- adpsp_matrix[,-1] 


## -----------------------------------------------------------------------------
adpsp_scaled <- adpsp_matrix %>% t %>% scale %>% t 


## -----------------------------------------------------------------------------
adpsp_dist <- adpsp_scaled %>% t %>%
  dist(., method = "euclidean", diag = FALSE, upper = FALSE) 


## -----------------------------------------------------------------------------
adpsp_dist_prot <- adpsp_scaled %>%
  dist(., method = "euclidean", diag = FALSE, upper = FALSE) 


## ----fig.asp = 1.5------------------------------------------------------------
adpsp_hclust <- hclust(adpsp_dist, method = "ward.D2", members = NULL)
adpsp_hclust_prot <- hclust(adpsp_dist_prot, method = "ward.D2", members = NULL)


png("./figs_prot/hclust_ad_psp.png", width = 2500, height = 3000, res = 300)
par(mfrow = c(2,1))
adpsp_hclust %>% plot(cex = 0.2, main = "Conditions (Number of samples: AD = 84, PSP = 85)")
adpsp_hclust_prot %>% plot(cex = 0.2, main = "Proteins (Number of samples: AD = 84, PSP = 85)")
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_prot/hclust_ad_psp.png") 


## ----fig.asp = 1--------------------------------------------------------------
# Set colours for heatmap, 25 increments
my_palette <- colorRampPalette(c("blue","white","red"))(n = 25)

# Plot heatmap with heatmap.2
png("./figs_prot/heatmap_ad_psp.png", width = 2000, height = 2000, res = 300)
par(cex.main=0.75) # Shrink title fonts on plot
adpsp_scaled %>%
  # Plot heatmap
  gplots::heatmap.2(.,                     # Tidy, normalised data
          Colv=as.dendrogram(adpsp_hclust),     # Experiments clusters in cols
          Rowv=as.dendrogram(adpsp_hclust_prot),     # Protein clusters in rows
          revC=TRUE,                  # Flip plot to match pheatmap
          density.info="histogram",   # Plot histogram of data and colour key
          trace="none",               # Turn of trace lines from heat map
          col = my_palette,           # Use my colour scheme
          cexRow=0.3, cexCol=0.3,     # Amend row and column label fonts
          xlab = paste0("Number of samples: AD = ", length(ad_cols), ", PSP = ", length(psp_cols))
          )
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_prot/heatmap_ad_psp.png") 


## ----fig.asp = 1--------------------------------------------------------------
venn_list <- list("AD vs Control" =  dat_fc_pv %>%
                    filter(ad_control_pv <= pv_cut) %>%
                    ## filter(ad_control_pv <= pv_cut & abs(fc_control_ad) > fc_cut) %>%
                    pull(PIDshort),
                  "PSP vs Control" = dat_fc_pv %>%
                    filter(psp_control_pv <= pv_cut) %>%
                    ## filter(psp_control_pv <= pv_cut & abs(fc_control_psp) > fc_cut) %>%
                    pull(PIDshort),
                  "AD vs PSP" = dat_fc_pv %>%
                    filter(ad_psp_pv <= pv_cut) %>%
                    ## filter(ad_psp_pv <= pv_cut & abs(fc_ad_psp) > fc_cut) %>%
                    pull(PIDshort)
                  )

# Prevent the output of a log file
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Create a venn diagram object
prot_venn <- venn.diagram(venn_list,
                          NULL,
                          col = "transparent",
                          fill = c("cornflowerblue", "green", "yellow"),
                          alpha = 0.50,
                          cex = 0.8,
                          fontfamily = "sans",
                          fontface = "bold",
                          cat.col = c("darkblue", "darkgreen", "orange"),
                          cat.cex = 0.8,
                          cat.fontfamily = "sans",
                          margin = 0.2,
                          main = paste0("Pairwise significantly DE proteins (adj.p<",pv_cut,") identified in the experiment"),
                          main.pos = c(0.5, 1.05),
                          main.fontfamily = "sans",
                          sub = "Number of samples: AD = 84, PSP = 85, Control = 31",
                          sub.pos = c(0.5, 0.92),
                          sub.cex = 0.8,
                          sub.fontfamily = "sans",
                          print.mode = c("raw","percent"), # Show both numbers and percent
                          )

# Plot the venn diagram using the gridExtra package
png("./figs_prot/venn.png", width = 2000, height = 2000, res = 300)
grid.arrange(gTree(children = prot_venn))
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_prot/venn.png") 


## ----fig.asp = 1--------------------------------------------------------------
venn_list2 <- list(
  "AD-C +" =  dat_fc_pv %>%
    filter(ad_control_pv <= pv_cut & fc_control_ad >= 0) %>%
    pull(PIDshort),
  "AD-C -" =  dat_fc_pv %>%
    filter(ad_control_pv <= pv_cut & fc_control_ad < 0) %>%
    pull(PIDshort),
  "PSP-C +" =  dat_fc_pv %>%
    filter(psp_control_pv <= pv_cut & fc_control_psp >= 0) %>%
    pull(PIDshort),
  "PSP-C -" =  dat_fc_pv %>%
    filter(psp_control_pv <= pv_cut & fc_control_psp < 0) %>%
    pull(PIDshort),
  "AD-PSP +" =  dat_fc_pv %>%
    filter(ad_psp_pv <= pv_cut & fc_ad_psp >= 0) %>%
    pull(PIDshort),
  "AD-PSP -" =  dat_fc_pv %>%
    filter(ad_psp_pv <= pv_cut & fc_ad_psp < 0) %>%
    pull(PIDshort)
) 


## -----------------------------------------------------------------------------
prot_venn_up <- venn.diagram(venn_list2[c(1,3,5)],NULL,
                          col = "transparent",
                          fill = c("cornflowerblue", "green", "yellow"),
                          alpha = 0.50,
                          cex = 0.8,
                          fontfamily = "sans",
                          fontface = "bold",
                          cat.col = c("darkblue", "darkgreen", "orange"),
                          cat.cex = 0.8,
                          cat.fontfamily = "sans",
                          margin = 0.2,
                          main = "Pairwise significantly upregulated proteins (adj.p<0.05) identified in the experiment",
                          main.fontfamily = "sans",
                          sub = "Number of samples: AD = 84, PSP = 85, Control = 31",
                          sub.pos = c(0.5, 0.92),
                          sub.cex = 0.8,
                          sub.fontfamily = "sans",
                          print.mode = c("raw","percent"), # Show both numbers and percent
                          main.pos = c(0.5, 1.05)
                          )

# Plot the venn diagram using the gridExtra package
png("./figs_prot/venn_up.png", width = 2000, height = 2000, res = 300)
grid.arrange(gTree(children = prot_venn_up))
dev.off()
 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_prot/venn_up.png") 


## -----------------------------------------------------------------------------
prot_venn_down <- venn.diagram(venn_list2[c(2,4,6)],NULL,
                          col = "transparent",
                          fill = c("cornflowerblue", "green", "yellow"),
                          alpha = 0.50,
                          cex = 0.8,
                          fontfamily = "sans",
                          fontface = "bold",
                          cat.col = c("darkblue", "darkgreen", "orange"),
                          cat.cex = 0.8,
                          cat.fontfamily = "sans",
                          margin = 0.2,
                          main = "Pairwise significantly downregulated proteins (adj.p<0.05) identified in the experiment",
                          main.pos = c(0.5, 1.05),
                          main.fontfamily = "sans",
                          sub = "Number of samples: AD = 84, PSP = 85, Control = 31",
                          sub.pos = c(0.5, 0.92),
                          sub.cex = 0.8,
                          sub.fontfamily = "sans",
                          print.mode = c("raw","percent") # Show both numbers and percent
                          )

# Plot the venn diagram using the gridExtra package
png("./figs_prot/venn_down.png", width = 2000, height = 2000, res = 300)
grid.arrange(gTree(children = prot_venn_down))
dev.off()
 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_prot/venn_down.png") 


## -----------------------------------------------------------------------------

EnhancedVolcano(dat_fc_pv,
                x = "fc_control_ad",
                y = "ad_control_pv",
                lab = dat_fc_pv$PIDshort,
                title = 'AD versus Control',
                subtitle = "Number of samples: AD = 84, PSP = 85, Control = 31",
                pCutoff = pv_cut,
                FCcutoff = fc_cut,
                pointSize = 1.0,
                labSize = 3.0) +
  ggplot2::scale_y_continuous(limits = c(0, 7)) +
  ggplot2::scale_x_continuous(limits = c(-1, 1))
ggsave("./figs_prot/ad_vs_control_volcano.png")

ad_vs_control_proteins <- dat_fc_pv %>%
  dplyr::select(PIDshort, fc = fc_control_ad, pv = ad_control_pv) %>%
  arrange(-abs(fc)) %>%
  filter(pv <= pv_cut & abs(fc) >= fc_cut)

require(openxlsx)
write.xlsx(ad_vs_control_proteins,
           file = paste0("./results_prot/significantly_de_proteins_pv_", pv_cut, "_log2FC_", fc_cut, ".xlsx"),
           sheetName = "AD vs Control", append=TRUE)

EnhancedVolcano(dat_fc_pv,
                x = "fc_control_psp",
                y = "psp_control_pv",
                lab = dat_fc_pv$PIDshort,
                title = 'PSP versus Control',
                subtitle = "Number of samples: AD = 84, PSP = 85, Control = 31",
                pCutoff = pv_cut,
                FCcutoff = fc_cut,
                pointSize = 1.0,
                labSize = 3.0) +
  ggplot2::scale_y_continuous(limits = c(0, 6)) +
  ggplot2::scale_x_continuous(limits = c(-1, 1))
ggsave("./figs_prot/psp_vs_control_volcano.png")

psp_vs_control_proteins <- dat_fc_pv %>%
  dplyr::select(PIDshort, fc = fc_control_psp, pv = psp_control_pv) %>%
  arrange(-abs(fc)) %>%
  filter(pv <= pv_cut & abs(fc) >= fc_cut)

require(openxlsx)
write.xlsx(psp_vs_control_proteins,
           file = paste0("./results_prot/significantly_de_proteins_pv_", pv_cut, "_log2FC_", fc_cut, ".xlsx"),
           sheetName = "PSP vs Control", append=TRUE)

EnhancedVolcano(dat_fc_pv,
                x = "fc_ad_psp",
                y = "ad_psp_pv",
                lab = dat_fc_pv$PIDshort,
                title = 'PSP versus AD',
                subtitle = "Number of samples: AD = 84, PSP = 85, Control = 31",
                pCutoff = pv_cut,
                FCcutoff = fc_cut,
                pointSize = 1.0,
                labSize = 3.0) +
  ggplot2::scale_y_continuous(limits = c(0, 20)) +
  ggplot2::scale_x_continuous(limits = c(-1, 1))
ggsave("./figs_prot/psp_vs_ad_volcano.png")

psp_vs_ad_proteins <- dat_fc_pv %>%
  dplyr::select(PIDshort, fc = fc_ad_psp, pv = ad_psp_pv) %>%
  arrange(-abs(fc)) %>%
  filter(pv <= pv_cut & abs(fc) >= fc_cut)

require(openxlsx)
## write.xlsx(psp_vs_ad_proteins, "significantly_de_proteins.xlsx", sheetName = "PSP vs AD", append=TRUE)
write.xlsx(psp_vs_ad_proteins,
           file = paste0("./results_prot/significantly_de_proteins_pv_", pv_cut, "_log2FC_", fc_cut, ".xlsx"),
           sheetName = "PSP vs AD", append=TRUE) 


## -----------------------------------------------------------------------------
## On the example of AD vs Control.

## At first i took the significantly DE genes for e.g. AD vs Control.
## Data with FC and p-values ​​has already been collected in dat_fc_pv. I
## had to select data from this table only for AD vs Control, since it
## was necessary to filter out significant genes for this particular
## variant.

ad_vs_control_names <- dat_fc_pv %>%
  ## Filtered out relevant proteins. In this case, with p-vales <= pv_cut and abs(log2FC) >= fc_cut
  filter(abs(fc_control_ad) >= fc_cut & ad_control_pv <= pv_cut) %>%
  ## Sorting the table
  arrange(-abs(fc_control_ad)) %>%
  slice_head(n=25) %>%
  ## The full form of protein names is needed to select them from the
  ## cra.log table.  The short form of protein names is needed for the
  ## graph.
  dplyr::select(ProteinIDs, PIDshort)

## Then, I filtered the required proteins from the cra.log dataframe,
## transposed the resulting table and changed its class from matrix to
## dataframe.
avc <- cra.log[ad_vs_control_names$ProteinIDs, ] %>% t %>% data.frame

## Renamed columns by protein names.
colnames(avc) <- ad_vs_control_names$PIDshort

## Created a "condition" column filled with NAs. It is needed to build boxplots.
avc$condition <- NA

## In the "condition" column, I replaced NAs with the variants names: Control, AD, PSP.
avc[control_cols, ]$condition <- "Control"
avc[ad_cols, ]$condition <- "AD"
avc[psp_cols, ]$condition <- "PSP"
 


## ----fig.show='hide'----------------------------------------------------------
avc %>% filter(condition != "PSP") %>%
  ## For ggplot, you need to convert the data from wide format to long format.
  pivot_longer(-condition, names_to = "protein", values_to = "level") %>%
  ## gglot uses the "condition" column to group boxplots: fill=condition in aes().
  ggplot(aes(protein, level, fill=condition)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_x_discrete(guide = guide_axis(n.dodge = 3)) +
  xlab("") +
  ylab("log2 expression level") +
  ggtitle("AD vs Control top 25 significant proteins (Number of samples: AD = 84, Control = 31)") +
  theme_light() +
  theme(plot.title = element_text(size = 16, face = "bold"))
ggsave("./figs_prot/boxplot_ad_vs_control.png", dpi = 300, scale = 2, width = 6, height = 3, units = "in") 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_prot/boxplot_ad_vs_control.png") 


## ----fig.show='hide'----------------------------------------------------------
psp_vs_control_names <- dat_fc_pv %>%
  filter(psp_control_pv <= pv_cut & abs(fc_control_psp) >= fc_cut) %>%
  arrange(-abs(fc_control_psp)) %>%
  slice_head(n=25) %>%
  dplyr::select(ProteinIDs, PIDshort)

pvc <- cra.log[psp_vs_control_names$ProteinIDs, ] %>% t %>% data.frame

colnames(pvc) <- psp_vs_control_names$PIDshort

pvc$condition <- NA

pvc[control_cols, ]$condition <- "Control"
pvc[ad_cols, ]$condition <- "AD"
pvc[psp_cols, ]$condition <- "PSP"

pvc %>% filter(condition != "AD") %>%
  pivot_longer(-condition, names_to = "protein", values_to = "level") %>%
  mutate(conditions = factor(condition, levels = c("PSP", "Control"))) %>%
  ggplot(aes(protein, level, fill=conditions)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_x_discrete(guide = guide_axis(n.dodge = 3)) +
  xlab("") +
  ylab("log2 expression level") +
  ggtitle("PSP vs Control top 25 significant proteins (Number of samples: PSP = 85, Control = 31)") +
  theme_light() +
  theme(plot.title = element_text(size = 16, face = "bold"))
ggsave("./figs_prot/boxplot_psp_vs_control.png", dpi = 300, scale = 2, width = 6, height = 3, units = "in") 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_prot/boxplot_psp_vs_control.png") 


## ----fig.show='hide'----------------------------------------------------------
psp_vs_ad_names <- dat_fc_pv %>%
  filter(abs(fc_ad_psp) >= fc_cut & ad_psp_pv <= pv_cut ) %>%
  arrange(-abs(fc_ad_psp)) %>%
  slice_head(n=25) %>%
  dplyr::select(ProteinIDs, PIDshort)


pva <- cra.log[psp_vs_ad_names$ProteinIDs, ]  %>% t %>% data.frame

## %>% t %>% data.frame

colnames(pva) <- psp_vs_ad_names$PIDshort

pva$condition <- NA

pva[control_cols, ]$condition <- "Control"
pva[ad_cols, ]$condition <- "AD"
pva[psp_cols, ]$condition <- "PSP"


pva %>% filter(condition != "Control") %>%
  pivot_longer(-condition, names_to = "protein", values_to = "level") %>%
  ggplot(aes(protein, level, fill=condition)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_x_discrete(guide = guide_axis(n.dodge = 3)) +
  xlab("") +
  ylab("log2 expression level") +
  ggtitle("PSP vs AD top 25 significant proteins (Number of samples: AD = 84, PSP = 85)") +
  theme_light() +
  theme(plot.title = element_text(size = 16, face = "bold"))
ggsave("./figs_prot/boxplot_ad_vs_psp.png", dpi = 300, scale = 2, width = 6, height = 3, units = "in") 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_prot/boxplot_ad_vs_psp.png") 


## ----cache=TRUE---------------------------------------------------------------
transcripts <- read.delim("./data_transcr/MayoRNAseq_RNAseq_TCX_transcriptCounts.tsv",
                          stringsAsFactors = FALSE) %>%
  `rownames<-`(.$ensembl_id)

ensembl_id <- data.frame(transcripts[1])

transcripts <- transcripts[-1] 


## ----cache=TRUE---------------------------------------------------------------
## transcripts_norm %>% names %>% sub("^.+_([0-9]+)_.+", "\\1", .)
names(transcripts) <- transcripts %>%
  names %>%
  sub("^X([0-9]+_.+)", "\\1", .)

id_keys_raw <- read.csv("./data_transcr/Mayo_Proteomics_ID_key.csv",
                        strip.white = TRUE,
                        header = TRUE,
                        stringsAsFactors = FALSE) %>%
  na.omit %>%
  dplyr::select(1:2)

id_keys <- read.csv("./data_transcr/Mayo_Proteomics_TC_traits.csv",
                    strip.white = TRUE,
                    header = TRUE,
                    stringsAsFactors = FALSE) %>%
  slice_head(n=200) %>%
  dplyr::select(c(1,4)) %>%
  cbind.data.frame(., id_keys_raw) %>%
  dplyr::select(c(1,4,2)) %>%
  filter(Samples_Simple != "b1_091_20b")

rm(id_keys_raw) 


## -----------------------------------------------------------------------------
samples_in_common <- intersect(names(transcripts),
                               id_keys$RNA_SampleID) 


## -----------------------------------------------------------------------------
id_keys <- id_keys %>%
  filter(RNA_SampleID %in% samples_in_common) 


## ----cache=TRUE---------------------------------------------------------------
transcripts <- transcripts %>%
  dplyr::select(all_of(id_keys$RNA_SampleID)) 


## -----------------------------------------------------------------------------
control_cols <- id_keys %>%
  mutate(n = rownames(.) %>% as.numeric) %>%
  filter(Diagnosis == "Control") %>%
  .$n

ad_cols <- id_keys %>%
  mutate(n = rownames(.) %>% as.numeric) %>%
  filter(Diagnosis == "AD") %>%
  .$n

psp_cols <- id_keys %>%
  mutate(n = rownames(.) %>% as.numeric) %>%
  filter(Diagnosis == "PSP") %>%
  .$n 


## ----cache=TRUE---------------------------------------------------------------
# Obtain CPMs
transcripts_cpm <- cpm(transcripts)
# Have a look at the output
head(transcripts_cpm)[,1:5] 


## ----cache=TRUE---------------------------------------------------------------
# Which values in myCPM are greater than 0.5?
thresh <- transcripts_cpm > 1 
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)[,1:5] 


## ----cache=TRUE---------------------------------------------------------------
# Summary of how many TRUEs there are in each row
# There are 11433 genes that have TRUEs in all 12 samples.
table(rowSums(thresh)) 


## ----cache=TRUE---------------------------------------------------------------
# we would like to keep genes that have at least 31 TRUES in each row of thresh
keep <- rowSums(thresh) >= 31
# Subset the rows of countdata to keep the more highly expressed genes
transcripts_keep <- transcripts[keep,]
summary(keep) 


## -----------------------------------------------------------------------------
dim(transcripts_keep) 


## -----------------------------------------------------------------------------
# Let's have a look and see whether our threshold of 1 does indeed correspond to a count of about 30-40
# We will look at the first sample
## plot(transcripts_cpm[,1],transcripts[,1],
##      ylim=c(0,50),xlim=c(0,3))
## abline(v=1) 


## ----cache=TRUE---------------------------------------------------------------
dgeObj <- DGEList(transcripts_keep)
# have a look at dgeObj
dgeObj[,1:5] 


## ----cache=TRUE---------------------------------------------------------------
dgeObj$samples$lib.size 


## ----cache=TRUE---------------------------------------------------------------
png("./figs_transcr/barplot_sample_sizes.png", width = 2000, height = 2000, res = 300)
# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(dgeObj$samples$lib.size, names=colnames(dgeObj), las=2, cex.names = .5)
# Add a title to the plot
## title("Barplot of library sizes")
title(paste("Barplot of library sizes\n", "Number of samples: AD = 84, PSP = 85, Control = 31"))
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_transcr/barplot_sample_sizes.png") 


## -----------------------------------------------------------------------------
# Get log2 counts per million
logcounts <- cpm(dgeObj, log = TRUE)
# Check distributions of samples using boxplot
png("./figs_transcr/boxplot_samples_unnorm.png", width = 4000, height = 2000, res = 300)
boxplot(logcounts, xlab="", ylab="Log2 counts per million", las=2, cex.axis=.4, outcex=.1)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised). Number of samples: AD = 82, PSP = 84, Control = 31")
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_transcr/boxplot_samples_unnorm.png") 


## -----------------------------------------------------------------------------
levels(id_keys$Diagnosis %>% as.factor) 


## -----------------------------------------------------------------------------
col.condition <- c("red", "green", "purple")[id_keys$Diagnosis %>% as.factor] 


## ----cache=TRUE---------------------------------------------------------------
png("./figs_transcr/mds_plot.png", width = 2000, height = 2000, res = 300)
plotMDS(dgeObj, col=col.condition)
legend("topleft", fill=c("red", "green", "purple"), legend = levels(id_keys$Diagnosis %>% as.factor))
# Add a title
title("MDS plot. Number of samples: AD = 82, PSP = 84, Control = 31")
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_transcr/mds_plot.png") 


## -----------------------------------------------------------------------------
# Apply normalisation to DGEList object
dgeObj <- calcNormFactors(dgeObj) 


## -----------------------------------------------------------------------------
dgeObj$samples %>% head 


## -----------------------------------------------------------------------------
png("./figs_transcr/mean_difference.png", width = 2000, height = 2000, res = 300)
par(mfrow=c(1,2))
plotMD(logcounts, column = 1)
abline(h=0,col="grey")
plotMD(logcounts,column = 2)
abline(h=0,col="grey")
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_transcr/mean_difference.png") 


## -----------------------------------------------------------------------------
png("./figs_transcr/mean_difference_2.png", width = 2000, height = 2000, res = 300)
par(mfrow=c(1,2))
plotMD(dgeObj,column = 1)
abline(h=0,col="grey")
plotMD(dgeObj,column = 2)
abline(h=0,col="grey")
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_transcr/mean_difference_2.png") 


## -----------------------------------------------------------------------------
## save(group, dgeObj,sampleinfo,file="./Robjects/preprocessing.Rdata") 


## -----------------------------------------------------------------------------
diagnosis <- as.character(id_keys$Diagnosis)
design <- model.matrix( ~ diagnosis + 0) 


## -----------------------------------------------------------------------------
## plotMDS(dgeObj, labels=diagnosis, cex=0.75, xlim=c(-4, 5)) 


## -----------------------------------------------------------------------------
dgeObj <- estimateCommonDisp(dgeObj) 


## -----------------------------------------------------------------------------
dgeObj <- estimateGLMTrendedDisp(dgeObj)
dgeObj <- estimateTagwiseDisp(dgeObj) 


## -----------------------------------------------------------------------------
png("./figs_transcr/est_dispersion.png", width = 2000, height = 2000, res = 300)
plotBCV(dgeObj)
mytitle = "Number of samples: AD = 82, PSP = 84, Control = 31"
mtext(side=3, line=2, at=-0.01, adj=0, cex=1, mytitle)
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_transcr/est_dispersion.png") 


## -----------------------------------------------------------------------------
# Fit the linear model
fit <- glmFit(dgeObj, design)
names(fit) 


## -----------------------------------------------------------------------------
head(coef(fit)) 


## -----------------------------------------------------------------------------
AvsC <- makeContrasts(diagnosisAD - diagnosisControl, levels=design) 


## -----------------------------------------------------------------------------
lrt.AvsC <- glmLRT(fit, contrast=AvsC)
## topTags(lrt.AvsC, sort.by = "logFC", p.value = 0.05, n = 61) 


## -----------------------------------------------------------------------------
## gene_ids_tr <- lrt.AvsC %>%
##   rownames %>%
##   gconvert(organism="hsapiens", target="ENSG", filter_na = FALSE) %>%
##   dplyr::select("GeneID" = name)

## saveRDS(gene_ids_tr, "./results_transcr/gene_ids_tr.rds")

gene_ids_tr <- readRDS("./results_transcr/gene_ids_tr.rds") 


## -----------------------------------------------------------------------------
## Set log-fold cutoff for transcripts
fc_cut_tr <- 2

ad_vs_control_top_tr  <- topTags(lrt.AvsC, sort.by = "logFC", p.value = 0.05, n = 100) %>%
  data.frame %>%
  filter(abs(logFC) >= fc_cut_tr) %>%
  rownames_to_column(var = "TranscriptIDs")

AD_vs_Control_GeneIDs <- ad_vs_control_top_tr$TranscriptIDs %>%
  noquote %>%
  gconvert(organism="hsapiens", target="ENSG") %>%
  dplyr::select("GeneID" = name)
## write.csv(ad_vs_control_prot, "./results/prot_ad_vs_control.csv")

write.csv(data.frame(AD_vs_Control_GeneIDs, ad_vs_control_top_tr),
          "./results_transcr/ad_vs_control_tr.csv")

knitr::kable(data.frame(AD_vs_Control_GeneIDs, ad_vs_control_top_tr),
             caption = "Top DE transcripts: AD vs Control") 


## -----------------------------------------------------------------------------
all_samples_ad_vs_control_tr <- cpm(dgeObj, normalized.lib.sizes = TRUE) %>%
  data.frame() %>%
  filter(rownames(.) %in% ad_vs_control_top_tr$TranscriptIDs) %>%
  dplyr::select(all_of(ad_cols)) %>%
  `rownames<-`(AD_vs_Control_GeneIDs$GeneID %>% make.unique) 


## -----------------------------------------------------------------------------
PvsC <- makeContrasts(diagnosisPSP - diagnosisControl, levels=design) 


## -----------------------------------------------------------------------------
lrt.PvsC <- glmLRT(fit, contrast=PvsC)
## topTags(lrt.PvsC, sort.by = "logFC", p.value = 0.05, n = 25) 


## -----------------------------------------------------------------------------
## fold-change cutoff here is 1.9 to make the numger of transcripts
## close to the number of proteins
psp_vs_control_top_tr  <- topTags(lrt.PvsC, sort.by = "logFC", p.value = 0.05, n = 200) %>%
  data.frame %>%
  filter(abs(logFC) >= 1.9) %>%
  rownames_to_column(var = "TranscriptIDs")

PSP_vs_Control_GeneIDs <- psp_vs_control_top_tr$TranscriptIDs %>%
  noquote %>%
  gconvert(organism="hsapiens", target="ENSG") %>%
  dplyr::select("GeneID" = name)

write.csv(data.frame(PSP_vs_Control_GeneIDs, psp_vs_control_top_tr),
          "./results_transcr/psp_vs_control_tr.csv")

knitr::kable(data.frame(PSP_vs_Control_GeneIDs, psp_vs_control_top_tr),
             caption = "Top DE transcripts: PSP vs Control") 


## -----------------------------------------------------------------------------
all_samples_psp_vs_control_tr <- cpm(dgeObj, normalized.lib.sizes = TRUE) %>%
  data.frame() %>%
  filter(rownames(.) %in% psp_vs_control_top_tr$TranscriptIDs) %>%
  dplyr::select(all_of(psp_cols)) %>%
  `rownames<-`(PSP_vs_Control_GeneIDs$GeneID %>% make.unique) 


## -----------------------------------------------------------------------------
AvsP <- makeContrasts(diagnosisAD - diagnosisPSP, levels=design) 


## -----------------------------------------------------------------------------
lrt.AvsP <- glmLRT(fit, contrast=AvsP)
## topTags(lrt.AvsP, sort.by = "logFC", p.value = 0.05, n = 25) 


## -----------------------------------------------------------------------------
## fold-change cutoff here is 1 to make the numger of transcripts
## close to the number of proteins
ad_vs_psp_top_tr  <- topTags(lrt.AvsP, sort.by = "logFC", p.value = 0.05, n = 500) %>%
  data.frame %>%
  filter(abs(logFC) >= 1) %>%
  rownames_to_column(var = "TranscriptIDs")

AD_vs_PSP_GeneIDs <- ad_vs_psp_top_tr$TranscriptIDs %>%
  noquote %>%
  gconvert(organism="hsapiens", target="ENSG", filter_na = FALSE) %>%
  dplyr::select("GeneID" = name)
## write.csv(ad_vs_psp_prot, "./results/prot_ad_vs_psp.csv")

write.csv(data.frame(AD_vs_PSP_GeneIDs, ad_vs_psp_top_tr),
          "./results_transcr/ad_vs_psp_tr.csv")

knitr::kable(data.frame(AD_vs_PSP_GeneIDs, ad_vs_psp_top_tr),
             caption = "Top DE transcripts: AD vs PSP") 


## -----------------------------------------------------------------------------
all_samples_ad_vs_psp_tr <- cpm(dgeObj, normalized.lib.sizes = TRUE) %>%
  data.frame() %>%
  filter(rownames(.) %in% ad_vs_psp_top_tr$TranscriptIDs) %>%
  dplyr::select(all_of(ad_cols)) %>%
  `rownames<-`(AD_vs_PSP_GeneIDs$GeneID %>% make.unique) 


## -----------------------------------------------------------------------------
write.xlsx(all_samples_ad_vs_control_tr, "./results_transcr/transcriptsdata.xlsx",
           sheetName = "all_samples_ad_vs_control_tr", append = TRUE)
write.xlsx(all_samples_ad_vs_psp_tr, "./results_transcr/transcriptsdata.xlsx",
           sheetName = "all_samples_ad_vs_psp_tr", append = TRUE)
write.xlsx(all_samples_psp_vs_control_tr, "./results_transcr/transcriptsdata.xlsx",
           sheetName = "all_samples_psp_vs_control_tr", append = TRUE)
write.xlsx(ad_vs_control_top_tr, "./results_transcr/transcriptsdata.xlsx",
           sheetName = "ad_vs_control_top_tr", append = TRUE)
write.xlsx(psp_vs_control_top_tr, "./results_transcr/transcriptsdata.xlsx",
           sheetName = "psp_vs_control_top_tr", append = TRUE)
write.xlsx(ad_vs_psp_top_tr, "./results_transcr/transcriptsdata.xlsx",
           sheetName = "ad_vs_psp_top_tr", append = TRUE) 


## -----------------------------------------------------------------------------
all_samples_ad_and_control_tr <- cpm(dgeObj, normalized.lib.sizes = TRUE) %>%
  data.frame() %>%
  filter(rownames(.) %in% ad_vs_control_top_tr$TranscriptIDs) %>%
  dplyr::select(all_of(c(ad_cols, control_cols))) %>%
  `colnames<-`(c(make.unique(rep("AD", length(ad_cols))),
                 make.unique(rep("Control", length(control_cols))))) %>%
  `rownames<-`(AD_vs_Control_GeneIDs$GeneID %>% make.unique) 


## -----------------------------------------------------------------------------
## data.matrix() leaves numeric data as is.
ad_tr_matrix <- data.matrix(all_samples_ad_and_control_tr)

## rownames(ad_tr_matrix) <- dat_filt_ad$PIDshort
## ad_tr_matrix <- ad_tr_matrix[,-1] 


## -----------------------------------------------------------------------------
ad_tr_scaled <- ad_tr_matrix %>% t %>% scale %>% t 


## -----------------------------------------------------------------------------
ad_dist_sampl <- ad_tr_scaled %>% t %>%
  dist(., method = "euclidean", diag = FALSE, upper = FALSE) 


## -----------------------------------------------------------------------------
ad_dist_tr <- ad_tr_scaled %>% 
  dist(., method = "euclidean", diag = FALSE, upper = FALSE) 


## ----fig.asp = 1.5------------------------------------------------------------
ad_hclust_sampl <- hclust(ad_dist_sampl, method = "ward.D2", members = NULL)
ad_hclust_tr <- hclust(ad_dist_tr, method = "ward.D2", members = NULL)


png("./figs_transcr/hclust_ad.png", width = 2500, height = 3000, res = 300)
par(mfrow = c(2,1))
ad_hclust_sampl %>% plot(cex = 0.4, main = "Conditions. Number of samples: AD = 82, Control = 31")
ad_hclust_tr %>% plot(cex = 0.6, main = "Proteins. Number of samples: AD = 82, Control = 31")
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_transcr/hclust_ad.png") 


## -----------------------------------------------------------------------------
# Set colours for heatmap, 25 increments
my_palette <- colorRampPalette(c("blue","white","red"))(n = 25)

png("./figs_transcr/heatmap_ad.png", width = 2000, height = 2000, res = 300)
# Plot heatmap with heatmap.2
par(cex.main=0.75) # Shrink title fonts on plot
ad_tr_scaled %>%
  # Plot heatmap
  gplots::heatmap.2(.,                     # Tidy, normalised data
                    Colv=as.dendrogram(ad_hclust_sampl),     # Experiments clusters in cols
                    Rowv=as.dendrogram(ad_hclust_tr),     # Protein clusters in rows
                    revC=TRUE,                  # Flip plot to match pheatmap
                    density.info="histogram",   # Plot histogram of data and colour key
                    trace="none",               # Turn of trace lines from heat map
                    col = my_palette,           # Use my colour scheme
                    cexRow=0.3, cexCol=0.2,     # Amend row and column label fonts
                    xlab = paste0("Number of samples: AD = ", length(ad_cols), ", Control = ", length(control_cols))
                    )
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_transcr/heatmap_ad.png") 


## -----------------------------------------------------------------------------
all_samples_psp_and_control_tr <- cpm(dgeObj, normalized.lib.sizes = TRUE) %>%
  data.frame() %>%
  filter(rownames(.) %in% psp_vs_control_top_tr$TranscriptIDs) %>%
  dplyr::select(all_of(c(psp_cols, control_cols))) %>%
  `colnames<-`(c(make.unique(rep("PSP", length(ad_cols))),
                 make.unique(rep("Control", length(control_cols))))) %>%
  `rownames<-`(PSP_vs_Control_GeneIDs$GeneID %>% make.unique) 


## -----------------------------------------------------------------------------
## data.matrix() leaves numeric data as is.
psp_tr_matrix <- data.matrix(all_samples_psp_and_control_tr) 


## -----------------------------------------------------------------------------
psp_tr_scaled <- psp_tr_matrix %>% t %>% scale %>% t 


## -----------------------------------------------------------------------------
psp_dist_sampl <- psp_tr_scaled %>% t %>%
  dist(., method = "euclidean", diag = FALSE, upper = FALSE) 


## -----------------------------------------------------------------------------
psp_dist_tr <- psp_tr_scaled %>% 
  dist(., method = "euclidean", diag = FALSE, upper = FALSE) 


## ----fig.asp = 1.5------------------------------------------------------------
psp_hclust_sampl <- hclust(psp_dist_sampl, method = "ward.D2", members = NULL)
psp_hclust_tr <- hclust(psp_dist_tr, method = "ward.D2", members = NULL)


png("./figs_transcr/hclust_psp.png", width = 2500, height = 3000, res = 300)
par(mfrow = c(2,1))
psp_hclust_sampl %>% plot(cex = 0.4, main = "Conditions. Number of samples: PSP = 84, Control = 31")
psp_hclust_tr %>% plot(cex = 0.6, main = "Proteins. Number of samples: PSP = 84, Control = 31")
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_transcr/hclust_psp.png") 


## ----fig.asp = 1--------------------------------------------------------------
# Set colours for heatmap, 25 increments
my_palette <- colorRampPalette(c("blue","white","red"))(n = 25)

png("./figs_transcr/heatmap_psp.png", width = 2000, height = 2000, res = 300)
# Plot heatmap with heatmap.2
par(cex.main=0.75) # Shrink title fonts on plot
psp_tr_scaled %>%
  # Plot heatmap
  gplots::heatmap.2(.,                     # Tidy, normalised data
                    Colv=as.dendrogram(psp_hclust_sampl),     # Experiments clusters in cols
                    Rowv=as.dendrogram(psp_hclust_tr),     # Protein clusters in rows
                    revC=TRUE,                  # Flip plot to match pheatmap
                    density.info="histogram",   # Plot histogram of data and colour key
                    trace="none",               # Turn of trace lines from heat map
                    col = my_palette,           # Use my colour scheme
                    cexRow=0.3, cexCol=0.2,     # Amend row and column label fonts
                    xlab = paste0("Number of samples: PSP = ", length(psp_cols), ", Control = ", length(control_cols))
                    )
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_transcr/heatmap_psp.png") 


## -----------------------------------------------------------------------------
all_samples_ad_and_psp_tr <- cpm(dgeObj, normalized.lib.sizes = TRUE) %>%
  data.frame() %>%
  filter(rownames(.) %in% ad_vs_psp_top_tr$TranscriptIDs) %>%
  dplyr::select(all_of(c(ad_cols, psp_cols))) %>%
  `colnames<-`(c(make.unique(rep("AD", length(ad_cols))),
                 make.unique(rep("PSP", length(control_cols))))) %>%
  `rownames<-`(AD_vs_PSP_GeneIDs$GeneID %>% make.unique) 


## -----------------------------------------------------------------------------
## data.matrix() leaves numeric data as is.
ad_psp_tr_matrix <- data.matrix(all_samples_ad_and_psp_tr) 


## -----------------------------------------------------------------------------
ad_psp_tr_scaled <- ad_psp_tr_matrix %>% t %>% scale %>% t 


## -----------------------------------------------------------------------------
ad_psp_dist_sampl <- ad_psp_tr_scaled %>% t %>%
  dist(., method = "euclidean", diag = FALSE, upper = FALSE) 


## -----------------------------------------------------------------------------
ad_psp_dist_tr <- ad_psp_tr_scaled %>% 
  dist(., method = "euclidean", diag = FALSE, upper = FALSE) 


## ----fig.asp = 1.5------------------------------------------------------------
ad_psp_hclust_sampl <- hclust(ad_psp_dist_sampl, method = "ward.D2", members = NULL)
ad_psp_hclust_tr <- hclust(ad_psp_dist_tr, method = "ward.D2", members = NULL)


png("./figs_transcr/hclust_ad_psp.png", width = 2500, height = 3000, res = 300)
par(mfrow = c(2,1))
ad_psp_hclust_sampl %>% plot(cex = 0.4, main = "Conditions. Number of samples: AD = 82, PSP = 84")
ad_psp_hclust_tr %>% plot(cex = 0.6, main = "Proteins. Number of samples: AD = 82, PSP = 84")
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_transcr/hclust_ad_psp.png") 


## -----------------------------------------------------------------------------
# Set colours for heatmap, 25 increments
my_palette <- colorRampPalette(c("blue","white","red"))(n = 25)

png("./figs_transcr/heatmap_ad_psp.png", width = 2000, height = 2000, res = 300)
# Plot heatmap with heatmap.2
par(cex.main=0.75) # Shrink title fonts on plot
ad_psp_tr_scaled %>%
  # Plot heatmap
  gplots::heatmap.2(.,                     # Tidy, normalised data
                    Colv=as.dendrogram(ad_psp_hclust_sampl),     # Experiments clusters in cols
                    Rowv=as.dendrogram(ad_psp_hclust_tr),     # Protein clusters in rows
                    revC=TRUE,                  # Flip plot to match pheatmap
                    density.info="histogram",   # Plot histogram of data and colour key
                    trace="none",               # Turn of trace lines from heat map
                    col = my_palette,           # Use my colour scheme
                    cexRow=0.3, cexCol=0.2,     # Amend row and column label fonts
                    xlab = paste0("Number of samples: AD = ", length(ad_cols), ", PSP = ", length(psp_cols))
                    )
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_transcr/heatmap_ad_psp.png") 


## ----fig.asp = 1--------------------------------------------------------------
## venn_list_tr <- list("AD vs Control" =  topTags(lrt.AvsC, sort.by = "logFC",
##                                                 p.value = 0.05, n = nrow(lrt.AvsC)) %>%
##                        rownames() %>%
##                        noquote() %>%
##                        gconvert(organism="hsapiens", target="ENSG") %>%
##                        pull(name),
##                      "PSP vs Control" = topTags(lrt.PvsC, sort.by = "logFC",
##                                                 p.value = 0.05, n = nrow(lrt.PvsC)) %>%
##                        rownames() %>%
##                        noquote() %>%
##                        gconvert(organism="hsapiens", target="ENSG") %>%
##                        pull(name),
##                      "AD vs PSP" = topTags(lrt.AvsP, sort.by = "logFC",
##                                            p.value = 0.05, n = nrow(lrt.AvsP)) %>%
##                        rownames() %>%
##                        noquote() %>%
##                        gconvert(organism="hsapiens", target="ENSG") %>%
##                        pull(name))

## saveRDS(venn_list_tr, "./results_transcr/venn_list_tr.rds")

venn_list_tr <- readRDS("./results_transcr/venn_list_tr.rds")

# Prevent the output of a log file
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Create a venn diagram object
tr_venn <- venn.diagram(venn_list_tr,
                          NULL,
                          col = "transparent",
                          fill = c("cornflowerblue", "green", "yellow"),
                          alpha = 0.50,
                          cex = 0.8,
                          fontfamily = "sans",
                          fontface = "bold",
                          cat.col = c("darkblue", "darkgreen", "orange"),
                          cat.cex = 0.8,
                          cat.fontfamily = "sans",
                          margin = 0.2,
                          main = "Pairwise significantly DE transctipts identified in the experiment",
                          main.pos = c(0.5, 1.05),
                          main.fontfamily = "sans",
                          sub = "Number of samples: AD = 82, PSP = 84, Control = 31",
                          sub.pos = c(0.5, 0.92),
                          sub.cex = 0.8,
                          sub.fontfamily = "sans",
                          print.mode = c("raw","percent") # Show both numbers and percent
                          )

# Plot the venn diagram using the gridExtra package
png("./figs_transcr/venn.png", width = 2000, height = 2000, res = 300)
grid.arrange(gTree(children = tr_venn))
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_transcr/venn.png") 


## ----fig.asp = 1--------------------------------------------------------------
## venn_list_tr_2 <- list(
##   "AD-C +"  =  topTags(lrt.AvsC, sort.by = "logFC",
##                        p.value = 0.05, n = nrow(lrt.AvsC)) %>%
##     data.frame() %>%
##     filter(logFC > 0) %>%
##     rownames() %>%
##     noquote() %>%
##     gconvert(organism="hsapiens", target="ENSG") %>%
##     pull(name),
##   "AD-C -"  =  topTags(lrt.AvsC, sort.by = "logFC",
##                        p.value = 0.05, n = nrow(lrt.AvsC)) %>%
##     data.frame() %>%
##     filter(logFC < 0) %>%
##     rownames() %>%
##     noquote() %>%
##     gconvert(organism="hsapiens", target="ENSG") %>%
##     pull(name),
##   "P-C +" = topTags(lrt.PvsC, sort.by = "logFC",
##                     p.value = 0.05, n = nrow(lrt.PvsC)) %>%
##     data.frame() %>%
##     filter(logFC > 0) %>%
##     rownames() %>%
##     noquote() %>%
##     gconvert(organism="hsapiens", target="ENSG") %>%
##     pull(name),
##   "P-C -" = topTags(lrt.PvsC, sort.by = "logFC",
##                     p.value = 0.05, n = nrow(lrt.PvsC)) %>%
##     data.frame() %>%
##     filter(logFC < 0) %>%
##     rownames() %>%
##     noquote() %>%
##     gconvert(organism="hsapiens", target="ENSG") %>%
##     pull(name),
##   "AD-PSP +" = topTags(lrt.AvsP, sort.by = "logFC",
##                        p.value = 0.05, n = nrow(lrt.AvsP)) %>%
##     data.frame() %>%
##     filter(logFC > 0) %>%
##     rownames() %>%
##     noquote() %>%
##     gconvert(organism="hsapiens", target="ENSG") %>%
##     pull(name),
##   "AD-PSP -" = topTags(lrt.AvsP, sort.by = "logFC",
##                        p.value = 0.05, n = nrow(lrt.AvsP)) %>%
##     data.frame() %>%
##     filter(logFC < 0) %>%
##     rownames() %>%
##     noquote() %>%
##     gconvert(organism="hsapiens", target="ENSG") %>%
##     pull(name)
## )

## saveRDS(venn_list_tr_2, "./results_transcr/venn_list_tr_2.rds")

venn_list_tr_2 <- readRDS("./results_transcr/venn_list_tr_2.rds") 


## -----------------------------------------------------------------------------
tr_venn_up <- venn.diagram(venn_list_tr_2[c(1,3,5)],NULL,
                          col = "transparent",
                          fill = c("cornflowerblue", "green", "yellow"),
                          alpha = 0.50,
                          cex = 0.8,
                          fontfamily = "sans",
                          fontface = "bold",
                          cat.col = c("darkblue", "darkgreen", "orange"),
                          cat.cex = 0.8,
                          cat.fontfamily = "sans",
                          margin = 0.2,
                          main = "Pairwise significantly upregulated transcripts (adj.p<0.05) identified in the experiment",
                          main.pos = c(0.5, 1.05),
                          main.fontfamily = "sans",
                          sub = "Number of samples: AD = 82, PSP = 84, Control = 31",
                          sub.pos = c(0.5, 0.92),
                          sub.cex = 0.8,
                          sub.fontfamily = "sans",
                          print.mode = c("raw","percent") # Show both numbers and percent
                          )

# Plot the venn diagram using the gridExtra package
png("./figs_transcr/venn_up.png", width = 2000, height = 2000, res = 300)
grid.arrange(gTree(children = tr_venn_up))
dev.off()
 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_transcr/venn_up.png") 


## -----------------------------------------------------------------------------
tr_venn_down <- venn.diagram(venn_list_tr_2[c(2,4,6)],NULL,
                          col = "transparent",
                          fill = c("cornflowerblue", "green", "yellow"),
                          alpha = 0.50,
                          cex = 0.8,
                          fontfamily = "sans",
                          fontface = "bold",
                          cat.col = c("darkblue", "darkgreen", "orange"),
                          cat.cex = 0.8,
                          cat.fontfamily = "sans",
                          margin = 0.2,
                          main = "Pairwise significantly downregulated transcripts (adj.p<0.05) identified in the experiment",
                          main.pos = c(0.5, 1.05),
                          main.fontfamily = "sans",
                          sub = "Number of samples: AD = 82, PSP = 84, Control = 31",
                          sub.pos = c(0.5, 0.92),
                          sub.cex = 0.8,
                          sub.fontfamily = "sans",
                          print.mode = c("raw","percent") # Show both numbers and percent
                          )

# Plot the venn diagram using the gridExtra package
png("./figs_transcr/venn_down.png", width = 2000, height = 2000, res = 300)
grid.arrange(gTree(children = tr_venn_down))
dev.off()
 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_transcr/venn_down.png") 


## -----------------------------------------------------------------------------

EnhancedVolcano(lrt.AvsC %>% data.frame,
                x = "logFC",
                y = "PValue",
                lab = gene_ids_tr$GeneID,
                title = 'AD versus Control',
                subtitle = "Number of samples: AD = 82, PSP = 84, Control = 31",
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 0.5,
                axisLabSize = 10,
                labSize = 2.0)
ggsave("./figs_transcr/ad_vs_control_volcano.png") 


## -----------------------------------------------------------------------------

EnhancedVolcano(lrt.PvsC %>% data.frame,
                x = "logFC",
                y = "PValue",
                lab = gene_ids_tr$GeneID,
                title = 'PSP versus Control',
                subtitle = "Number of samples: AD = 82, PSP = 84, Control = 31",
                pCutoff = 0.05,
                FCcutoff = 1.9,
                pointSize = 0.5,
                axisLabSize = 10,
                labSize = 2.0)
ggsave("./figs_transcr/psp_vs_control_volcano.png")
 


## -----------------------------------------------------------------------------
EnhancedVolcano(lrt.AvsP %>% data.frame,
                x = "logFC",
                y = "PValue",
                lab = gene_ids_tr$GeneID,
                title = 'AD versus PSP',
                subtitle = "Number of samples: AD = 82, PSP = 84, Control = 31",
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 0.5,
                axisLabSize = 10,
                labSize = 2.0)
ggsave("./figs_transcr/ad_vs_psp_volcano.png") 


## ----fig.show='hide'----------------------------------------------------------
cpm(dgeObj, log = TRUE, normalized.lib.sizes = TRUE) %>%
  data.frame() %>%
  filter(rownames(.) %in% ad_vs_control_top_tr$TranscriptIDs) %>%
  dplyr::select(all_of(c(ad_cols, control_cols, psp_cols))) %>%
  `rownames<-`(AD_vs_Control_GeneIDs$GeneID %>% make.unique) %>%
  t %>% data.frame %>%
  dplyr::select(1:25) %>%
  mutate(condition = c(rep("AD", length(ad_cols)),
                       rep("Control", length(control_cols)),
                       rep("PSP", length(psp_cols)))) %>%
  filter(condition != "PSP") %>%
  ## For ggplot, you need to convert the data from wide format to long format.
  pivot_longer(-condition, names_to = "protein", values_to = "level") %>%
  ## gglot uses the "condition" column to group boxplots: fill=condition in aes().
  ggplot(aes(protein, level, fill=condition)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_x_discrete(guide = guide_axis(n.dodge = 3)) +
  xlab("") +
  ylab("log2 expression level") +
  ggtitle("AD vs Control top 25 significant genes. Number of samples: AD = 82, Control = 31") +
  theme_light() +
  theme(plot.title = element_text(size = 16, face = "bold"))
ggsave("./figs_transcr/boxplot_ad_vs_control.png", dpi = 300, scale = 2, width = 6, height = 3, units = "in") 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_transcr/boxplot_ad_vs_control.png") 


## ----fig.show='hide'----------------------------------------------------------
cpm(dgeObj, log = TRUE, normalized.lib.sizes = TRUE) %>%
  data.frame() %>%
  filter(rownames(.) %in% psp_vs_control_top_tr$TranscriptIDs) %>%
  dplyr::select(all_of(c(ad_cols, control_cols, psp_cols))) %>%
  `rownames<-`(PSP_vs_Control_GeneIDs$GeneID %>% make.unique) %>%
  t %>% data.frame %>%
  dplyr::select(1:25) %>%
  mutate(condition = c(rep("AD", length(ad_cols)),
                       rep("Control", length(control_cols)),
                       rep("PSP", length(psp_cols)))) %>%
  filter(condition != "AD") %>%
  ## For ggplot, you need to convert the data from wide format to long format.
  pivot_longer(-condition, names_to = "protein", values_to = "level") %>%
  ## gglot uses the "condition" column to group boxplots: fill=condition in aes().
  mutate(conditions = factor(condition, levels = c("PSP", "Control"))) %>%
  ggplot(aes(protein, level, fill=conditions)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_x_discrete(guide = guide_axis(n.dodge = 3)) +
  xlab("") +
  ylab("log2 expression level") +
  ggtitle("PSP vs Control top 25 significant genes. Number of samples: PSP = 84, Control = 31") +
  theme_light() +
  theme(plot.title = element_text(size = 16, face = "bold"))
ggsave("./figs_transcr/boxplot_psp_vs_control.png", dpi = 300, scale = 2, width = 6, height = 3, units = "in") 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_transcr/boxplot_psp_vs_control.png") 


## ----fig.show='hide'----------------------------------------------------------
cpm(dgeObj, log = TRUE, normalized.lib.sizes = TRUE) %>%
  data.frame() %>%
  filter(rownames(.) %in% ad_vs_psp_top_tr$TranscriptIDs) %>%
  dplyr::select(all_of(c(ad_cols, control_cols, psp_cols))) %>%
  `rownames<-`(AD_vs_PSP_GeneIDs$GeneID %>% make.unique) %>%
  t %>% data.frame %>%
  dplyr::select(1:25) %>%
  mutate(condition = c(rep("AD", length(ad_cols)),
                       rep("Control", length(control_cols)),
                       rep("PSP", length(psp_cols)))) %>%
  filter(condition != "Control") %>%
  ## For ggplot, you need to convert the data from wide format to long format.
  pivot_longer(-condition, names_to = "protein", values_to = "level") %>%
  ## gglot uses the "condition" column to group boxplots: fill=condition in aes().
  ggplot(aes(protein, level, fill=condition)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_x_discrete(guide = guide_axis(n.dodge = 3)) +
  xlab("") +
  ylab("log2 expression level") +
  ggtitle("AD vs PSP top 25 significant genes. Number of samples: AD = 82, PSP = 84") +
  theme_light() +
  theme(plot.title = element_text(size = 16, face = "bold"))
ggsave("./figs_transcr/boxplot_ad_vs_psp.png", dpi = 300, scale = 2, width = 6, height = 3, units = "in") 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_transcr/boxplot_ad_vs_psp.png") 


## -----------------------------------------------------------------------------
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
} 


## -----------------------------------------------------------------------------
ccol <- c(1:31) # Control columns in 'data_all'
acol <- c(32:115) # AD columns in 'data_all'
pcol <- c(116:200) # PSP columns in 'data_all'

## We remove 3 samples that are absent in the transcriptomics.
acol <- acol[-c(29, 38)]
pcol <- pcol[-4] 


## -----------------------------------------------------------------------------
## AD proteins in AD samples
ad_pr_top <- data_all %>%
  `rownames<-`(.$PIDshort %>% make.unique) %>%
  filter(rownames(.) %in% rownames(AvsC_top)) %>%
  ## filter(rownames(.) %in% ad_vs_control_prot$PIDshort) %>%
  dplyr::select(all_of(acol))

## AD proteins in PSP samples
psp_pr_top <- data_all %>%
  `rownames<-`(.$PIDshort %>% make.unique) %>%
  filter(rownames(.) %in% rownames(AvsC_top)) %>%
  ## filter(rownames(.) %in% ad_vs_control_prot$PIDshort) %>%
  dplyr::select(all_of(pcol))

## AD proteins in Control samples
control_pr_top <- data_all %>%
  `rownames<-`(.$PIDshort %>% make.unique) %>%
  filter(rownames(.) %in% rownames(AvsC_top)) %>%
  ## filter(rownames(.) %in% ad_vs_control_prot$PIDshort) %>%
  dplyr::select(all_of(ccol)) 


## -----------------------------------------------------------------------------
## data_all contains log2 transformed data
ad_pr_top_in_all <- data_all %>%
  `rownames<-`(.$PIDshort %>% make.unique) %>%
  filter(rownames(.) %in% rownames(AvsC_top)) %>%
  ## filter(rownames(.) %in% ad_vs_control_prot$PIDshort) %>%
  dplyr::select(all_of(c(acol, pcol, ccol))) 


## -----------------------------------------------------------------------------
ad_pr_top_cor_all <- rcorr(t(as.matrix(ad_pr_top_in_all)), t(as.matrix(ad_pr_top_in_all)),
                     type = "spearman") 


## -----------------------------------------------------------------------------
png("./figs_cor/prot_ad_in_all.png", width = 2000, height = 2000, res = 300)
corrplot(ad_pr_top_cor_all$r[1:nrow(ad_pr_top), 1:nrow(ad_pr_top)],
         p.mat = ad_pr_top_cor_all$P[1:nrow(ad_pr_top), 1:nrow(ad_pr_top)],
         sig.level = 0.05, insig = "blank",
         ## type="upper", order="hclust",
         col = colorRampPalette(c("darkblue", "white", "green"))(200),
         tl.col = "black", tl.cex = 0.4, tl.srt = 45)
mytitle = "Number of samples: AD = 82, PSP = 84, Control = 31"
mtext(side=3, line=2, at=-0.07, adj=0, cex=1, mytitle)
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_cor/prot_ad_in_all.png") 


## -----------------------------------------------------------------------------
png("./figs_cor/prot_ad_in_all_numb.png", width = 2000, height = 2000, res = 300)
corrplot(ad_pr_top_cor_all$r[1:nrow(ad_pr_top), 1:nrow(ad_pr_top)],
         p.mat = ad_pr_top_cor_all$P[1:nrow(ad_pr_top), 1:nrow(ad_pr_top)],
         sig.level = 0.05, insig = "blank",
         ## type="upper", order="hclust",
         method = "number", number.cex = 0.2,
         col = colorRampPalette(c("darkblue", "white", "green"))(200),
         tl.col = "black", tl.cex = 0.4, tl.srt = 45)
mytitle = "Number of samples: AD = 82, PSP = 84, Control = 31"
mtext(side=3, line=2, at=-0.07, adj=0, cex=1, mytitle)
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_cor/prot_ad_in_all_numb.png") 


## -----------------------------------------------------------------------------
ad_top_r_all <- ad_pr_top_cor_all$r[1:nrow(ad_pr_top), 1:nrow(ad_pr_top)] %>%
  data.frame() %>%
  filter(rownames(.) %in% rownames(AvsC_top))
  ## filter(rownames(.) %in% ad_vs_control_prot$PIDshort)

write.csv(ad_top_r_all, "./results_cor/ad_prot_cor_coeffs_all.csv") 


## -----------------------------------------------------------------------------
write.csv(flattenCorrMatrix(cormat = ad_pr_top_cor_all$r[1:nrow(ad_pr_top),
                                                     1:nrow(ad_pr_top)],
                            pmat = ad_pr_top_cor_all$P[1:nrow(ad_pr_top),
                                                   1:nrow(ad_pr_top)]),
          "./results_cor/ad_prot_cor_coeffs+p-values_all.csv") 


## -----------------------------------------------------------------------------
transcripts_norm <- cpm(dgeObj, log = TRUE,
                        normalized.lib.sizes = TRUE) %>%
  data.frame() 


## -----------------------------------------------------------------------------
ad_all_tr_top <- transcripts_norm[c(ad_cols, psp_cols, control_cols)] %>%
  filter(rownames(.)  %in%  ad_vs_control_top_tr$TranscriptIDs) %>%
  `rownames<-`(rownames(all_samples_ad_vs_control_tr))

ad_tr_top_cor_all <- rcorr(t(as.matrix(ad_all_tr_top)),
                     t(as.matrix(ad_all_tr_top)),
                     type = "spearman") 


## -----------------------------------------------------------------------------
png("./figs_cor/transcr_ad_in_all.png", width = 2000, height = 2000, res = 300)
corrplot(ad_tr_top_cor_all$r[1:nrow(all_samples_ad_vs_control_tr),
                         1:nrow(all_samples_ad_vs_control_tr)],
         p.mat = ad_tr_top_cor_all$P[1:nrow(all_samples_ad_vs_control_tr),
                                 1:nrow(all_samples_ad_vs_control_tr)],
         ## type="upper", order="hclust",
         sig.level = 0.05, insig = "blank",
         col = colorRampPalette(c("darkblue", "white", "green"))(200),
         tl.col = "black", tl.cex = 0.4, tl.srt = 45)
mytitle = "Number of samples: AD = 82, PSP = 84, Control = 31"
mtext(side=3, line=2, at=-0.07, adj=0, cex=1, mytitle)
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_cor/transcr_ad_in_all.png") 


## -----------------------------------------------------------------------------
png("./figs_cor/transcr_ad_in_all_numb.png", width = 2000, height = 2000, res = 300)
corrplot(ad_tr_top_cor_all$r[1:nrow(all_samples_ad_vs_control_tr),
                         1:nrow(all_samples_ad_vs_control_tr)],
         p.mat = ad_tr_top_cor_all$P[1:nrow(all_samples_ad_vs_control_tr),
                                 1:nrow(all_samples_ad_vs_control_tr)],
           sig.level = 0.05, insig = "blank",
           ## type="upper", order="hclust",
           method = "number", number.cex = 0.2,
           col = colorRampPalette(c("darkblue", "white", "green"))(200),
           tl.col = "black", tl.cex = 0.4, tl.srt = 45)
mytitle = "Number of samples: AD = 82, PSP = 84, Control = 31"
mtext(side=3, line=2, at=-0.07, adj=0, cex=1, mytitle)
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_cor/transcr_ad_in_all_numb.png") 


## -----------------------------------------------------------------------------
write.csv(ad_tr_top_cor_all$r[1:nrow(all_samples_ad_vs_control_tr),
                          1:nrow(all_samples_ad_vs_control_tr)],
          "./results_cor/ad_transcripts_cor_r_all.csv")

## knitr::kable(ad_tr_top_cor$r[1:nrow(all_samples_ad_vs_control_tr), 1:nrow(all_samples_ad_vs_control_tr)], caption = "Correlation coefficients for top 25 AD genes") 


## -----------------------------------------------------------------------------

write.csv(flattenCorrMatrix(cormat = ad_tr_top_cor_all$r,
                            pmat = ad_tr_top_cor_all$P),
          "./results_cor/ad_transcripts_cor_coeffs+p-values_all.csv")

## knitr::kable(ad25P, caption = "Correlation coefficient's p-values for top 25 AD proteins") 


## -----------------------------------------------------------------------------
i <- nrow(ad_pr_top_in_all)
j <- nrow(ad_all_tr_top)

ad_pr_tr_cor_all <- rcorr(t(as.matrix(ad_all_tr_top)),
      t(as.matrix(ad_pr_top_in_all)),
      type = "spearman")


## write_rds(list(A=A, B=B), "./ab.rds")
## ab <- readRDS("./ab.rds")
## rcorr(t(T).t(P))$r[number_of_transcripts+1 : number_of_transcripts+number_of_proteins,
##                    1: number_of_transcripts] 


## -----------------------------------------------------------------------------
png("./figs_cor/trpr_ad_in_all.png", width = 2000, height = 3300, res = 300)
corrplot(ad_pr_tr_cor_all$r[(j+1) : (i+j), 1:j],
         p.mat = ad_pr_tr_cor_all$P[(j+1) : (i+j), 1:j],
         ## type="upper", order="hclust",
         sig.level = 0.05, insig = "blank",
         col = colorRampPalette(c("darkblue", "white", "green"))(200),
         tl.col = "black", tl.cex = 0.3, tl.srt = 45)
mytitle = "Number of samples: AD = 82, PSP = 84, Control = 31"
mtext(side=3, line=2, at=-0.07, adj=0, cex=1, mytitle)
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_cor/trpr_ad_in_all.png") 


## -----------------------------------------------------------------------------
png("./figs_cor/trpr_ad_in_all_numb.png", width = 2000, height = 3300, res = 300)
corrplot(ad_pr_tr_cor_all$r[(j+1) : (i+j), 1:j],
         p.mat = ad_pr_tr_cor_all$P[(j+1) : (i+j), 1:j],
         ## type="upper", order="hclust",
         sig.level = 0.05, insig = "blank",
         method = "number", number.cex = 0.2,
         col = colorRampPalette(c("darkblue", "white", "green"))(200),
         tl.col = "black", tl.cex = 0.3, tl.srt = 45)
mytitle = "Number of samples: AD = 82, PSP = 84, Control = 31"
mtext(side=3, line=2, at=-0.07, adj=0, cex=1, mytitle)
dev.off() 


## ----echo=FALSE, label="", fig.show='asis', fig.cap=""------------------------
include_graphics("./figs_cor/trpr_ad_in_all_numb.png") 


## -----------------------------------------------------------------------------
write.csv(ad_pr_tr_cor_all$r[(j+1) : (i+j), 1:j],
          "./results_cor/ad_prot_transcr_cor_r_all.csv") 


## -----------------------------------------------------------------------------
write.csv(flattenCorrMatrix(cormat = ad_pr_tr_cor_all$r[(j+1) : (i+j), 1:j],
                            pmat = ad_pr_tr_cor_all$P[(j+1) : (i+j), 1:j]),
          "./results_cor/ad_prot_transcr_cor_coeffs+p-values_all.csv") 

