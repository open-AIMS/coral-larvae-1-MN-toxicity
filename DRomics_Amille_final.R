# DRomics using RNASeq count data from 1-MN larval Acropora millepora study
# See: vignette("DRomics_vignette")

library(limma)
library(DRomics)
library(sva)
library(tidyverse)
library(readxl)
library(vegan)
library(ssdtools)

#---- import count data with 1st row indicating treatment concentration ----
# As per DESeq2 recommendations, genes with >10 reads were removed prior to formatting for DRomics

RNASeq.dat <- read.table(file = "allcounts4DRomics.txt", header = FALSE)

o.RNAseq <- RNAseqdata(RNASeq.dat, transfo.method = "rlog")

plot(o.RNAseq, cex.main = 0.8, col = "green", range = 1e6)

#---- correct for RNA extraction batch ----

metaD <- read.csv("Metadata.csv", row.names = 1, header = TRUE)
o.RNAseq$batch <- as.factor(metaD$Extraction_Batch)
pca <- PCAdataplot(o.RNAseq, batch = o.RNAseq$batch)
pca
BECcounts <- ComBat_seq(as.matrix(o.RNAseq$raw.counts), 
                        batch = as.factor(o.RNAseq$batch), 
                        group = as.factor(o.RNAseq$dose)) 
BECdata4DRomics <- formatdata4DRomics(signalmatrix = BECcounts, 
                                      dose = o.RNAseq$dose)
#write.csv(BECdata4DRomics, file = "BECdata.csv")
o.BEC <- RNAseqdata(BECdata4DRomics)
p <- PCAdataplot(o.BEC, batch = o.RNAseq$batch) #name so pc's can be extracted/plotted
p
plot(p$data$dose, p$data$PC1) #dose is a factor so use metaD to get scatterplot vs boxplot
plot(p$data$dose, p$data$PC2)
plot(metaD$WAF, p$data$PC1)
plot(metaD$WAF, p$data$PC2)

#---- Make nicer PCA plots ----

ggplot(pca$data, aes(PC1, PC2)) + 
  geom_point(aes(shape=pca$data$batch, color=pca$data$dose), size = 4) +
  scale_color_manual(values = c("#b3d4ff", "#3d3dff", "#00008b", "#64ff64" , "#008b00", "#005000", "#ff8cff", "#b300b3", "#4d004d")) +
  labs(shape = "RNA\nBatch", color = expression("1-MN"~ (mu*g~L^-1)), x = "PC1 (21.93%)", y = "PC2 (4.69%)" ) +
  theme_classic()
#ggsave("PCA_noBatchCorrection.pdf", width = 8, height = 5, dpi = 500)
#ggsave("PCA_noBatchCorrection.jpeg", width = 8, height = 5, dpi = 500)

ggplot(p$data, aes(PC1, PC2)) + 
  geom_point(aes(shape=p$data$batch, color=p$data$dose), size = 4) +
  scale_color_manual(values = c("#b3d4ff", "#3d3dff", "#00008b", "#64ff64" , "#008b00", "#005000", "#ff8cff", "#b300b3", "#4d004d")) +
  labs(shape = "RNA\nBatch", color = expression("1-MN"~ (mu*g~L^-1)), x = "PC1 (9.11%)", y = "PC2 (4.06%)" ) +
  theme_classic()
#ggsave("PCA_BatchCorrected.pdf", width = 8, height = 5, dpi = 500)
#ggsave("PCA_BatchCorrected.jpeg", width = 8, height = 5, dpi = 500)

#---- select significantly responding items and Fit dose-response models (uses drc) ----

s_quad <- itemselect(o.BEC, select.method = "quadratic", FDR = 0.05) #0.01 

f <- drcfit(s_quad, progressbar = FALSE)

head(f$fitres, 20)

head(f$unfitres) #look at table of unsuccessful fits
table(f$unfitres$cause)

#---- Plot all fit genes and save drc object ----

plotfit2pdf(f, dose_log_transfo = TRUE, path2figs = "C:/Users/hluter/OneDrive - Australian Institute of Marine Science/Environmental Impacts/Current Projects/Molecular and Microbial Ecotox/1-MN coral larval exp/SeqProc_R")

save(f, file = "FittedGenes_FDR0.05.RData")

#---- BMD calculation ----

r <- bmdcalc(f, z = 1, x = 10)

head(r$res, 10)

#write.csv(r$res, file = "fittedGenes_wBMD.csv") 

plot(r, BMDtype = "zSD", plottype = "ecdf")

plot(r, BMDtype = "zSD", plottype = "ecdf", dose_log_transfo = TRUE)

plot(r, BMDtype = "zSD", plottype = "density", by = "trend")

bmdplotwithgradient(r$res, BMDtype = "zSD",
                    facetby = "trend",
                    shapeby = "model",
                    line.size = 1.2) +
  scale_x_log10(limits=c(0.2, 1000))

#---- Calculate 95% confidence intervals using bootstraps min niter value of 1000 ----

b <- bmdboot(r, niter = 1000, progressbar = FALSE)

plot(b, BMDtype = "zSD", by = "trend", dose_log_transfo = TRUE)

plot(f, BMDoutput=b, dose_log_transfo = TRUE) +
  scale_x_log10(limits=c(0.2, 1000))

# Save DRomics workflow results with BMDs and CIs

plotfit2pdf(f, BMDoutput=b, dose_log_transfo = TRUE, path2figs = "C:/Users/hluter/OneDrive - Australian Institute of Marine Science/Environmental Impacts/Current Projects/Molecular and Microbial Ecotox/1-MN coral larval exp/SeqProc_R") +
  scale_x_log10(limits=c(0.2, 1000))

save(b, file = "FullDRomicsResults_FDR0.05.RData")
