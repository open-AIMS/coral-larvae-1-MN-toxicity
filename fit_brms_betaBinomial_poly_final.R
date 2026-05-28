
library(rlang)
require(brms)
library(tidyverse)

####---- select DRomics fitted genes from batch corrected input data ----####

BECdat <- read.csv('BatchCorrectedData.csv')
fitGenes <- read.csv('DromicsFittedGenesFDR0.05.csv') |> 
  dplyr::inner_join(BECdat, fitGenes, by = "GeneID") |> 
  dplyr::mutate(max = rowSums(across(where(is.numeric))))

dat_abundant <- fitGenes |> 
  tidyr::pivot_longer(cols = starts_with("X"), names_to = "SampleID", values_to = "y")
dat <- read.csv('Metadata_1-MN.csv') |> 
  dplyr::left_join(dat_abundant) |> 
  dplyr::mutate(x = as.numeric(as.character(Measured)), 
                x = ifelse(x==0, 0.3, x),     
                log_x = log(x))
head(dat)
head(dat_abundant)

summary_dat <- data.frame(gene =unique(paste("gene", dat$GeneID, sep = "_")), 
                          n = unique(paste("n", dat$GeneID, sep = "_")))

dat_max <- dat|> 
  dplyr::select(-y) |> 
  tidyr::pivot_wider(names_from = GeneID, values_from = max, names_prefix = "n_")
head(dat_max)

dat_use <- dat|> 
  dplyr::select(-max) |> 
  tidyr::pivot_wider(names_from = GeneID, values_from = y, names_prefix = "gene_") |> 
  dplyr::left_join(dat_max)
colnames(dat_use)

bf_3k <- paste0("bf(", summary_dat$gene, " | trials(", summary_dat$n, ") ~ poly(log_x, 3, raw=TRUE))",
                collapse = " + ") |>
  parse(text = _) |>
  eval()

save(dat_use, bf_3k, file = "DRomicsGenes_brms_dat_3k_poly.RData") 

#### create data for 4th and 5th order polynomials to see which is best

# to run on on HPC ------------------------------------------------------------
load(file = "DRomicsGenes_brms_dat_3k_poly.RData")

#load(file = "DRomicsGenes_brms_dat_4k_poly.RData")
#load(file = "DRomicsGenes_brms_dat_5k_poly.RData")

iter <- 2000
fit3k <- brm(bf_3k + set_rescor(FALSE), data = dat_use, chains = 4, cores = 4, 
           #threads = threading(4),
           init = 0, seed = 100,
           iter = iter, warmup = iter-1000, family = "beta_binomial")

#save(fit3k, file = "fit3k_bb_DRomicsGenes_poly_iter1000.RData")
# end hpc ----------------------------------------------------------------------

load(file = "fit3k_bb_DRomicsGenes_poly_iter1000.RData")  

#load(file = "condEffects_bb_5k_TITAN.indi_TSV") #HPC
cond_effects5k <- conditional_effects(fit5k, resolution = 50)
  
plot_dat5k <- lapply(cond_effects5k, FUN = function(x){
    x |> dplyr::select(log_x, estimate__, lower__, upper__)
  })

plot_dat5k <- plot_dat5k |> 
  dplyr::bind_rows(.id = "GeneID") |> 
  tidyr::separate_wider_delim(GeneID, ".", names = c("GeneID", "junk")) |> 
  dplyr::select(GeneID, log_x, estimate__, lower__, upper__)  |>
  dplyr::mutate(GeneID=stringr::str_sub(GeneID, 5))

#save(plot_dat5k, dat_use, file = "plot_dat_5k_bb_DRomicsGenes.RData") 

load(file = "fit4k_bb_DRomicsGenes_poly_iter1000.RData")  

cond_effects4k <- conditional_effects(fit4k, resolution = 50)

plot_dat4k <- lapply(cond_effects4k, FUN = function(x){
  x |> dplyr::select(log_x, estimate__, lower__, upper__)
})

plot_dat4k <- plot_dat4k |> 
  dplyr::bind_rows(.id = "GeneID") |> 
  tidyr::separate_wider_delim(GeneID, ".", names = c("GeneID", "junk")) |> 
  dplyr::select(GeneID, log_x, estimate__, lower__, upper__)  |>
  dplyr::mutate(GeneID=stringr::str_sub(GeneID, 5))

save(plot_dat4k, dat_use, file = "plot_dat_4k_bb_DRomicsGenes.RData") 

load(file = "fit3k_bb_DRomicsGenes_poly_iter1000.RData")  

cond_effects3k <- conditional_effects(fit3k, resolution = 50)

plot_dat3k <- lapply(cond_effects3k, FUN = function(x){
  x |> dplyr::select(log_x, estimate__, lower__, upper__)
})

plot_dat3k <- plot_dat3k |> 
  dplyr::bind_rows(.id = "GeneID") |> 
  tidyr::separate_wider_delim(GeneID, ".", names = c("GeneID", "junk")) |> 
  dplyr::select(GeneID, log_x, estimate__, lower__, upper__)  |>
  dplyr::mutate(GeneID=stringr::str_sub(GeneID, 5))

save(plot_dat3k, dat_use, file = "plot_dat_3k_bb_DRomicsGenes.RData") 

################ NSECs ##############

library(toxval)

toxval_nsec_lower <- nsec_multi(object=fit3k,
                            
                                resolution = 20,
                                
                                type = "lower", #both
                                
                                posterior = FALSE,
                                
                                xform = exp,
                                
                                x_var="log_x",  trials_var = "n_", multi_var = "gene_")
save(toxval_nsec_lower, file = "toxval_nsec_lower_bb_3k_poly_DRomicsGenes_backtransformed.RData")

################# Plots - log scale ##############

library(stringr)

#### load fit/nsec data
#load("plot_dat_3k_bb_DRomicsGenes.RData")
#load("toxval_nsec_lower_bb_3k_poly_DRomicsGenes_backtransformed.RData")

# remove 'GeneID' prefix from plot_dat5k and toxval_nsec dataframes to match dat for plotting

####### lowest nsec 
plot_dat3k$GeneID <- str_replace(plot_dat3k$GeneID, "GeneID", "")

plot_dat3k_52 <- plot_dat3k |> 
      filter(GeneID %in% toxval_nsec_lw_F2$var)

dat52 <- dat |> 
      filter(GeneID %in% toxval_nsec_lw_F2$var)

toxval_nsec_lw_F2$var <- str_replace(toxval_nsec_lw_F2$var, "gene_", "") #toxval_nsec_lower for all
toxval_nsec_lw_F2 <- toxval_nsec_lw_F2 |> 
  mutate_at(c(1:4), as.numeric) |> #for some reason this data frame coming up as characters - 'both' was numeric to begin with
  mutate(log_NSEC = log(val),
         log_lw = log(lw),
         log_up = log(up),
         GeneID = var)

plot_dat3k_52 <- plot_dat3k_52 |>
    mutate(x = exp(log_x))

plot.brm.fits.lower <-   ggplot(plot_dat3k_52, 
                                aes(x = log_x, 
                                    y = estimate__)) +
  geom_ribbon(aes(x=log_x, 
                  ymin=lower__, ymax = upper__), alpha=0.3) + 
  geom_line(aes(x=log_x, y=estimate__), color = 'grey30') + 
  #scale_x_continuous(trans = "log", labels = scales::number_format(accuracy = 0.01)) +
  facet_wrap(~ GeneID, scales="free_y")  + 
  geom_point(data = dat52, 
             aes(x = log_x, y = y), size = 0.5) +
  labs(x = expression("log 1-MN (ug/L)"), y = "gene expression") +
  
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(margin = unit(c(6, 6, 0, 0), "mm")),
    axis.title.x = element_text(margin = unit(c(6, 6, 0, 0), "mm")),
    legend.key.size = unit(0.5, 'cm'),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)) 

plot.brm.fits.lower

# Add NSECs
plot.brm.fits.lowerNSEC <- plot.brm.fits.lower + 
  geom_vline(data = toxval_nsec_lw_F2,
             aes(xintercept = log_NSEC, color = "NSEC", linetype = "dotdash")) +
  geom_vline(data = toxval_nsec_lw_F2,
             aes(xintercept = log_lw, color = "LW", linetype = "")) +
  geom_vline(data = toxval_nsec_lw_F2,
             aes(xintercept = log_up, color = "UP", linetype = "")) +
  facet_wrap(~ GeneID, scales = "free_y") +
 scale_color_manual(
   name = "",
  values = c("NSEC" = "#4DB6D0", "LW" = "#D9717D", "UP" = "#D9717D"))

plot.brm.fits.lowerNSEC

#ggsave("plot.brm.fits.bb.lowerNSEC_poly3_LogScale.jpg", plot.brm.fits.lowerNSEC, width = 15, height = 10, limitsize = FALSE)

################# Plots - backtransformed ##############

plot.brm.fits.lower.bt <-   ggplot(plot_dat3k_52, 
                                   aes(x = x, 
                                       y = estimate__)) +
  geom_ribbon(aes(x=x, 
                  ymin=lower__, ymax = upper__), alpha=0.3) + 
  geom_line(aes(x=x, y=estimate__), color = 'grey') + 
  #scale_x_continuous(trans = "log", labels = scales::number_format(accuracy = 1)) +
  scale_x_continuous(trans = "log", breaks = c(1, 36, 882)) +
  facet_wrap(~ GeneID, scales="free_y")  + 
  geom_point(data = dat52, 
             aes(x = x, y = y), size = 0.5) +
  labs(x = expression("1-methylnaphthalene"~(μg~"L"^{-1})), y = "gene expression") +
  
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(margin = unit(c(6, 6, 0, 0), "mm")),
    axis.title.x = element_text(margin = unit(c(6, 6, 0, 0), "mm")),
    legend.key.size = unit(0.5, 'cm'),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)) 

plot.brm.fits.lower.bt

# Add NSECs
plot.brm.fits.lowerNSEC.bt <- plot.brm.fits.lower.bt + 
  geom_vline(data = toxval_nsec_lw_F2,
             aes(xintercept = val, color = "NSEC", linetype = "dotdash")) +
  geom_vline(data = toxval_nsec_lw_F2,
             aes(xintercept = lw, color = "LW", linetype = "")) +
  geom_vline(data = toxval_nsec_lw_F2,
             aes(xintercept = up, color = "UP", linetype = "")) +
  facet_wrap(~ GeneID, scales = "free_y") +
  scale_color_manual(
    name = "",
    values = c("NSEC" = "#4DB6D0", "LW" = "#D9717D", "UP" = "#D9717D"))

plot.brm.fits.lowerNSEC.bt

ggsave("plot.brm.fits.bb.lowerNSEC_poly3_backtransformed_52genes.pdf", plot.brm.fits.lowerNSEC.bt, width = 15, height = 10, limitsize = FALSE)

####### Compare models with varying k values ######

# https://paul-buerkner.github.io/brms/articles/brms_multivariate.html

fit5k_test <- add_criterion(fit5k, "loo")
summary(fit5k_test)

fit4k_test <- add_criterion(fit4k, "loo")
fit3k_test <- add_criterion(fit3k, "loo")


loo(fit5k_test, fit4k_test, fit3k_test)

loo_model_weights(fit5k_test, fit4k_test, fit3k_test)

loo_compare(fit5k_test, fit4k_test, fit3k_test) # condensed version of loo command above

####### Selecting best model info ######

# The model with the lower LOOIC is “better” (LOOIC performs like AIC)
# Rather than looking at LOOIC, you should look at the elpd (expected log predictive density), as LOOIC is just -2 * elpd
# For elpd, the larger the better, but there is no intrinsic meaning in that number on its own
# You want the elpd_diff to be > 4-5 times bigger than the se_diff

#------ Above agrees with weights, 3k best (also most simple) for this dataset


####### Filter toxval_nsec_lower to only exclude genes with NSEC at at highest concentration ######

#load(file = "toxval_nsec_lower_bb_5k_DRomicsGenes_backtransformed.RData")

toxval_nsec_lw_F2 <- filter(toxval_nsec_lower, val != 882) #52 genes

length(which(toxval_nsec_lw_F2$direction == "dec")) #27 decreasing genes
length(which(toxval_nsec_lw_F2$direction == "inc")) #25 increasing genes

#write.csv(toxval_nsec_lw_F2, file = "lower_NSECs_filtered_52genes.csv")
#save(toxval_nsec_lw_F2, file = "Lower_NSECs_filteredData_52genes.RData")

#### SSD #### 

library(ssdtools)

#load(file = "Lower_NSECs_filteredData_52genes.RData")

ssd_dat <- toxval_nsec_lw_F2
names(ssd_dat)[names(ssd_dat) == 'val'] <- 'Conc'
ssd_dat <- ssd_dat |> 
  mutate_at(vars(lw, Conc, up, ref), as.numeric)

## Fit models

fit <- ssd_fit_bcanz(ssd_dat, left = 'Conc') #updated to default models agreed upon in new release of ssdtools

ssd_plot_cdf(fit, delta = 20) 

## Calculate goodness of fit 

ssd_gof(fit)

## Calculate HC5 values

set.seed(99)
hc1 <- ssd_hc(fit, proportion = 0.01, ci=TRUE, nboot = 1000L, min_pboot = 0.8)
hc5 <- ssd_hc(fit, proportion = 0.05, ci=TRUE, nboot = 1000L, min_pboot = 0.8)
hc10 <- ssd_hc(fit, proportion = 0.10, ci=TRUE, nboot = 1000L, min_pboot = 0.8)
hc50 <- ssd_hc(fit, proportion = 0.50, ci=TRUE, nboot = 1000L, min_pboot = 0.8)
print(hc1) #52 genes = 1.13 (0.29-4.16)
print(hc5) # 4.75 (1.31-15.3); 52 genes = 6.44 (2.65-14.4)
print(hc10) # 8.75 (2.9-24.4); 52 genes = 13.8 (6.64-27.9)
print(hc50) # 41.5 (25.2-70.8); 52 genes = 101 (73.5-140)

## Calculate model-averaged predictions

set.seed(99)
pred <- predict(fit, ci = TRUE, nboot = 1000, min_pboot = 0.8) 

#NOTE: if min_pboot left as default, no upper or lower confidence limits are given!

## save fit, hc5, pred 

save(fit, hc1, hc5, hc10, hc50, pred, file = "SSDfitData_52genes.RData")

## Plot

ssd_plot(ssd_dat, pred,
         #shape = "direction",
         color = "direction",
         hc = NULL,
         xlab = "1-methylnaphthalene (μg/L)",
         ylab = "Genes affected", ribbon = TRUE) +
         expand_limits(x = 3000) +
  #scale_shape_manual(values = c(25,17)) +
  scale_color_manual(values = c("#156082", "#ED9C3F")) +
  #geom_hcintersect(xintercept = c(4.75, 8.75, 41.5), yintercept = c(0.05, 0.1, 0.5)) +
  geom_hcintersect(xintercept = c(1.13, 6.44, 13.8, 101), yintercept = c(0.01, 0.05, 0.1, 0.5)) +
  theme_classic()
ggsave('SSD_52genes_HC1_5_10_50_genesByColor.pdf', width = 15, height = 15, units = 'cm') 
