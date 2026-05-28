#### Repeat settlement stats for 1-MN experiment in updated bayesnec workflow

library(tidyverse)
library(bayesnec)
library(purrr)
library(furrr)
source("functions.R")


dat <- read_csv("2021.12_Amille_1-MN_meta.csv") |> 
  data.frame() |> 
  dplyr::mutate(raw.x=as.numeric(as.character(raw.x)),
                raw.x=ifelse(raw.x==0, 0.1, raw.x), # deal with 0 concentration values
                nom.x=as.factor(as.character(nom.x)),
                suc=as.integer(suc), # ensure integers
                log.x=log10(raw.x), # create log values (doesn't really matter if it's log or log10() unless there's a difference in the scaling of the x-data)
                tot=as.integer(tot),
                prop=as.numeric(suc/tot))

datF <- dat |> 
  dplyr::filter(use_in_model=="Yes") # remove samples with only dead material

plot(datF$log.x, datF$prop)

mod <- bnec(suc|trials(tot) ~ crf(log.x, model = "nec3param"),
            data = datF,   family = beta_binomial(link="identity"),
          #  seed = 333, control = list(adapt_delta = 0.99),
            iter = 20000, warmup = 19000)

# gives an error in expand_manec(mod_fits, formula = formulas, x_range = x_range) no models fit? - try defaults

mod <- bnec(suc|trials(tot) ~ crf(log.x, model = "decline"),
            data = datF)

check_dvts(mod) # need to remove ecxexp and ecxll5 which had high divergent transitions

mod.clean <- amend(mod, drop = c("ecxexp", "ecxll5"))

rhat(mod.clean, rhat_cutoff = 1.01) # all good

check_chains(mod.clean)

summary(mod.clean) #binomial, 8000 posterior draws per mod

save(mod.clean, file = "1-MN_cleanFits.RData")

EC10 <- ecx(mod.clean, ecx_val = 10, xform = function(x){10^x})
EC50 <- ecx(mod.clean, ecx_val = 50, xform = function(x){10^x})

autoplot(mod.clean, xform = function(x){10^x}, nec = FALSE) +
  scale_x_continuous(trans = "log10", labels = scales::label_number(drop0trailing = TRUE)) +
  scale_y_continuous(labels = function(x) x * 100) +
  labs (y = "Larval metamorphosis (%)", x = "1-methylnaphthalene (µg/L)") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) 



