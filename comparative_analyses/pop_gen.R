library(ggplot2)
library(cowplot)
library(readr)
library(tidyr)
library(dplyr)
library(patchwork)
theme_set(theme_cowplot())

setwd("~/Dropbox/Encelia/manuscripts/Encelia_Phylogeny/figures/")

#############################################
# functions
#############################################

read_div_files <- function(d) {
  dd = read.csv(d, stringsAsFactors = F)
  dd1 = dd[dd$denom > 0, ]
  pi1 = sum(dd1[dd1$metric == "pi1", "value"]) / sum(dd1[dd1$metric == "pi1", "denom"])
  pi2 = sum(dd1[dd1$metric == "pi2", "value"]) / sum(dd1[dd1$metric == "pi2", "denom"])
  dxy = sum(dd1[dd1$metric == "diff", "value"]) / sum(dd1[dd1$metric == "diff", "denom"])
  da = dxy - (pi1 + pi2) * 0.5
  denom = sum(dd1[dd1$metric == "diff", "denom"])
  dfst = dd1[dd1$metric == 'Fst', ]
  fstval = sum(dfst$denom * dfst$value) / sum(dfst$denom)
  res = c(unique(dd1$sp1), unique(dd1$sp2), dxy, da, pi1, pi2, denom, fstval, sum(dfst$denom))
  return(res)
}

################################################
# divergence
################################################

divfiles = list.files("~/Dropbox/Encelia/analysis/pop_gen/div/", pattern = "divergence", full.names = T)
totres = lapply(divfiles, read_div_files)
totres1 = as.data.frame(do.call(rbind, totres), stringsAsFactors = F)
names(totres1) = c("sp1", "sp2", "dxy", "da", 
                   "pi1", "pi2", "dxy_denom", "fst", "fst_denom")
totres1 = totres1[totres1$sp1 != totres1$sp2, ]

outs = c("Geraea_canescens", "Enceliopsis_covillei", "Xylorhiza_tortifolia")
totres2 = totres1[!(totres1$sp1 %in% outs), ]
totres2 = totres2[!(totres2$sp2 %in% outs), ]
fst = ggplot(totres2, aes(as.numeric(fst))) + 
  geom_histogram(binwidth = 0.05) + xlab(expression(F[ST]))
dxy = ggplot(totres2, aes(as.numeric(dxy))) + 
  geom_histogram(binwidth = 0.0008) + xlab(expression(d[xy]))
da = ggplot(totres2, aes(as.numeric(da))) + 
  geom_histogram(binwidth = 0.001) + xlab(expression(d[a]))
fig = (fst | dxy | da) + plot_annotation(tag_levels = 'A')
save_plot("pop_gen.pdf", fig, ncol = 3, base_width = 3, base_height = 2.5)

################################################
# diversity
################################################

pifiles = list.files("~/Dropbox/Encelia/analysis/pop_gen/pi/", pattern = "pi", full.names = T)
totpi = lapply(pifiles, read_csv)
totpi2 = do.call(rbind, totpi)
totpi3 = totpi2 %>% dplyr::filter(metric == "pi")
hist(totpi3$value)

# not worth pursuing for this paper