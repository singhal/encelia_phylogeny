library(ape)
library(dplyr)
library(tidyr)

t = read.tree("~/Dropbox/Encelia/analysis/phylogeny/concatenated/RAxML_bestTree.concatenated0.6")
d = read.csv("~/Dropbox/Encelia/analysis/phylogeny/morphological_data/Encelia_Samples - GENERAL.csv", stringsAsFactors = F)
phy = table(d[ match(t$tip.label, d$PLANT_ID), "PUTATIVE_SPECIES" ])

m = read.csv("~/Dropbox/Encelia/analysis/phylogeny/morphological_data/morphological_data.csv", stringsAsFactors = F)
m2 = m %>% group_by(PUTATIVE_SPECIES, measurement) %>% summarize(count = n()) %>% ungroup()

keep = c("leaf_Area", "leaf_Round", 
         "BTSA",
         "Kshoot_mmol/s/m2/MPa", "leaf_gray_mean",
         "top_trichomes_per_mm2", "vessel_diam_um", "wood_density",
         "leaf_LMA")
m3 = m2 %>% filter(measurement %in% keep, 
                   ! PUTATIVE_SPECIES %in% c("Artemesia_tridentata", "Baiylea_multiradiata"))
m4 = m3 %>% spread(measurement, count)
m4$genetics = phy[m4$PUTATIVE_SPECIES]
write.csv(m4, "~/Dropbox/Encelia/manuscripts/Encelia_Phylogeny/figures/sample_counts.csv",
          row.names = F, quote = F)
