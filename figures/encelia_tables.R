
t1 = read.tree("~/Dropbox/Encelia/analysis/phylogeny/concatenated/RAxML_bestTree.concatenated0.6")
dropped = c("ART")
tips = c(t1$tip.label, dropped)
d = read_csv("~/Dropbox/Encelia/analysis/phylogeny/morphological_data/Encelia_Samples - GENERAL.csv")
cg = read_csv("~/Dropbox/Encelia/analysis/phylogeny/morphological_data/Encelia_Samples - COMMON_GARDEN.csv")

dd = cbind(tips, d[match(tips, d$PLANT_ID), 
  c("PUTATIVE_SPECIES", "NOTES", "LOCALITY", "LATITUDE", "LONGITUDE")])

d = read.csv("~/Desktop/samples.csv", stringsAsFactors = F)
a = read.csv("~/Desktop/site_counts.csv", stringsAsFactors = F)
b = read.csv("~/Desktop/read_counts.csv", stringsAsFactors = F)
c = read.csv("~/Desktop/contig_counts.csv", stringsAsFactors = F)
d1 = left_join(d, a)
d2 = left_join(d1, b)
d3 = left_join(d2, c)
write.csv(d3, "~/Desktop/samples2.csv")
