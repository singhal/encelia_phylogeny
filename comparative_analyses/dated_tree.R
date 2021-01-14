library(ape)
library(phytools)
library(phangorn)

t = read.tree("~/Dropbox/Encelia/analysis/phylogeny/external_dating/v0.1/ALLMB.tre")
tips1 = t$tip.label[grep("Encelia_", t$tip.label)]
tips2 = t$tip.label[grep("Enceliopsis_", t$tip.label)]
tips3 = t$tip.label[grep("Geraea_", t$tip.label)]
# tips4 = t$tip.label[grep("Xylorhiza", t$tip.label)]
mrca = findMRCA(t, c(tips1, tips2, tips3))
# climb up
anc1 = Ancestors(t, mrca, type = "parent")
anc2 = Ancestors(t, anc1, type = "parent")
desc = Descendants(t, anc2, type = "tips")[[1]]
t1 = drop.tip(t, setdiff(t$tip.label, t$tip.label[desc]))
write.tree(t1, "~/Desktop/encelia.tre")
