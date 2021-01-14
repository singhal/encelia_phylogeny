library(ape)

trees = list.files("/Volumes/heloderma4/sonal/encelia/phylogeny_v2/trees_40_0.05_0.6_0.3_cov2/", pattern = "bestTree", full.names=T)
num_inds = 72

miss = c(0.6, 0.7, 0.8)
drop = c("ART", "BAI", "FAR-FAR-DVSC-2", "GER-1")

for (j in 1:length(miss)) {
	for (i in 1:length(trees)) {
		t = read.tree(trees[i])
		t = drop.tip(t, drop)
		tips = length(t$tip.label)
		# print(Nnode(t))
		t1 = di2multi(t, 2e-6)
		# print(Nnode(t1))
		if (tips / num_inds > miss[j]) {
			write.tree(t1, file=paste("astral.miss", miss[j],  ".trees", sep=""), append = T)
		}
	}
}
