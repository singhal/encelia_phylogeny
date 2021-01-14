library(ape)

d = read.csv("/Volumes/heloderma4/sonal/encelia/encelia_samples_v4.csv", 
	stringsAsFactors = F, na.string = c("", "NA"))
d = d[complete.cases(d$lineage), ]
indir = '/Volumes/heloderma4/sonal/encelia/phylogeny/trees_40_0.05_0.6_0.3_cov2/'
trees = list.files(indir, pattern = "bestTree", full.names = T)

parse_tree <- function(t) {
	t = read.tree(t)
	t = unroot(t)

	t = drop.tip(t, t$tip.label[!(t$tip.label %in% d$sample)])
	t1 = di2multi(t, 2e-6)
	lins = table(d[match(t1$tip.label, d$sample), 'lineage'])
	drop = vector('list', length(lins))
	for (i in 1:length(lins)) {
		if (lins[i] > 1) {
			sp = names(lins)[i]
			tips = d[d$lineage == sp, 'sample']
			lintips = t1$tip.label[t1$tip.label %in% tips]
			drop[[i]] = sample(lintips, size= (lins[i] - 1), replace = F)
			
		}
	}
	drop = unlist(drop)
	t2 = drop.tip(t1, drop)
	t2$tip.label = d[match(t2$tip.label, d$sample), 'lineage']

	outs = c("Xylorhiza_tortifolia", "Geraea_canescens")
	t2 = drop.tip(t2, outs)
	return(t2)
}


trees2 = lapply(trees, parse_tree)
nnode = unlist(lapply(trees2, function(x) {return(x$Nnode)}))
trees3 = trees2[which(nnode > 7)]
ntips = unlist(lapply(trees3, function(x) {return(length(x$tip.label))}))
trees4 = trees3[which(ntips > 13)]
class(trees4) <- "multiPhylo"
write.tree(trees4, "snaq.t14_n8.trees")

