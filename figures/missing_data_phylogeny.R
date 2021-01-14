
t6 = read.tree("~/Dropbox/Encelia/analysis/phylogeny/astral/astral.miss0.6.tre")
t6 = root(t6, c("Enceliopsis_covillei", "Xylorhiza_tortifolia"),
         resolve.root = T)
t6 = drop.tip(t6, c("Enceliopsis_covillei", "Xylorhiza_tortifolia"))
t6$tip.label = gsub("Encelia_", "E. ", t6$tip.label)

t7 = read.tree("~/Dropbox/Encelia/analysis/phylogeny/astral/astral.miss0.7.tre")
t7 = root(t7, c("Enceliopsis_covillei", "Xylorhiza_tortifolia"),
          resolve.root = T)
t7 = drop.tip(t7, c("Enceliopsis_covillei", "Xylorhiza_tortifolia"))
t7$tip.label = gsub("Encelia_", "E. ", t7$tip.label)

t8 = read.tree("~/Dropbox/Encelia/analysis/phylogeny/astral/astral.miss0.8.tre")
t8 = root(t8, c("Enceliopsis_covillei", "Xylorhiza_tortifolia"),
          resolve.root = T)
t8 = drop.tip(t8, c("Enceliopsis_covillei", "Xylorhiza_tortifolia"))
t8$tip.label = gsub("Encelia_", "E. ", t8$tip.label)

png("~/Dropbox/Encelia/manuscripts/Encelia_Phylogeny/figures/missing_data_phylogeny.png",
    width = 6, height = 3, units = "in", res = 300)
par(mfrow= c(1,3), mar = c(0, 0, 3, 0), xpd = T)

plot(t6)
nodelabs = as.numeric(t6$node.label)
nodelabs[!complete.cases(nodelabs)] = 1.0
for (i in 1:length(nodelabs)) {
  if (nodelabs[i] > 0.95) {
    nodelabels("", i + Ntip(t6), frame = "none", pch = 16, cex = 0.7)
  } else {
    nodelabels("", i + Ntip(t6), frame = "none", pch = 16, 
               cex = 2.5, col = "lightblue") 
    nodelabels(nodelabs[i], i + Ntip(t6), frame = "none", cex = 0.5) 
  }
}
title(">60% complete")

plot(t7)
nodelabs = as.numeric(t7$node.label)
nodelabs[!complete.cases(nodelabs)] = 1.0
for (i in 1:length(nodelabs)) {
  if (nodelabs[i] > 0.95) {
    nodelabels("", i + Ntip(t7), frame = "none", pch = 16, cex = 0.7)
  } else {
    nodelabels("", i + Ntip(t7), frame = "none", pch = 16, 
               cex = 2.5, col = "lightblue") 
    nodelabels(nodelabs[i], i + Ntip(t7), frame = "none", cex = 0.5) 
  }
}
title(">70% complete")

plot(t8)
nodelabs = as.numeric(t8$node.label)
nodelabs[!complete.cases(nodelabs)] = 1.0
for (i in 1:length(nodelabs)) {
  if (nodelabs[i] > 0.95) {
    nodelabels("", i + Ntip(t8), frame = "none", pch = 16, cex = 0.7)
  } else {
    nodelabels("", i + Ntip(t8), frame = "none", pch = 16, 
               cex = 2.5, col = "lightblue") 
    nodelabels(nodelabs[i], i + Ntip(t8), frame = "none", cex = 0.5) 
  }
}

title(">80% complete")
dev.off()