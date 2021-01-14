library(ape)

d = read.csv("~/Dropbox/Encelia/ddRAD/analysis/encelia_samples_v4.csv", stringsAsFactors = F)

rename <- function(t) {
  tips = gsub("_", "-", t$tip.label)
  t$tip.label = paste(d[match(tips, d$sample), "lineage"], t$tip.label)
  return(t)
}

reroot <- function(t) {
  t = unroot(t)
  outs = c("ENC_1", "ENC_2")
  outs = outs[outs %in% t$tip.label]
  t = root(t, outs, resolve.root = T)
}

tfiles = list.files("Desktop/untitled folder/", full.names = T)
tfiles = tfiles[grep("species", tfiles, invert= T)]
trees = lapply(tfiles, read.nexus)
trees = lapply(trees, reroot)
trees = lapply(trees, rename)


pdf("~/Desktop/trees.pdf", height = 10, width = 6)
for (i in 1:length(trees)) {
  par(mar = c(0, 0, 3, 0))
  plot(trees[[i]], cex = 0.7)
  title(tfiles[[i]])
}
dev.off()