library(ape)
library(phangorn)
library(readr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(gridGraphics)
library(geiger)
theme_set(theme_cowplot())

setwd("~/Dropbox/Encelia/manuscripts/Encelia_Phylogeny/figures/")

inds = read_csv("~/Dropbox/Encelia/analysis/hybrid_zones/Encelia_Samples - GENERAL.csv")
x = read_csv("~/Dropbox/Encelia/ddRAD/analysis/encelia_samples_v4.csv")

###################################
# species tree
###################################

# prep tree
t = read.tree("~/Dropbox/Encelia/analysis/phylogeny/astral/astral.miss0.6.brlen.dated.tre")
t = root(t, c("Enceliopsis_covillei", "Xylorhiza_tortifolia"),
         resolve.root = T)
t1 = read.tree("~/Dropbox/Encelia/analysis/phylogeny/astral/astral.miss0.6.tre")
t1 = minRotate(t1, setNames(1:Ntip(t), t$tip.label))
t$node.label = as.numeric(t1$node.label)
td = root(t, c("Enceliopsis_covillei", "Xylorhiza_tortifolia"),
         resolve.root = T)
td = read.tree(text = write.tree(ladderize(td, right = F)))

pdf("Encelia_species_phylogeny.astral_miss0.6.pdf", height = 6,
    width = 7)
par(mar = c(2, 0.5, 0.5, 12), xpd = T)
plot(td, show.tip.label = F)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
xscale = range(lastPP$xx)
tscale <- c(xscale[2] - xscale[1], 0)
beta <- diff(xscale)/diff(tscale)
alpha <- xscale[1] - beta * tscale[1]
labs = pretty(tscale)
labs = labs[ labs >= min(xscale) & labs <= max(xscale) ]
x <- beta * labs + alpha
axis(x, side = 1, labels = NA, tck = -0.01)
axis(side = 1, lwd = 0, line = -.8, at = x,
     labels = labs, cex.axis = 0.7)

tipnames = gsub("_", " ", td$tip.label)
tiplabels(tipnames, frame = "none", adj = 0, font = 3)
nodelabs = as.numeric(td$node.label)
nodelabs[!complete.cases(nodelabs)] = 1.0
for (i in 1:length(nodelabs)) {
  if (nodelabs[i] > 0.95) {
    nodelabels("", i + Ntip(td), frame = "none", pch = 16, cex = 0.7)
  } else {
    nodelabels("", i + Ntip(td), frame = "none", pch = 16, 
               cex = 2.5, col = "lightblue") 
    nodelabels(nodelabs[i], i + Ntip(td), frame = "none", cex = 0.5) 
  }
}
dev.off()

###################################
# LTT
###################################

tdno = drop.tip(td, c("Xylorhiza_tortifolia",
                      "Enceliopsis_covillei",
                      "Encelia_frutescens_glandulosa",
                      "Encelia_farinosa_phenicodonta"))
tdno1 = drop.tip(tdno, c("Encelia_californica2", "Encelia_virginensis1"))

ltt1 = ltt(tdno)

b<-phytools:::qb(tdno)
lttexp = data.frame(times = ltt1$times,
                    actual = ltt1$ltt,
                    `pure birth` = 2*exp(b * ltt1$times))
lttexp = lttexp %>% tidyr::gather(key = "exp_actual", value = "lineage", -times)
xx = ggplot(lttexp, aes(times, lineage)) +
  geom_line(aes(color = exp_actual)) + scale_y_log10() +
  scale_color_manual(values = c("black", "red")) +
  xlab("time") + ylab("lineages") + theme(legend.title = element_blank())
save_plot("LTT.png", xx, base_height = 3, base_width = 4)


###################################
# div rate
###################################

# + 3 for conspersa, halimifolia, hispida
n = Ntip(tdno) + 3
bd.ms(tdno, missing = 3, epsilon = 0.0)
bd.ms(tdno, missing = 3, epsilon = 0.1)
bd.ms(tdno, missing = 3, epsilon = 0.3)
bd.ms(tdno, missing = 3, epsilon = 0.5)
bd.ms(tdno, missing = 3, epsilon = 0.9)

###################################
# SNAQ tree
###################################

t = read.tree("~/Dropbox/Encelia/analysis/phylogeny/SNAQ/run1.OUT3.SNAQ.tre")
t = root(t, "Enceliopsis_covillei")
t$edge.length = NULL

t$tip.label = tipnames = gsub("_", " ", t$tip.label)
plt_t = ggtree(td) + geom_tiplab(fontface = "italic") + xlim(0, 6.5) +
  geom_taxalink('Encelia californica1', 'Encelia asperifolia',
                color='red', linetype = 'dashed', curvature = -0.5)
save_plot("SNAQ_network.pdf", plt_t)

###################################
# individual tree
###################################

t = ggtree::read.tree("~/Dropbox/Encelia/analysis/phylogeny/concatenated/RAxML_bipartitions.boot")
outs = c("XYL", "ENC-1", "ENC-2")
t1 = root(t, outs, resolve.root = T)
t2 = drop.tip(t1, outs)
t3 = read.tree(text = write.tree(ladderize(t2)))

lins = unique(dplyr::pull(x[match(t3$tip.label, x$sample), "lineage"]))
lins2 = data.frame(lineage = x[match(t3$tip.label, x$sample), "lineage"],
                   tips = t3$tip.label, stringsAsFactors = F)
cols = rep(brewer.pal(6, "Set3"), 3)
tt = ggtree(t3)
for (i in 1:length(lins)) {
  tips = lins2[lins2$lineage == lins[i], "tips"]
  node = findMRCA(t3, tips, type = 'node')
  tt = tt + geom_hilight(node= node, fill = cols[i], alpha=0.5) + 
    geom_cladelabel(node, linnames[i], fontface = "italic", offset=0, 
                    barsize = NA, angle=0, offset.text=0.0005, 
                    align = T, size = 1)
}
tt = tt + xlim(0, 0.05) + 
  geom_point2(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) >= 95),
              size = 0.2) +
  geom_point2(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) < 95),
              size = 0.7, fill = "white", shape = 21) 
save_plot("Encelia_individual_phylogeny.concat_miss0.6.pdf", 
          tt, base_height = 5, base_width = 6)

x = 1
ts1 = gsub("^.*/", "", ts1)
res = data.frame(t1 = rep(NA, 3), t2 = rep(NA, 3), rf = rep(NA, 3))
for (i in 1:2) {
  for (j in (i + 1):3) {
    res[x, "t1"] = ts1[[i]]
    res[x, "t2"] = ts1[[j]]
    res[x, "rf"] = RF.dist(ts3[[i]], ts3[[j]], normalize = T)
    x = x + 1
  }
} 

ts1 = list.files("~/Dropbox/Encelia/analysis/phylogeny/concatenated/", pattern = "best", full.names = T)
ts2 = lapply(ts1, read.tree)
ts3 = lapply(ts2, unroot)

############
# make astral mapping file
############

x1 = x[x$sample %in% t2$tip.label, ]
lins = unique(x1$lineage)
for (i in 1:length(lins)) {
  tinds = paste(pull(x1[x1$lineage == lins[i], "sample"]), collapse=",")
  cat(paste(lins[i], ":", tinds, "\n", sep=""))
}

############
# plot astral
############

t = read.tree("/Users/Sonal/Dropbox/Encelia/analysis/phylogeny/astral/astral.miss0.6.tre")
t1 = root(t, c("Xylorhiza_tortifolia", "Enceliopsis_covillei"))
t2 = read.tree(text = write.tree(ladderize(t1)))

pdf("~/Desktop/Encelia_phylogeny.astral.pdf", height = 5, width = 8)
par(xpd = T, mar = c(0.5, 0.5, 0.5, 12))
plot(t2, show.tip.label = F)
tiplabels(gsub("_", " ", t2$tip.label), 
          cex = 1.2, frame = "none", adj = 0, font = 3)
nodelabels("", frame = "none", pch = 16, col="skyblue", cex = 3)
nodelabels(t2$node.label, frame = "none", col="black", cex = 0.7)
dev.off()

ts1 = list.files("~/Dropbox/Encelia/analysis/phylogeny/astral/", pattern = "", full.names = T)
ts2 = lapply(ts1, read.tree)
ts3 = lapply(ts2, unroot)

x = 1
ts1 = gsub("^.*/", "", ts1)
res = data.frame(t1 = rep(NA, 3), t2 = rep(NA, 3), rf = rep(NA, 3))
for (i in 1:2) {
  for (j in (i + 1):3) {
    res[x, "t1"] = ts1[[i]]
    res[x, "t2"] = ts1[[j]]
    res[x, "rf"] = RF.dist(ts3[[i]], ts3[[j]], normalize = T)
    x = x + 1
  }
} 


##########
# plot uncertainty
##########

chisq = function(DF1, DF2, N){
  tryCatch({
    # converts percentages to counts, runs chisq, gets pvalue
    chisq.test(c(round(DF1*N)/100, round(DF2*N)/100))$p.value
  },
  error = function(err) {
    # errors come if you give chisq two zeros
    # but here we're sure that there's no difference
    return(1.0)
  })
}

t = read.tree("~/Dropbox/Encelia/analysis/phylogeny/concatenated/RAxML_bipartitions.boot")
g = read.table("~/Dropbox/Encelia/analysis/phylogeny/concordance/concord.cf.stat", 
               header = T, 
               skip = 12, sep = "\t")
s = read.table("~/Dropbox/Encelia/analysis/phylogeny/concordance/scf.cf.stat", header = T, 
               skip = 11, 
               sep = "\t")
cf = inner_join(g, s)
cf$bootstrap = as.numeric(t$node.label)[2:70]
cf2 = cf %>% 
  group_by(ID) %>%
  mutate(gEF_p = chisq(gDF1, gDF2, gN)) %>%
  mutate(sEF_p = chisq(sDF1, sDF2, sN))

a = ggplot(cf2, aes(x = gCF, y = sCF)) + 
  geom_point(aes(colour = bootstrap)) + 
  scale_colour_viridis(direction = -1) 
b = ggplot(cf2, aes(x = Length, y = sCF)) + 
  geom_point(aes(colour = gCF)) + xlab("branch len.") +
  scale_colour_viridis(direction = -1) 

t$node.label = as.numeric(t$node.label)
outs = c("XYL", "ENC-1", "ENC-2")
t1 = root(t, outs, resolve.root = T)
t2 = drop.tip(t1, outs)
tips2 = pull(x[match(t2$tip.label, x$sample), "lineage"])
t2$tip.label = gsub("_", " ", tips2)

plot_tree <- function() {
  par(xpd = T, mar = c(0, 0, 0, 6))
  plot.phylo(t2, cex = 0.5, show.tip.label = F)
  tiplabels(t2$tip.label, 
            cex = 0.8, frame = "none", adj = 0, font = 3)
  nodes = as.numeric(t2$node.label)
  for (i in 1:length(nodes)) {
    if (nodes[i] >= 95) {
      nodelabels("", i + Ntip(t2), pch = 16, col = "black",
                 frame = "none", cex  = 0.5)
    } else {
      nodelabels(nodes[i], i + Ntip(t2), col = "darkblue", 
                 frame = "none", adj = c(1, 0), cex = 0.5)
    }
  }
}
plot_tree()



# 
# c = ggtree(t2) + geom_point2(aes(subset=!isTip, 
#                              fill=as.numeric(label)),
#                          shape=21, size=3) +  
#   geom_tiplab(fontface = "italic", size  = 2.5) +
#   scale_fill_viridis(direction = -1, name = "bootstrap") + xlim(0, 0.035)
# layout = "
# AB
# AC
# "
# abc = ~plot_tree + a + b + plot_layout(design = layout, widths = c(2.2, 1))

ab = plot_grid(a, b, ncol = 1)
abc = plot_grid(plot_tree, ab)
save_plot("~/Dropbox/Encelia/manuscripts/Encelia_Phylogeny/figures/Encelia_individual_phylogeny.concordance.png", abc, base_width = 10, base_height = 6)

############
# compare species trees
############

compare_trees <- function(t1, t2) {
  t1tips = t1$tip.label
  t2tips = t2$tip.label
  
  nodelabs = rep(NA, t1$Nnode)
  for (i in 1:t1$Nnode) {
    nodenum = i + length(t1tips)
    alldesc = t1tips[Descendants(t1, nodenum, "tips")[[1]]]
    alldesc2 = alldesc[alldesc %in% t2tips]
    if (length(alldesc2) > 1) {
      intnode = getMRCA(t2, alldesc2)
      intdesc = t2tips[Descendants(t2, intnode, "tips")[[1]]]
      intdesc2 = intdesc[intdesc %in% t1tips]
      if (setequal(intdesc2, alldesc2)) {
        nodelabs[i] = 1
      } else{
        nodelabs[i] = 0
      }
    } else {
      nodelabs[i] = 1
    }
  }
  return(nodelabs)
}

t1 = read.tree("~/Dropbox/Encelia/analysis/phylogeny/astral/astral.miss0.6.tre")
t1 = root(t1, c("Enceliopsis_covillei", "Xylorhiza_tortifolia"),
         resolve.root = T)
t1 = drop.tip(t1, c("Enceliopsis_covillei", "Xylorhiza_tortifolia"))

t2 = read.nexus("~/Dropbox/Encelia/analysis/phylogeny/SVDquartets/encelia.miss0.5.MAC2.thinned.species.tre")
t2 = root(t2, c("Enceliopsis_covillei"),
          resolve.root = T)
t2 = drop.tip(t2, c("Enceliopsis_covillei"))
t2 = minRotate(t2, setNames(1:Ntip(t1), t1$tip.label))
t2 = read.tree(text = write.tree(t2))

pdf("species.ASTRAL_vs_SVD.pdf", width = 10, height = 3.5)
par(mfrow=c(1, 2), mar = c(1, 1, 1, 8), xpd = T)
plot(t1, show.tip.label = F)
tiplabels(gsub("_", " ", t1$tip.label), 
          cex = 0.9, frame = "none", adj = 0, font = 3)
nodelabs1 = compare_trees(t1, t2)
nodelabs2 = compare_trees(t2, t1)
nodelabels("", bg = ifelse(nodelabs1 == 1, "black", "red"), 
           pch = 21, cex = ifelse(nodelabs1 == 1, 0.2, 1), frame = "none")
plot(t2, show.tip.label = F)
tiplabels(gsub("_", " ", t2$tip.label), 
          cex = 0.9, frame = "none", adj = 0, font = 3)
nodelabels("", bg = ifelse(nodelabs2 == 1, "black", "red"), 
           pch = 21, cex = ifelse(nodelabs2 == 1, 0.2, 1), frame = "none")
dev.off()

############
# compare individual trees
############

t1 = read.tree("~/Dropbox/Encelia/analysis/phylogeny/concatenated/RAxML_bestTree.concatenated0.6")
outs = c("XYL", "ENC-1", "ENC-2")
t1 = root(t1, outs,
          resolve.root = T)
t1 = drop.tip(t1, outs)

t2 = read.nexus("~/Dropbox/Encelia/analysis/phylogeny/SVDquartets/encelia.miss0.5.MAC2.thinned.tre")
t2$tip.label = gsub("_", "-", t2$tip.label)
outs = c("ENC-1", "ENC-2")
t2 = root(t2, outs,
          resolve.root = T)
t2 = drop.tip(t2, outs)
t2 = drop.tip(t2, c("RAV-5", "RAV-4",
                    "RAV-6", "RAV-7",
                    "RAV-8", "RAV-9",
                    "RAV-3", "DEN-13-C",
                    "DEN-18-B", "DEN-12-B",
                    "DEN-13-D", "DEN-13-E",
                    "DEN-12-E", "DEN-13-F",
                    "DEN-12-C", "DEN-12-D",
                    "FAR-PHE-60-1", "FAR-PHE-71-2",
                    "FAR-PHE-73-2"))

t1 = read.tree(text = write.tree(ladderize(t1)))
t2 = minRotate(t2, setNames(1:Ntip(t1), t1$tip.label))

png("individual.concat_vs_SVD.png", width = 10, 
    height = 8, units = "in", res = 200)
par(mfrow=c(1, 2), mar = c(1, 1, 1, 10), xpd = T)
plot(t1, show.tip.label = F)
tips1 = dplyr::pull(x[match(t1$tip.label, x$sample), "lineage"])
tips1 = paste(gsub("Encelia_", "", tips1), t1$tip.label)
tiplabels(gsub("_", " ", tips1), 
          cex = 0.7, frame = "none", adj = 0, font = 3)
nodelabs1 = compare_trees(t1, t2)
nodelabs2 = compare_trees(t2, t1)
nodelabels("", bg = ifelse(nodelabs1 == 1, "black", "red"), 
           pch = 21, cex = ifelse(nodelabs1 == 1, 0.2, 1), frame = "none")
plot(t2, show.tip.label = F)
tips2 = dplyr::pull(x[match(t2$tip.label, x$sample), "lineage"])
tips2 = paste(gsub("Encelia_", "", tips2), t2$tip.label)
tiplabels(gsub("_", " ", tips2), 
          cex = 0.7, frame = "none", adj = 0, font = 3)
nodelabels("", bg = ifelse(nodelabs2 == 1, "black", "red"), 
           pch = 21, cex = ifelse(nodelabs2 == 1, 0.2, 1), frame = "none")
dev.off()
