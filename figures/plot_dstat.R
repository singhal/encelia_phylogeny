library(ape)
library(maptools)
library(phytools)
library(RColorBrewer)
library(scales)
library(viridis)

setwd("~/Dropbox/Encelia/manuscripts/Encelia_Phylogeny/figures/")

t = read.tree("~/Dropbox/Encelia/analysis/phylogeny/astral/astral.miss0.6.brlen.dated.tre")
t = drop.tip(t, "Xylorhiza_tortifolia" )
t = drop.tip(t, "Enceliopsis_covillei" )
t = read.tree(text = write.tree(ladderize(t)))
bts = branching.times(t)
names(bts) = seq(1 + Ntip(t), Ntip(t) + Nnode(t))

d = read.csv("~/Dropbox/Encelia/analysis/pop_gen/dstat/dvals.csv", stringsAsFactors = F)
d$X = NULL
names(d) = c("sp1", "sp2", "sp3", "out", "num_sites", "D", "d_sd", "d_z")
d$p = 2*pnorm(-abs(as.numeric(d$d_z)))

# only look at ones using covillei as outgroup
# because it is the true outgroup
d = d[d$out == 'Enceliopsis_covillei', ]

d$focal_sp = ifelse(d$D < 0, d$sp1, d$sp2)

sps = unique(c(d$sp1, d$sp2, d$sp3))
res = data.frame(matrix(NA, sum(seq(1, length(sps))), 6))
names(res) = c("sp1", "sp2", "sp3", "Dstat", "z", "p")
track = 1
for (i in 1:length(sps)) {
  for (j in (i + 1):length(sps)) {
    # for every two species
    sp1 = sps[i]
    sp2 = sps[j]
    
    res[track, "sp1"] = sp1
    res[track, "sp2"] = sp2
    
    d1 = d[which(d$sp3 == sp2), ]
    d1 = d1[which(d1$focal_sp == sp1), ]
    
    d2 = d[which(d$sp3 == sp1), ]
    d2 = d2[which(d2$focal_sp == sp2), ]
    xx = rbind(d1, d2)
    if (nrow(xx) > 0){
      # get the MRCA of the comparisons
      xx$bt = apply(xx, 1, function(x) { bts[as.character(findMRCA(t, as.character(x[1:3])))]})
      xx = xx[complete.cases(xx$D), ]
      xx = xx[with(xx, order(bt, -abs(D))), ]
      # take the min
      # very conservative approach
      min_bt = xx[xx$bt == min(xx$bt), ]
      min_bt2 = min_bt[order(min_bt$d_z), ]
      res[track, "Dstat"] = abs(min_bt2[1, "D"])
      res[track, "z"] = min_bt2[1, "d_z"]
      res[track, "p"] = min_bt2[1, "p"]
      
      # this is the bystander species!!!
      poss_sps = as.character(xx[1, 1:2])
      res[track, "sp3"] = poss_sps[which(poss_sps != xx[1, "focal_sp"])]
      
      # take the max
      # res[track, "Dstat"] = abs(xx[1, "D"])
      # res[track, "z"] = xx[1, "d_z"]
      # res[track, "p"] = xx[1, "p"]
      # take the mean -- should this be absolute - yes
      # but this is also statistically incomprehensible
      # res[track, "Dstat"] = mean(xx$D, na.rm = T)
      # res[track, "z"] =  mean(xx$d_z, na.rm =T)
      # res[track, "p"] =  mean(xx$p, na.rm = T)
    }
    track = track + 1
  }
}

res = res[complete.cases(res$sp1), ]
res = res[complete.cases(res$sp2), ]
write.csv(res, "~/Dropbox/Encelia/manuscripts/Encelia_Phylogeny/figures/dstat.csv", row.names = F, quote = F)
# View(res)

tips = t$tip.label
m = matrix(NA, nrow=length(tips), ncol=length(tips), dimnames=list(tips, tips))
for (i in 1:length(tips)) {
  for (j in 1:length(tips)) {
    sp1 = tips[i]
    sp2 = tips[j]
    res2 = res[(res$sp1 == sp1 & res$sp2 == sp2) | (res$sp1 == sp2 & res$sp2 == sp1), ]
    if (nrow(res2) > 0) {
      sps = c(sp1, sp2)
      sps = sps[order(match(sps, tips))]
      m[sps[2], sps[1]] = res2$p
    }
  }
}

range(res$p, na.rm = T)
range(log(m), na.rm = T)
range(m, na.rm = T)
# vals = seq(1.36, 3.6, length = 99)
# vals = c(-3.2, vals)
# breaks = -1 * qnorm(exp(vals) / 2)
# breaks = rev(seq(0, 40, by = 0.5))
# breaks = exp(vals)
breaks = c(1e-300, 1e-200, 1e-100, 1e-50, 1e-10, 1e-5, 5e-4, 1)

# newcol <- colorRampPalette(brewer.pal(9, 'Greys'))
# col <- newcol(length(breaks) - 1)
col = viridis(length(breaks) - 1)

t$tip.label = gsub("Encelia_", "E. ", t$tip.label)
# t$tip.label = gsub("_\\S+", "", t$tip.label)
png("Dstat.png", height = 7, width = 10, units = "in", res = 200)
op <- par(no.readonly = TRUE)
figs <- matrix(c(0, 0.2, 0.9, 1, 
                 0.4, 0.9, 0.6, 1, # phylo 1
                 0, 0.4, 0, 0.6, # phylo2
                 0.4, 0.9, 0, 0.6,  # grid
                 0.9, 1, 0.15, 0.85), byrow = TRUE, nrow = 5, ncol = 4)
bt = max(branching.times(t)) 
par(fig = figs[2, ], new = FALSE, mar = c(6, 0, 1, 0), xpd = T)
plot.phylo(t, direction = "downwards", 
           show.tip.label = TRUE, 
           x.lim = c(0.5, length(t$tip.label)),
           y.lim = c(0, bt), cex = 0.7)
par(fig = figs[3, ], new = TRUE, mar = c(0, 1, 0, 6), xpd = T)
plot.phylo(t, direction = "rightwards", 
           show.tip.label = TRUE, 
           y.lim = c(0.5, length(t$tip.label)), 
           x.lim = c(0, bt), cex = 0.7)
par(fig = figs[4, ], new = TRUE, mar = c(0, 0, 0, 0), xpd = F)
gl <- 1:(length(t$tip.label) + 1)
plot(0, 0, type = "n", 
     axes = FALSE, ann = FALSE, 
     xlim = c(1, length(gl) - 0.5 ), 
     ylim = c(1, length(gl) - 0.5 ))
image(gl, gl, m, axes = FALSE, xlab = "",
      ylab = "", col = col, 
      xlim = c(1, length(gl) - 1), 
      ylim = c(1, length(gl) - 1), breaks = breaks, add = TRUE)
# pvals = c(1e-200, 1e-100, 1e-10, 1e-5)
# locs = qnorm(pvals / 2) * -1 
# locs = sapply(locs, function(x) {which(abs(breaks-x)==min(abs(breaks-x)))}) / length(breaks)
pvals = breaks
locs = seq(0, 1, 1 / (length(breaks) - 1))

par(fig = figs[5, ], new = TRUE, mar = c(2, 2, 2, 2), xpd = F)

barLegend <- function(pal, locs, labels, fig, side, mar = rep(0,4), colpalette = NULL, ...) {
  n <- length(pal);
  x <- seq(0, n, 1) / n;
  x <- rep(x, each = 2);
  x <- x[-c(1, length(x))];
  x <- matrix(x, ncol = 2, byrow = TRUE);
  # par(fig = fig, mar = mar, new = TRUE);
  plot.new();
  
  if (side == 2 || side == 4) {
    xlim <- c(-0.1, 0.1);
    ylim <- c(0, 1);
    plot.window(xlim, ylim);
    segments(x0 = 0, y0 = x[,1], x1 = 0, y1 = x[,2], col = pal, lwd = 8, lend = 2);
  }
  
  else {
    xlim <- c(0,1);
    ylim <- c(-0.1, 0.1);
    plot.window(xlim, ylim);
    segments(x0 = x[,1], y0 = 0, x1 = x[,2], y1 = 0, col = pal, lwd = 8, lend = 2);
  }
  mtext(labels, side, line = 0.001, at = locs, las = 1, cex = 0.7)
  # axis(side, at = locs, labels = labels, las=1, ...);
}

barLegend(col, locs, pvals, 
          fig = figs[5, ], side = 4, cex.axis=0.8, tck=-0.1)
par(op)
dev.off()