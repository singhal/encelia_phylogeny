library(rgeos)
library(rgdal)
library(raster)
library(maptools)
library(readr)
library(RColorBrewer)
library(ape)
library(ggspatial)
library(rworldmap)

# get map
worldmap <- getMap(resolution = "coarse")

# sampled points
spoints = read_csv("~/Desktop/georef_samples.csv")
x = read_csv("~/Dropbox/Encelia/ddRAD/analysis/encelia_samples_v4.csv")
spoints = inner_join(spoints, x, by = c("sample" = "sample"))

# pts
pts = read_csv("~/Dropbox/Encelia/analysis/spatial_analyses/encelia/all_points_thinned.csv")
pts2 = pts[grep("Enc", pts$species), ]
pts3 = pts2[which(pts2$species == "Encelia canescens"), ]
pts4 = pts2[which(pts2$species != "Encelia canescens"), ]

# shp = readOGR("~/Dropbox/Encelia/spatial_rasters/NA_CEC_Eco_Level3/NA_CEC_Eco_Level3.shp")
# 
# # project it
# shp2 = spTransform(shp, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
# 
# # warm deserts
# mojave = shp2[shp2$NA_L3CODE == "10.2.1", ]
# mojave2 = gSimplify(mojave, tol = 0.01, topologyPreserve = T)
# son = shp2[shp2$NA_L3CODE == "10.2.2", ]
# son2 = gSimplify(son, tol = 0.01, topologyPreserve = T)
# baja = shp2[shp2$NA_L3CODE == "10.2.3", ]
# baja2 = gSimplify(baja, tol = 0.01, topologyPreserve = T)
# gila = shp2[shp2$NA_L3CODE == "13.1.1", ]
# gila2 = gSimplify(gila, tol = 0.01, topologyPreserve = T)
# # med
# med = shp2[grep("^11.1.", shp2$NA_L3CODE), ]
# med2 = gSimplify(med, tol = 0.01, topologyPreserve = T)
# # cold des
# colddes = shp2[grep("^10.1.", shp2$NA_L3CODE), ]
# colddes2 = gSimplify(colddes, tol = 0.01, topologyPreserve = T)
# 
# cols = brewer.pal(6, "Set3")
# shps = list(mojave2, son2, baja2, gila2, colddes2, med2)
# names(shps) = c("Mojave Desert", "Sonoran Desert", "Baja California Desert",
#                 "Upper Gila Mountains", "Great Basin Desert", "Mediterranean")

# prep tree
t = read.tree("~/Dropbox/Encelia/analysis/phylogeny/astral/astral.miss0.6.brlen.tre")
t1 = read.tree("~/Dropbox/Encelia/analysis/phylogeny/astral/astral.miss0.6.tre")
t1 = minRotate(t1, setNames(1:Ntip(t), t$tip.label))
t$node.label = as.numeric(t1$node.label)
t = root(t, c("Enceliopsis_covillei", "Xylorhiza_tortifolia"),
         resolve.root = T)
# our divergence dating says PAL - VEN = 1 myr
# smith et al tree says root height is 1.36 myr
rootenc1 = getMRCA(t, t$tip.label[grep("Encelia_", t$tip.label)])
td = chronopl(t, 0.0001, node = rootenc1, age.min = 1.36)
td = read.tree(text = write.tree(ladderize(td)))
td_nom = drop.tip(td, c("Encelia_californica1", "Encelia_virginensis2",
                        "Encelia_frutescens_glandulosa"))
td_nom$tip.label = gsub("\\d", "", td_nom$tip.label)
t5 = root(td_nom, c("Enceliopsis_covillei", "Xylorhiza_tortifolia"),
              resolve.root = T)
t5 = drop.tip(t5, c("Enceliopsis_covillei", "Xylorhiza_tortifolia",
                    "Encelia_farinosa_phenicodonta"))
tips = gsub("_", " ", t5$tip.label)
tips[tips == "Encelia farinosa farinosa"] = "Encelia farinosa"
tips[tips == "Encelia frutescens frutescens"] = "Encelia frutescens"
tips[tips == "Encelia actoni"] = "Encelia actonii"

# ggplot() + geom_raster(data = chi2, aes(x = x, y = y))

tips2 = rev(tips)
tips2 = c(tips2[1:2], "blank",
          tips2[3:5],
          tips2[6:7], "blank",
          tips2[8:10],
          tips2[11:12], "blank")
# pdf("~/Desktop/Encelia_maps.pdf", width = 8, height = 10)
# layout(matrix(data = c(1, 2, 3, 4,
#                   1, 5, 6, 7,
#                   1, 8, 9, 10,
#                   1, 11, 12, 13,
#                   1, 14, 15, 16,
#                   17, 17, 17, 17), 6, byrow = T), 
#          widths = c(0.4, 0.2, 0.2, 0.2),
#          heights = c(0.18, 0.18, 0.18, 0.18, 0.18, 0.1))
# #  matrix(c(rep(1, 12), seq(2, 13), rep(14, 12)), 3, byrow = T),
# #       heights = c(0.25, 0.55, 0.2))
# # layout.show(17)
# par(mar = c(0.5, 0.5, 2, 0.5))
# for (n in 1:length(tips2)) {
#   if (tips2[n] == "blank") {
#     plot(1, type="n", axes=F, xlab="", ylab="")
#   }
#   else if (tips2[n] != "Encelia canescens") {
#     plot(worldmap, col = "lightgrey", 
#          border = NA,
#          xlim = c(-121, -108.5), 
#          ylim = c(24, 38),
#          bg = "white",
#          asp = 1, main = tips2[n])
#     for (i in 1:length(shps)) {
#       plot(shps[[i]], col = cols[[i]], border = NA, add = T)
#      }
#     sppts = pts2[pts2$species == tips2[n], ]
#     points(sppts$Longitude, sppts$Latitude, pch = 16, cex = 0.5)
#     sppts2 = spoints %>% filter(lineage.x == tips2[n])
#     points(sppts2$Longitude, sppts2$Latitude, pch = 16, col = "red")
#   } else {
#     plot(worldmap, col = "lightgrey", 
#          border = NA,
#          xlim = c(-80, -68), 
#          ylim = c(-36, -7),
#          bg = "white",
#          asp = 1, main = tips2[n])
#     sppts = pts2[pts2$species == tips2[n], ]
#     points(sppts$Longitude, sppts$Latitude, pch = 16, cex = 0.5)
#     sppts2 = spoints %>% filter(lineage.x == tips2[n])
#     points(sppts2$Longitude, sppts2$Latitude, pch = 16, col = "red")
#   }
# }
# par(mar = c(1, 1, 1, 1))
# plot(1, type = "n", axes=FALSE, xlab="", ylab="")
# legend(x = "top", inset = 0,
#        legend = names(shps),
#        fill = cols, cex = 1.5, ncol = 3, bty= "n")
# dev.off()

ptsalpha = 1 - log(table(pts$species)) / 10
tips2n = gsub("Encelia ", "E. ", tips2)

pdf("~/Dropbox/Encelia/manuscripts/Encelia_Phylogeny/figures/Encelia_maps_gray.pdf", width = 6, height = 6)
layout(matrix(data = c(1, 2, 3, 4,
                       1, 5, 6, 7,
                       1, 8, 9, 10,
                       1, 11, 12, 13,
                       1, 14, 15, 16,
                       1, 17, 17, 17), 6, byrow = T), 
       widths = c(0.55, 0.15, 0.15, 0.15), 
       heights = c(rep(0.19, 5), 0.05))
par(mar = c(2, 0.5, 1.5, 12), xpd = T)
tdo = drop.tip(td, c("Xylorhiza_tortifolia", "Enceliopsis_covillei"))
plot(tdo, show.tip.label = F)
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

tipnames = gsub("Encelia_", "E. ", tdo$tip.label)
tipnames = gsub("_", " ", tdo$tip.label)
tiplabels(tipnames, frame = "none", adj = 0, font = 3)
nodelabs = as.numeric(tdo$node.label)
nodelabs[!complete.cases(nodelabs)] = 1.0
for (i in 1:length(nodelabs)) {
  if (nodelabs[i] > 0.95) {
    nodelabels("", i + Ntip(tdo), frame = "none", pch = 16, cex = 0.7)
  } else {
    nodelabels("", i + Ntip(tdo), frame = "none", pch = 16, 
               cex = 2.5, col = "lightblue") 
    nodelabels(nodelabs[i], i + Ntip(tdo), frame = "none", cex = 0.5) 
  }
}
par(mar = c(0.5, 0.5, 2, 0.5), xpd = F)
## show the regions that have been allocated to each plot
for (n in 1:length(tips2)) {
  # plot(shp3)
  if (tips2[n] == "blank") {
    plot(1, type="n", axes=F, xlab="", ylab="")
  }
  else if (tips2[n] != "Encelia canescens") {
    plot(worldmap, col = "lightgrey", 
         border = NA,
         xlim = c(-120.51, -108.7), 
         ylim = c(24.33, 37.9),
         bg = "white",
         asp = 1, main = tips2n[n], font.main = 3,
         cex = 0.8)
    sppts = pts2[pts2$species == tips2[n], ]
    points(sppts$Longitude, sppts$Latitude, pch = 16, cex = 0.5,
           col = alpha("black", ptsalpha[tips2[n]]))
    sppts2 = spoints %>% filter(lineage.x == tips2[n])
    points(sppts2$LONGITUDE, sppts2$LATITUDE, pch = 16, col = "red")
  } else {
    plot(worldmap, col = "lightgrey", 
         border = NA,
         xlim = c(-80, -68), 
         ylim = c(-36, -7),
         bg = "white",
         asp = 1, main = tips2n[n], font.main = 3,
         cex = 0.8)
    sppts = pts2[pts2$species == tips2[n], ]
    points(sppts$Longitude, sppts$Latitude, pch = 16, cex = 0.5,
           col = alpha("black", ptsalpha[tips2[n]]))
    sppts2 = spoints %>% filter(lineage.x == tips2[n])
    points(sppts2$LONGITUDE, sppts2$LATITUDE, pch = 16, col = "red")
  }
}
dev.off()








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



ggplot(worldmap, col = "lightgrey") + geom_polygon() + xlim(-120.51, -108.7) 
ylim = c(24.33, 37.9)