library(ggplot2)
library(cowplot)
library(rnaturalearth)
library("rnaturalearthdata")
theme_set(theme_cowplot())
library(readr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggtree)

world <- ne_countries(scale = "small", returnclass = "sf")

# sampled points
spoints = read_csv("~/Dropbox/Encelia/analysis/spatial_analyses/georef_samples.csv")
x = read_csv("~/Dropbox/Encelia/ddRAD/analysis/encelia_samples_v4.csv")
spoints = inner_join(spoints, x, by = c("sample" = "sample"))

# pts
pts = read_csv("~/Dropbox/Encelia/analysis/spatial_analyses/encelia/all_points_thinned.csv")
pts2 = pts[grep("Enc", pts$species), ]
pts3 = pts2[which(pts2$species == "Encelia canescens"), ]
pts4 = pts2[which(pts2$species != "Encelia canescens"), ]


########################
# tree plot
########################

t = ggtree::read.tree("~/Dropbox/Encelia/analysis/phylogeny/concatenated/RAxML_bipartitions.boot")
outs = c("XYL", "ENC-1", "ENC-2")
t1 = root(t, outs, resolve.root = T)
t2 = drop.tip(t1, outs)
t3 = read.tree(text = write.tree(ladderize(t2)))

lins = unique(dplyr::pull(x[match(t3$tip.label, x$sample), "lineage"]))
lins2 = data.frame(lineage = x[match(t3$tip.label, x$sample), "lineage"],
                   tips = t3$tip.label, stringsAsFactors = F)
getPalette = colorRampPalette(brewer.pal(12, "Set3"))
cols = getPalette(length(lins))
names(cols) = sort(lins)

linnames = gsub("Encelia_", "E. ", lins)
linnames = gsub("_", " ", linnames)
tt = ggtree(t3)
for (i in 1:length(lins)) {
  tips = lins2[lins2$lineage == lins[i], "tips"]
  node = findMRCA(t3, tips, type = 'node')
  tt = tt + geom_hilight(node= node, fill = cols[lins[i]], alpha=0.5) + 
    geom_cladelabel(node, linnames[i], fontface = "italic", offset=0, 
                    barsize = NA, angle=0, offset.text=0.0005, 
                    align = T, size = 1)
}
tt = tt + xlim(0, 0.05) + 
  geom_point2(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) >= 95),
              size = 0.2) +
  geom_point2(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) < 95),
              size = 0.7, fill = "white", shape = 21) 


spmaps = lins
spmaps = spmaps[!spmaps %in% c("Encelia_farinosa_phenicodonta",
                     "Encelia_californica2",
                     "Encelia_virginensis2",
                     "Encelia_frutescens_glandulosa")]
spmaps = gsub("\\d+", "", spmaps)
spmaps = gsub("_", " ", spmaps)
spmaps[which(spmaps == "Encelia actoni")] = "Encelia actoni"
spmaps[which(spmaps == "Encelia frutescens frutescens")] = "Encelia frutescens"
spmaps[which(spmaps == "Encelia farinosa farinosa")] = "Encelia farinosa"
spplots = vector("list", length(spmaps))

tips2n = gsub("Encelia ", "E. ", spmaps)
spoints[spoints$lineage.x == "Encelia actonii", "lineage.x"] = "Encelia actoni"
pts2[pts2$species == "Encelia actonii", "species"] = "Encelia actoni"
ptsalpha = 1 - log(table(pts2$species)) / 8
names(cols) = gsub("_", " ", names(cols))
spoints$lineage.y = gsub("_", " ", spoints$lineage.y)

for (i in 1:length(spmaps)) {
  if (spmaps[i] == "Encelia canescens") {
    xlim = c(-80, -50)
    ylim = c(-40, 10)
  } else {
    xlim = c(-120.51, -108.7)
    ylim = c(24.33, 37.9)
  }

  sub = spoints %>% filter(lineage.x == spmaps[i])
  
   spplots[[i]] = ggplot(data = world) + 
    geom_sf(color = "gray80", fill = "gray80") + 
    xlim(xlim) +
    ylim(ylim) + 
    geom_point(data = pts2 %>% filter(species == spmaps[i]),
               aes(Longitude, Latitude), size = 0.5, 
               alpha = ptsalpha[spmaps[i]]) +
    geom_point(data = sub,
               aes(LONGITUDE, LATITUDE, fill = lineage.y), 
               size = 1.8, shape = 21) + ggtitle(tips2n[i]) +
              scale_fill_manual(values = cols[unique(sub$lineage.y)]) +
     
    theme_void() +
    theme(plot.title = element_text(size=10, face="italic", hjust = 0.5),
          legend.position = "none")
}

# special ones
# farinosa, californica, virginensis, frutescens


layout <- "
ADBC
AEF#
AGHI
AJK#
ALM#
"
png("~/Dropbox/Encelia/manuscripts/Encelia_Phylogeny/figures/Encelia_maps_gray.png", width = 8, height = 6, units = "in", res = 200)
tt + spplots + plot_layout(design = layout, widths = c(4, 1, 1, 1))
dev.off()

