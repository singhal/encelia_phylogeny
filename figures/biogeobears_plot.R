library(rgeos)
library(rgdal)
library(raster)
library(maptools)
library(readr)
library(RColorBrewer)
library(ape)
library(rworldmap)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)

# get map
worldmap <- getMap(resolution = "coarse")

shp = readOGR("~/Dropbox/Encelia/spatial_rasters/NA_CEC_Eco_Level3/NA_CEC_Eco_Level3.shp")

## project it
shp2 = spTransform(shp, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

## warm deserts
mojave = shp2[shp2$NA_L3CODE == "10.2.1", ]
mojave2 = gSimplify(mojave, tol = 0.01, topologyPreserve = T)
son = shp2[shp2$NA_L3CODE == "10.2.2", ]
son2 = gSimplify(son, tol = 0.01, topologyPreserve = T)
baja = shp2[shp2$NA_L3CODE == "10.2.3", ]
baja2 = gSimplify(baja, tol = 0.01, topologyPreserve = T)
gila = shp2[shp2$NA_L3CODE == "13.1.1", ]
gila2 = gSimplify(gila, tol = 0.01, topologyPreserve = T)
## med
med = shp2[grep("^11.1.", shp2$NA_L3CODE), ]
med2 = gSimplify(med, tol = 0.01, topologyPreserve = T)
## cold des
colddes = shp2[grep("^10.1.", shp2$NA_L3CODE), ]
colddes2 = gSimplify(colddes, tol = 0.01, topologyPreserve = T)

cols = brewer.pal(6, "Set3")
shps = list(baja2, med2, mojave2, son2, colddes2)
geonames =  c("Baja California Desert", "Mediterranean",
              "Mojave Desert", "Sonoran Desert",
              "cold deserts", "Peru")
names(shps) = geonames[1:length(shps)]


###################
# prepare biogeobears object
###################
wd = np("~/Dropbox/Encelia/analysis/phylogeny/BioGeoBEARS/")

trfn = "~/Dropbox/Encelia/analysis/phylogeny/BioGeoBEARS/nominal.tre"
resfn = "~/Dropbox/Encelia/analysis/phylogeny/BioGeoBEARS/encelia_DEC.Rdata"
geogfn = "~/Dropbox/Encelia/analysis/phylogeny/BioGeoBEARS/regions.txt"
load(resfn)
results_object = res

# States
td_nom = read.tree(trfn)
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
class(td_nom) = "phylo"

states = results_object$ML_marginal_prob_each_state_at_branch_top_AT_node
maxstate = apply(states, 1, function(x) { which(x == max(x))})
geos = results_object$inputs$all_geog_states_list_usually_inferred_from_areas_maxareas
maxgeo = geos[maxstate]

pievals = matrix(0, nrow = length(maxstate), ncol = length(geonames))
for (i in 1:length(maxstate)) {
  curgeo = maxgeo[[i]] + 1
  per = 100 / length(curgeo)
  pievals[i, curgeo] = per
}

#############################
# plot
############################

png("~/Dropbox/Encelia/manuscripts/Encelia_Phylogeny/figures/biogeobears.png", height = 4, width = 8, units = "in", res = 200)
layout(matrix(c(1,1,2,3), 2, 2, byrow = FALSE), 
       heights = c(0.8, 0.2), widths = c(0.65, 0.35))
par(mar = c(0, 0, 0, 10), xpd = T)

tipnames = gsub("Encelia_", "E. ", td_nom$tip.label)
tipnames = gsub("_", " ", tipnames)
plot.phylo(td_nom, show.tip.label = F)
tiplabels(tipnames, frame = "none", adj = 0, font = 3, offset = 0.1)
tiplabels(pie = pievals[1:Ntip(td_nom), ], piecol=cols, 
          frame = "none", cex = 1.7)
nodelabels(pie = pievals[(Ntip(td_nom) + 1):nrow(pievals), ], 
           piecol=cols, cex = 1.7)

par(mar = c(1, 1, 1, 1), xpd = F)
plot(worldmap, col = "lightgrey", 
     border = NA,
     xlim = c(-120.51, -108.7), 
     ylim = c(24.33, 37.9),
     bg = "white",
     asp = 1, 
     cex = 0.8)
for (i in 1:length(shps)) {
  plot(shps[[i]], col = cols[i], border = NA,
       add = T)
}
plot.new()
par(mar = c(1, 1, 1, 1), xpd = T)
legend("center", legend = geonames, fill = cols,
       ncol = 2, bty = "n", cex = 0.8)
dev.off()


###########################
# make distribution for anc node
###########################

ancnode = findMRCA(td_nom, c("Encelia_densifolia", "Encelia_californica"))
probs = states[ancnode, ]
top10 = sort(probs, index.return=TRUE, decreasing=TRUE)$ix[1:7]

geo10 = geos[ top10 ]
probs2 = probs[ top10 ]

bar = c()
geo = c()
prob = c()

for (i in 1:length(probs2)) {
  tot = probs2[i]
  totdiv = probs2[i] / length(geo10[[i]])
  vals2 = rep(0, length(geonames))
  vals2[geo10[[i]] + 1] = totdiv
  
  bar = c(bar, rep(i, length(geonames2)))
  prob = c(prob, vals2)
  geo = c(geo, geonames)
}

vals = data.frame(bar = bar,
                  geo = geo,
                  prob = prob, 
                  stringsAsFactors = F)

names(cols) = geonames
aa = ggplot(data=vals, aes(x=bar, y=prob, fill=geo)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = cols[order(names(cols))]) + 
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("probability")



###########################
# make distribution for cold node
###########################

ancnode = findMRCA(td_nom, c("Encelia_actoni", "Encelia_resinifera"))
probs = states[ancnode, ]
top10 = sort(probs, index.return=TRUE, decreasing=TRUE)$ix[1:4]

geo10 = geos[ top10 ]
probs2 = probs[ top10 ]

bar = c()
geo = c()
prob = c()

for (i in 1:length(probs2)) {
  tot = probs2[i]
  totdiv = probs2[i] / length(geo10[[i]])
  vals2 = rep(0, length(geonames))
  vals2[geo10[[i]] + 1] = totdiv
  
  bar = c(bar, rep(i, length(geonames2)))
  prob = c(prob, vals2)
  geo = c(geo, geonames)
}

vals = data.frame(bar = bar,
                  geo = geo,
                  prob = prob, 
                  stringsAsFactors = F)


bb = ggplot(data=vals, aes(x=bar, y=prob, fill=geo)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = cols[order(names(cols))]) + 
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("probability")

ab = plot_grid(aa, bb)
save_plot("~/Dropbox/Encelia/manuscripts/Encelia_Phylogeny/figures/biogeobears2.png", ab, base_height = 3, base_width = 4, ncol = 2)
