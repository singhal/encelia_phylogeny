library(ape)
library(raster)
library(pcaMethods)
library(ggplot2)
library(tidyr)
library(dplyr)
library(rgeos)
library(readr)
library(cowplot)
theme_set(theme_cowplot())

####################
# get phylogeny
####################

# prep tree
t = read.tree("~/Dropbox/Encelia/analysis/phylogeny/astral/astral.miss0.6.brlen.dated.tre")
td = read.tree(text = write.tree(ladderize(t)))
td_nom = drop.tip(td, c("Encelia_californica1", "Encelia_virginensis2",
                        "Encelia_frutescens_glandulosa"))
td_nom$tip.label = gsub("\\d", "", td_nom$tip.label)

tips = td_nom$tip.label
tips = tips[! tips %in% c("Xylorhiza_tortifolia", "Enceliopsis_covillei")]
com = combn(tips, 2)

res = data.frame(species1 = com[1, ],
                 species2 = com[2, ],
                 phy_dist = NA,
                 morph_dist = NA,
                 clim_dist = NA,
                 soil_dist = NA,
                 geo_dist = NA,
                 per_overlap = NA,
                 stringsAsFactors = F)

###################
# calc phy dist
###################

phydist = cophenetic(td_nom)
for (i in 1:nrow(res)) {
  res[i, "phy_dist"] = phydist[res[i, "species1"], res[i, "species2"]]
}

##########################
# calc morph dist
##########################

m = read_csv("~/Dropbox/Encelia/analysis/phylogeny/morphological_data/morphological_data_spread.csv")
keep1 = c("BTSA", "top_trichomes_per_mm2",
          "Kshoot_mmol/s/m2/MPa", "leaf_gray_mean",
          "leaf_Area", 
          "leaf_Round", "leaf_LMA",
          "wood_density", "vessel_diam_um")
mm = m %>% dplyr::select("PUTATIVE_SPECIES", keep1) %>% 
  slice(grep("Encelia", m$PUTATIVE_SPECIES))
mm1 = mm %>% dplyr::select(-wood_density)
mm2 = mm1 %>% filter(PUTATIVE_SPECIES != "Encelia_ravenii", 
                     PUTATIVE_SPECIES != "Encelia_resinifera")
pcres = pca(mm2 %>% dplyr::select(-PUTATIVE_SPECIES), 
            scale = "vector", center = T, nPcs = 3)
morph = data.frame(pcres@scores)
rownames(morph) = mm2$PUTATIVE_SPECIES
morphdist = as.matrix(dist(morph, method = "euclidean", diag = T, upper = T))
for (i in 1:nrow(res)) {
  sp1 = res[i, "species1"]
  sp2 = res[i, "species2"]
  if (sp1 %in% labels(morphdist)[[1]] && sp2 %in% labels(morphdist)[[2]]) {
    res[i, "morph_dist"] = morphdist[sp1, sp2]
  }
}

####################
# calc clim dist
####################

d = read_csv("~/Dropbox/Encelia/analysis/spatial_analyses/encelia/all_points_thinned.csv")
climr = list.files("~/Dropbox/Encelia/spatial_rasters/wc2.0_30s_bio/", pattern=".tif", full.names = T)
climr = stack(lapply(climr, raster))

# need to make points unique
d1 = unique(d)
# for now, no artesmia or bailey
d1 = d1 %>% filter(!species %in% c("Artemisia tridentata", 
                                   "Baileya multiradiata",
                                   "Enceliopsis covillei",
                                   "Xylorhiza tortifolia"))
d1 = as.data.frame(d1) %>% filter(Longitude != -74.37556)
d1$species = gsub(" ", "_", d1$species)
d1[d1$species == "Encelia_actonii", "species"] = "Encelia_actoni"
d1[d1$species == "Encelia_frutescens", "species"] = "Encelia_frutescens_frutescens"

d1[(d1$species == "Encelia_farinosa" & d1$Latitude < 31.6), "species"] = "Encelia_farinosa_farinosa"
d1[(d1$species == "Encelia_farinosa" & d1$Latitude >= 31.6), "species"] = "Encelia_farinosa_phenicodonta"

pts = d1[, c("Longitude", "Latitude")]

clim = raster::extract(climr, pts)
climpts = cbind(d1, clim)
climpts2 = climpts[complete.cases(climpts), ]
climpca = prcomp(climpts2 %>% 
                   dplyr::select(-species, -Longitude, -Latitude), 
                 center = T, scale. = T)
climpcax = as.data.frame(climpca[['x']])
allpts1 = cbind(climpts2 %>% dplyr::select(species), climpcax)
clim2 = allpts1 %>% group_by(species) %>% summarise_all(mean) %>% ungroup()
clim3 = as.matrix(dist(clim2 %>% 
                         dplyr::select(-species, PC1),
               method = "euclidean", diag = T, upper = T))
rownames(clim3) = clim2$species
colnames(clim3) = clim2$species

for (i in 1:nrow(res)) {
  sp1 = res[i, "species1"]
  sp2 = res[i, "species2"]
  if (sp1 %in% labels(clim3)[[1]] && sp2 %in% labels(clim3)[[1]]) {
    res[i, "clim_dist"] = clim3[sp1, sp2]
  }
}

####################
# calc soil dist
####################

soilr = list.files("~/Dropbox/Encelia/spatial_rasters/NACP_MSTMIP_UNIFIED_NA_SOIL_MA_1242/", pattern=".tif", full.names = T)
soilr = stack(lapply(soilr, raster))

soil = raster::extract(soilr, pts)
soilpts = cbind(d1, soil)
soilpts2 = soilpts[complete.cases(soilpts), ]

soilpca = prcomp(soilpts2 %>% 
                   dplyr::select(-species, -Longitude, -Latitude), 
                 center = T, scale. = T)
soilpcax = as.data.frame(soilpca[['x']])
allpts1 = cbind(soilpts2 %>% dplyr::select(species), soilpcax)
soil2 = allpts1 %>% group_by(species) %>% summarise_all(mean) %>% ungroup()
soil3 = as.matrix(dist(soil2 %>% dplyr::select(-species, PC1),
                       method = "euclidean", diag = T, upper = T))
rownames(soil3) = soil2$species
colnames(soil3) = soil2$species

for (i in 1:nrow(res)) {
  sp1 = res[i, "species1"]
  sp2 = res[i, "species2"]
  
  if (sp1 %in% labels(soil3)[[1]] && sp2 %in% labels(soil3)[[1]]) {
    res[i, "soil_dist"] = soil3[sp1, sp2]
  }
}


####################
# calc geo dist
####################

ranges = list.files("~/Dropbox/Encelia/analysis/spatial_analyses/encelia_ranges/", pattern = "shp" ,full.names = T)
ranges2 = lapply(ranges, rgdal::readOGR)
centroids = lapply(ranges2, rgeos::gCentroid)
spranges = gsub(".*\\/", "", ranges)
spranges = gsub(".shp", "", spranges)
names(ranges2) = spranges
names(centroids) = spranges

for (i in 1:nrow(res)) {
  sp1 = res[i, "species1"]
  sp2 = res[i, "species2"]
  if (sp1 == "Encelia_farinosa_farinosa") {
    sp1 = "Encelia_farinosa"
  }
  if (sp2 == "Encelia_farinosa_farinosa") {
    sp2 = "Encelia_farinosa"
  }
  if (sp1 == "Encelia_farinosa_phenicodonta") {
    sp1 = "Encelia_farinosa"
  }
  if (sp2 == "Encelia_farinosa_phenicodonta") {
    sp2 = "Encelia_farinosa"
  }
  if (sp1 %in% names(centroids) && sp2 %in% names(centroids)) {
    res[i, "geo_dist"] = sp::spDistsN1(centroids[[sp1]], centroids[[sp2]])
    sp1r = as(ranges2[[sp1]], "SpatialPolygons")
    sp2r = as(ranges2[[sp2]], "SpatialPolygons")
   
    r_inter = gIntersection(sp1r, sp2r)
    if (is.null(r_inter)) {
      res[i, "per_overlap"] = 0
    } else {
      minarea = min(gArea(sp1r), gArea(sp2r))
     res[i, "per_overlap"] = gArea(r_inter) / minarea
    }
  }
}


####################
# add in adjacency data
####################

adj = read.csv("~/Dropbox/Encelia/analysis/spatial_analyses/species_adjacency.csv",
               stringsAsFactors = F)
res2 = inner_join(res, adj)


#########################
# some analysis
##########################

res2 = res2 %>% mutate(overlap = ifelse(per_overlap > 0, TRUE, FALSE))
res2 %>% filter(species1 != "Encelia_canescens", species2 != "Encelia_canescens") %>% group_by(overlap) %>% summarize_all(mean, na.rm = T)

a =
  ggplot(res2 %>% filter(complete.cases(morph_dist)), 
           aes(clim_dist, morph_dist, fill = overlap)) + 
  geom_point(shape = 21) + xlab("climatic dist.") +
  ylab("morph. dist.") + 
  scale_fill_manual(values = c("white", "black"))
b =
  ggplot(res2 %>% filter(complete.cases(morph_dist)), 
         aes(soil_dist, morph_dist, fill = overlap)) + 
  geom_point(shape = 21) + xlab("soil dist.") +
  ylab("morph. dist.") + 
  scale_fill_manual(values = c("white", "black"))
ab = plot_grid(a + theme(legend.position="none"), 
               b + theme(legend.position="none"), labels = c("A", "B"))
legend <- get_legend(
  # create some space to the left of the legend
  a + theme(legend.box.margin = margin(t = 0, b = 0, r = 5, l = 12))
)
# the width of one plot (via rel_widths).
abl = plot_grid(ab, legend, rel_widths = c(3, .4))

save_plot("~/Dropbox/Encelia/manuscripts/Encelia_Phylogeny/figures/species_comparisons.png", abl, ncol = 2, base_height = 3, base_width = 4)

M = cor(res2 %>% dplyr::select(-species1, -species2, -adjacent),
    method = "spearman", use = "pairwise.complete.obs")
corrplot::corrplot.mixed(M, tl.cex = 0.7, tl.pos = "lt")

###############################
# simulations
###############################

do_sim <- function(sps, 
                   phylodists.right.order,
                   var, trait.data) {
                   
  # get trait dist
  trait.dist = matrix(NA, nrow = length(sps), ncol = length(sps))
  for (i in 1:length(sps)) {
    for (j in 1:length(sps)) {
      d1 = res2[(res2$species1 == sps[i] & res2$species2 == sps[j]), ]
      d2 = res2[(res2$species2 == sps[i] & res2$species1 == sps[j]), ]
      d = rbind(d1, d2)
      trait.dist[i, j] = d[1, var]
    }
  }

  # observed means
  observed.temp<-split(trait.dist[upper.tri(trait.dist)]/phylodists.right.order[upper.tri(phylodists.right.order)],sym.allo[upper.tri(sym.allo)])
  observed.temp<-unlist(lapply(observed.temp, mean, na.rm = T))
  obs.diff.temp<-observed.temp[1]-observed.temp[2]
  obs.diff = obs.diff.temp
  
  # simulate difference in means
  td.temp<-treedata(multi2di(spp.tree), trait.data, sort = T)
  s.temp<-ratematrix(td.temp$phy, td.temp$data)
  sims.temp <-sim.char(td.temp$phy, s.temp, 10000) # makes an array of sims
  
  sim.means.temp <- vector("numeric", 10000)
  for (n in 1:length(sim.means.temp)){
    
    trait.temp<-sims.temp[,,n]
    trait.temp<-trait.temp[sps]
    trait.dist.matrix.temp<-as.matrix(dist(trait.temp)) 
    #make distance matrix (using appropriate distance function for each dataset)
      
    simed<-split(trait.dist.matrix.temp[upper.tri(trait.dist.matrix.temp)]/phylodists.right.order[upper.tri(phylodists.right.order)],sym.allo[upper.tri(sym.allo)])
    simed<-unlist(lapply(simed, mean, na.rm = T))
    sim.means.temp[n]<-simed[1]-simed[2]
  }
  
  return(list(sim.means.temp, obs.diff))
  }


# tree 
spp.tree = drop.tip(td_nom, c("Enceliopsis_covillei", 
                              "Xylorhiza_tortifolia"))
sps = spp.tree$tip.label

# get sym allo
sym.allo = matrix(NA, nrow = length(sps), ncol = length(sps))
for (i in 1:length(sps)) {
  for (j in 1:length(sps)) {
    d1 = res2[(res2$species1 == sps[i] & res2$species2 == sps[j]), ]
    d2 = res2[(res2$species2 == sps[i] & res2$species1 == sps[j]), ]
    d = rbind(d1, d2)
    sym.allo[i, j] = d[1, "adjacent"]
  }
}
sym.allo[sym.allo == TRUE] = 1
sym.allo[sym.allo == FALSE] = 0


phylodist.this.tree<-cophenetic(spp.tree)
#re-order phylo.dist.temp to the right order
dim.number<-length(sps)
phylodists.right.order<-matrix(NA,nrow=dim.number,ncol=dim.number)
for (x in 1:dim.number){
  for (j in 1:dim.number){
    phylodists.right.order[x,j]<-phylodist.this.tree[sps[x], sps[j]]
  }
}
rownames(phylodists.right.order)<-colnames(phylodists.right.order)<-sps


morph.trait = morph$PC1
names(morph.trait) = rownames(morph)
clim.trait = clim2$PC1
names(clim.trait) = clim2$species
soil.trait = soil2$PC1
names(soil.trait) = soil2$species

res2$morph_dist2 = res2$morph_dist / res2$phy_dist
res2$clim_dist2 = res2$clim_dist / res2$phy_dist
res2$soil_dist2 = res2$soil_dist / res2$phy_dist


morphres = do_sim(sps, phylodists.right.order, "morph_dist", morph.trait)
ecdf(morphres[[1]])(morphres[[2]])

soilres = do_sim(sps, phylodists.right.order, "soil_dist", soil.trait)
ecdf(soilres[[1]])(soilres[[2]])

climres = do_sim(sps, phylodists.right.order, "clim_dist", clim.trait)
ecdf(climres[[1]])(climres[[2]])
ecdf(climres[[1]])(-0.41)

res2$sym = ifelse(res2$adjacent == TRUE, "sympatric", "allopatric")
a = ggplot(res2, aes(sym, clim_dist2)) + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha= 0.5) +
  geom_boxplot(outlier.shape = NA, fill = NA) + 
  xlab("range overlap") +
  ylab(expression(frac("clim. dist.", "phy. dist.")))
b = ggplot(res2, aes(sym, soil_dist2)) + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha= 0.5) +
  geom_boxplot(outlier.shape = NA, fill = NA) + 
  xlab("range overlap") + 
  ylab(expression(frac("soil dist.", "phy. dist.")))
c = ggplot(res2, aes(sym, morph_dist2)) + 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha= 0.5) +
  geom_boxplot(outlier.shape = NA, fill = NA) + 
  xlab("range overlap") + 
  ylab(expression(frac("morph. dist.", "phy. dist.")))

plot_grid(a, b, c, ncol = 3)


res2 %>% group_by(sym) %>% filter(clim_dist2 < 20) %>% summarize(mean(clim_dist2))