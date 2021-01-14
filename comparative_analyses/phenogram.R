library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(cowplot)
library(raster)
library(ape)
library(phytools)
theme_set(theme_cowplot())

########################
# phylogeny data
#########################

# prep tree
t = read.tree("~/Dropbox/Encelia/analysis/phylogeny/astral/astral.miss0.6.brlen.dated.tre")
td = read.tree(text = write.tree(ladderize(td)))

td_nom = drop.tip(td, c("Encelia_californica1", "Encelia_virginensis2",
                        "Encelia_frutescens_glandulosa"))
td_nom$tip.label = gsub("\\d", "", td_nom$tip.label)
td_nom = root(td_nom, c("Enceliopsis_covillei", "Xylorhiza_tortifolia"),
              resolve.root = T)

######################
# morph data
######################

m = read_csv("~/Dropbox/Encelia/analysis/phylogeny/morphological_data/morphological_data_spread.csv")
keep1 = c("BTSA", "top_trichomes_per_mm2",
          "Kshoot_mmol/s/m2/MPa", "leaf_gray_mean",
          "leaf_Area", 
          "leaf_Round", "leaf_LMA",
          "wood_density", "vessel_diam_um")
mm = m %>% dplyr::select("PUTATIVE_SPECIES", keep1) %>% 
  filter(PUTATIVE_SPECIES != "Artemesia_tridentata",
         PUTATIVE_SPECIES != "Baiylea_multiradiata")

######################
# env data
######################

climr = list.files("~/Dropbox/Encelia/spatial_rasters/wc2.0_30s_bio/", pattern=".tif", full.names = T)
soilr = list.files("~/Dropbox/Encelia/spatial_rasters/NACP_MSTMIP_UNIFIED_NA_SOIL_MA_1242/", pattern=".tif", full.names = T)
climr = stack(lapply(climr, raster))
soilr = stack(lapply(soilr, raster))

d = read_csv("~/Dropbox/Encelia/analysis/spatial_analyses/encelia/all_points_thinned.csv")
# need to make points unique
d1 = unique(d)
# for now, no artesmia or bailey
d1 = d1 %>% filter(!species %in% c("Artemisia tridentata", 
                                   "Baileya multiradiata")) %>% 
  filter(Longitude != -74.37556)
pts = d1[, c("Longitude", "Latitude")]
pts = as.data.frame(pts)

clim = raster::extract(climr, pts)
colnames(clim) =  gsub("wc2.0_bio_30s_", "bc", colnames(clim))
colnames(clim) =  gsub("_Encelia", "", colnames(clim))
clim2 = cbind(d1, clim)

soil = raster::extract(soilr, pts)
soilpts = cbind(d1, soil)
soilpts2 = soilpts[complete.cases(soilpts), ]
colnames(soilpts2) =  gsub("Unified_NA_Soil_Map_", "", colnames(soilpts2))
soilpca = prcomp(soilpts2 %>% 
                   dplyr::select(-species, -Longitude, -Latitude), 
                 center = T, scale. = T)
soilpcax = as.data.frame(soilpca[['x']])
names(soilpcax) = paste("soil", names(soilpcax), sep = "_")
soil2 = cbind(soilpts2, soilpcax)

env = full_join(clim2, soil2)
env2 = env %>% group_by(species) %>% 
  dplyr::select(-Longitude, -Latitude) %>% 
  summarize_all(mean, na.rm = T) %>% ungroup()
env2$species = gsub(" ", "_", env2$species)
env2[env2$species == "Encelia_actonii", "species"] = "Encelia_actoni"
env2[env2$species == "Encelia_frutescens", "species"] = "Encelia_frutescens_frutescens"
env2[env2$species == "Encelia_farinosa", "species"] = "Encelia_farinosa_farinosa"

######################
# prep data
######################

me = full_join(env2, mm, by = c("species" = "PUTATIVE_SPECIES"))
me2 = me %>% rename(Kshoot_mmol_s_m2_MPa = 'Kshoot_mmol/s/m2/MPa')

traits = c(keep1, "soil_PC1", 
           "soil_PC2", "bc05", 
           "bc06","bc13", "bc14")
ntraits = c("BTSA~(no.~tips~cm^-2)", 
            "trichome~density~(mm^-2)",
            "shoot~hydr.~(mmol/s/m2/MPa)", "leaf~color",
            "leaf~area~(cm^2)", 
            "leaf~roundness", "LMA~(g~cm^2)",
            "wood~density~(g~cm^-3)", "vessel~diam.~(um)", 
            "soil~PC1", "soil~PC2", "warmest~temp.~(C)", 
            "coldest~temp.~(C)","wettest~precip.~(mm)", "driest~precip.~(mm)")
traits[which(traits == 'Kshoot_mmol/s/m2/MPa')] = 'Kshoot_mmol_s_m2_MPa'
names(ntraits) = traits

shortnames = c("CAL", "ASP", "CAN", "FARFAR", "FARPHE",
               "VIR", "ACT", "FRUFRU", "VEN", "DEN",
               "PAL", "RES", "RAV", "XYL", "COV", "FAR")
fullnames = c("Encelia_californica", "Encelia_asperifolia", 
              "Encelia_canescens", "Encelia_farinosa_farinosa", 
              "Encelia_farinosa_phenicodonta", 
              "Encelia_virginensis", "Encelia_actoni", 
              "Encelia_frutescens_frutescens", 
              "Encelia_ventorum", "Encelia_densifolia",
              "Encelia_palmeri", "Encelia_resinifera",
              "Encelia_ravenii", "Xylorhiza_tortifolia",
              "Enceliopsis_covillei", "Encelia_farinosa")
names(shortnames) = fullnames

######################
# phenogram data
######################

me_no = me2 %>% filter(species != "Enceliopsis_covillei",
                      species != "Xylorhiza_tortifolia")
clades = list(c("FARFAR", "FARPHE", "PAL", "CAN", "CAL", "ASP"),
              c("RAV", "FRUFRU", "ACT", "VIR", "RES"),
              c("VEN", "DEN"))

intraits = c("top_trichomes_per_mm2", "Kshoot_mmol_s_m2_MPa",
             "leaf_Area", "leaf_LMA", 
             "bc05", "bc06", "bc14",
             "soil_PC1")
png("~/Dropbox/Encelia/manuscripts/Encelia_Phylogeny/figures/phenograms_main.png", height = 3, width = 7,
    units = "in", res = 200)
par(mfrow = c(2, 4), mar = c(2, 2, 1, 1))
for (i in 1:length(intraits)) {
  trait = pull(me_no[ match(td_nom$tip.label, me_no$species), intraits[i]])
  names(trait) = td_nom$tip.label
  trait = trait[complete.cases(trait)]
  
  keepsp = intersect(names(trait), td_nom$tip.label)
  td_nomx = drop.tip(td_nom, setdiff(td_nom$tip.label, keepsp))
  names(trait) = shortnames[names(trait)]
  td_nomx$tip.label = shortnames[td_nomx$tip.label]
  
  painted = td_nomx
  ncols = c("#8dd3c7", "#fb8072", "#80b1d3")
  for (c in 1:length(clades)) {
    tips = intersect(clades[[c]], td_nomx$tip.label)
    if (length(tips) > 1) {
      node = findMRCA(td_nomx, tips)
      painted = paintSubTree(painted, node=node, state=ncols[c])
    } else if (length(tips) == 1) {
      painted = paintBranches(painted, 
                              which(painted$tip.label == tips), state=ncols[c])
    }
  }
  cols<-colnames(painted$mapped.edge)
  names(cols)<-cols
  
  if (i %in% c()) {
    spread = FALSE 
    } else {
    spread = TRUE
  }
  
  ytcks = pretty(trait, 3)
  par(xaxt="n",yaxt="n",mar=c(0.5, 2.5, 2.5, 0.5))
  phenogram(painted, trait[td_nomx$tip.label], 
            spread.labels = spread, ylab = "", xlab = "",
            colors = cols, fsize = 0.7, tipcol = "red")
  par(yaxt="s")
  axis(2, at = ytcks, labels = NA, las = 2, 
       tck=-0.02, lwd = 1)
  axis(2, lwd = 0, line = -0.5, at = ytcks,
       labels = ytcks, cex.axis = 0.8, las = 1, col.axis = "gray40")
  title(parse(text = ntraits[intraits[i]]), font.main = 2, cex.main = 0.8, col.main = "gray40")
}
dev.off()


outtraits = c("BTSA", "leaf_gray_mean", "leaf_Round",
           "wood_density", "vessel_diam_um", "soil_PC2", "bc13")
png("~/Dropbox/Encelia/manuscripts/Encelia_Phylogeny/figures/phenograms_SI.png", height = 4, width = 8,
    units = "in", res = 200)
par(mfrow = c(2, 4), mar = c(2, 2, 1, 1))
for (i in 1:length(outtraits)) {
  trait = pull(me_no[ match(td_nom$tip.label, me_no$species), outtraits[i]])
  names(trait) = td_nom$tip.label
  trait = trait[complete.cases(trait)]
  
  keepsp = intersect(names(trait), td_nom$tip.label)
  td_nomx = drop.tip(td_nom, setdiff(td_nom$tip.label, keepsp))
  names(trait) = shortnames[names(trait)]
  td_nomx$tip.label = shortnames[td_nomx$tip.label]
  
  painted = td_nomx
  ncols = c("#8dd3c7", "#fb8072", "#80b1d3")
  for (c in 1:length(clades)) {
    tips = intersect(clades[[c]], td_nomx$tip.label)
    if (length(tips) > 1) {
      node = findMRCA(td_nomx, tips)
      painted = paintSubTree(painted, node=node, state=ncols[c])
    } else if (length(tips) == 1) {
      painted = paintBranches(painted, 
                              which(painted$tip.label == tips), state=ncols[c])
    }
  }
  cols<-colnames(painted$mapped.edge)
  names(cols)<-cols
  
  if (i %in% c(4)) {
    spread = FALSE 
  } else {
    spread = TRUE
  }
  
  ytcks = pretty(trait, 3)
  par(xaxt="n",yaxt="n",mar=c(0.5, 2.5, 2.5, 0.5))
  phenogram(painted, trait[td_nomx$tip.label], 
            spread.labels = spread, ylab = "", xlab = "",
            colors = cols, fsize = 0.7)
  par(yaxt="s")
  axis(2, at = ytcks, labels = NA, las = 2, 
       tck=-0.02, lwd = 1)
  axis(2, lwd = 0, line = -0.5, at = ytcks,
       labels = ytcks, cex.axis = 0.8, las = 1, col.axis = "gray40")
  title(parse(text = ntraits[outtraits[i]]), font.main = 2, cex.main = 1, col.main = "gray40")
}
dev.off()


######################
# correlations data
######################

library(GGally)

keep1a = c("BTSA", "top_trichomes_per_mm2",
           "Kshoot_mmol_s_m2_MPa", "leaf_gray_mean",
          "leaf_Area", "leaf_Round", 
          "leaf_LMA", "vessel_diam_um")
xx = me2 %>% dplyr::select(keep1a) %>% ggpairs()
plots2 = vector("list", length = length(keep1a) * length(keep1a))
track = 1
for (i in 1:length(keep1a)) {
  cat(i)
  for (j in 1:length(keep1a)) {
    cat(j, "\n")
    plots2[[track]] = getPlot(xx, i = i, j = j)
    
    if (j > i) {
      X = me2 %>% dplyr::select(species, keep1a[i], keep1a[j]) %>%
        na.omit()
      keepsp = intersect(X$species, td_nom$tip.label)
      td_nomx = drop.tip(td_nom, setdiff(td_nom$tip.label, keepsp))
      X = X[match(td_nomx$tip.label, X$species), ]
      
      obj = phyl.vcv(as.matrix(log(X[, 2:3])), vcv(td_nomx), 1)
      
      corvar = cov2cor(obj$R)[1, 2]
      t.xy = corvar * sqrt((Ntip(td_nomx)-2)/(1-corvar^2))
      pval = 2*pt(abs(t.xy),df=Ntip(td_nomx)-2,lower.tail=F)
      
      label = paste0("r = ", round(corvar, 3), ";\np-val = ", 
                     round(pval, 2))
      tcol = "black"
      if (pval < 0.05) {
        tcol = "red"
      }
      
      plots2[[track]] = ggplot(me2, aes_string(keep1a[j], keep1a[i])) + 
        geom_point(size = 0, col = "white") + 
        annotation_custom(grid::textGrob(label),
                          xmin = -Inf, xmax = Inf, 
                          ymin = -Inf, ymax = Inf) 
    }
    track = track + 1
  }
}

keep1an = c("BTSA", "trichomes",
  "shoot hydr.", "leaf color",
  "leaf area", "leaf roundness", 
  "LMA", "vessel diam.")
png("~/Dropbox/Encelia/manuscripts/Encelia_Phylogeny/figures/trait_correlation.png", height = 15, width = 15, res = 200, units = "in")
ggplot2::theme_set(ggplot2::theme_bw())
# wrap_plots(plots, ncol = 9)
ggmatrix(plots2, 
          nrow = length(keep1a), 
          ncol= length(keep1a), 
          xAxisLabels = keep1an, 
          yAxisLabels = keep1an)
dev.off()