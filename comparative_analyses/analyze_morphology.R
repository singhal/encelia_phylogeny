library(ggplot2)
library(corrplot)
library(cowplot)
library(readr)
library(ape)
library(phytools)
library(tidyr)
library(dplyr)
library(geiger)
library(raster)
library(nlme)
library(pcaMethods)
library(patchwork)
theme_set(theme_cowplot())

setwd("~/Dropbox/Encelia/manuscripts/Encelia_Phylogeny/figures/")

#############################
# get morphological data
#############################

m = read_csv("~/Dropbox/Encelia/analysis/phylogeny/morphological_data/morphological_data_spread.csv")
keep = c("leaf_Area", "leaf_Round", "bottom_trichomes_per_mm2", 
         "BTSA",
         "Kshoot_mmol/s/m2/MPa", "leaf_PC1", "leaf_PC2", "leaf_gray_mean",
         "top_trichomes_per_mm2", "vessel_diam_um", "wood_density",
         "leaf_LMA")


#############################
# PCA
#############################

keep1 = c("BTSA", "top_trichomes_per_mm2",
          "Kshoot_mmol/s/m2/MPa", "leaf_gray_mean",
          "leaf_Area", 
          "leaf_Round", "leaf_LMA",
          "wood_density", "vessel_diam_um")
mm = m %>% dplyr::select("PUTATIVE_SPECIES", keep1) %>% 
  slice(grep("Encelia", m$PUTATIVE_SPECIES))
apply(mm, 1, function(x) { sum(complete.cases(x))})
mm1 = mm %>% dplyr::select(-wood_density, -vessel_diam_um)
mm2 = mm1 %>% filter(PUTATIVE_SPECIES != "Encelia_ravenii", 
                     PUTATIVE_SPECIES != "Encelia_resinifera")
pcres = pca(mm2 %>%  dplyr::select(-PUTATIVE_SPECIES), 
            scale = "vector", center = T, nPcs = 7)

mm3 = mm2 %>%  dplyr::select(-`Kshoot_mmol/s/m2/MPa`)
td_nom2 = drop.tip(td_nom, setdiff(td_nom$tip.label, mm3$PUTATIVE_SPECIES))
mm4 = mm3[match(td_nom2$tip.label, mm3$PUTATIVE_SPECIES), ] %>% 
  dplyr::select(-PUTATIVE_SPECIES)
mm4a = scale(mm4)
pcres2 = phyl.pca(td_nom2, mm4a)

################################################
# plot traits on phylogeny & calculate phylo sig
################################################

# prep tree
t = read.tree("~/Dropbox/Encelia/analysis/phylogeny/astral/astral.miss0.6.brlen.dated.tre")
t1 = read.tree("~/Dropbox/Encelia/analysis/phylogeny/astral/astral.miss0.6.tre")
t1 = minRotate(t1, setNames(1:Ntip(t), t$tip.label))
t$node.label = as.numeric(t1$node.label)
td = read.tree(text = write.tree(ladderize(t)))

# rooten21 = getMRCA(to, to$tip.label[grep("Encelia_", to$tip.label)])
# tod = chronopl(to, 0.0001, node = rootenc1, age.min = 1.36)
# tod = read.tree(text = write.tree(ladderize(tod)))
td_nom = drop.tip(td, c("Encelia_californica1", "Encelia_virginensis2",
                        "Encelia_frutescens_glandulosa"))
td_nom$tip.label = gsub("\\d", "", td_nom$tip.label)
td_nom = root(td_nom, c("Enceliopsis_covillei", "Xylorhiza_tortifolia"),
              resolve.root = T)

get_vals <- function(dd, var, name) {
  dd[ , var] = as.numeric(pull(dd[, var]))
  par(mar=c(5, 1, 1, 1), tck=-0.02)
  xvals = pretty(range(dd[, var], na.rm = T), 2)
  plot.new()
  plot.window(xlim=range(xvals, na.rm = T), ylim=c(1, nrow(dd)))
  
  for (i in 1:nrow(dd)) {
    # lines(x = c(dd[i, varmin], dd[i, varmax]), y = c(i, i))
    points(dd[i, var], i, pch=21, bg="#377eb8", cex=1)
  }
  
  axis(1, at=xvals, labels=NA, line=1)
  axis(1, at=xvals, labels=xvals, line=0.4, lwd = 0, cex.axis=1)
  mtext(name, side=1, line=3.3, cex=1)
}

mt = m[match(td_nom$tip.label, m$PUTATIVE_SPECIES), ]

# pdf("Encelia_morphology.pdf", width = 10, height = 5)
# par(mfrow=c(2, 5))
# par(mar=c(5,0,1,0))
# plot(td_nom, show.tip.label = F, edge.width=1.5, cex=1.2)
# # PC1 is mainly capturing area
# get_vals(mt, "leaf_PC1", "leaf PC1")
# get_vals(mt, "leaf_PC2", "leaf PC2")
# get_vals(mt, "leaf_LMA", "LMA")
# get_vals(mt, "leaf_gray_mean", "leaf color")
# plot(td_nom, show.tip.label = F, edge.width=1.5, cex=1.2)
# get_vals(mt, "BTSA", "shoot ramification")
# get_vals(mt, "Kshoot_mmol/s/m2/MPa", "shoot hydr.")
# get_vals(mt, "vessel_diam_um", "vessel diam.")
# # top & bottom trichomes are so similar
# get_vals(mt, "top_trichomes_per_mm2", "trichome dens.")
# dev.off()

physig = vector('list', length(keep1))
names(physig) = keep1
for (i in 1:length(keep1)) {
  mtmp = pull(m[, keep1[i]])
  names(mtmp) = m$PUTATIVE_SPECIES
  mtmp = mtmp[complete.cases(mtmp)]
  stmp = intersect(names(mtmp), td_nom$tip.label)
  td_nom1 = drop.tip(td_nom, setdiff(td_nom$tip.label, stmp))
  physig[[i]] = phylosig(td_nom1, mtmp[td_nom1$tip.label],
                         method = "lambda", test = T)  
  physig[[i]][['n']] = Ntip(td_nom1)
}
physig2 = data.frame(traits = keep1,
                     lambda = round(unlist(lapply(physig, function(x) {x$lambda})), 2),
                     pval = round(unlist(lapply(physig, function(x) {x$P})), 2),
                     ntip = unlist(lapply(physig, function(x) {x$n})))
write.csv(physig2, "phylogenetic_signal_lambda.csv", row.names = F)

###########################
# rates
###########################

rates = vector('list', length(keep1))
for (i in 1:length(keep1)) {
  trait = pull(mt[ , keep1[i]])
  names(trait) = mt$PUTATIVE_SPECIES
  trait = trait[complete.cases(trait)]
  keepsp = intersect(names(trait), td_nom$tip.label)
  td_nomx = drop.tip(td_nom, setdiff(td_nom$tip.label, keepsp))

  raw = fitContinuous(td_nomx, log(trait))
  scale = fitContinuous(td_nomx, scale(trait))
  pic = mean(pic(log(trait), td_nomx) ^ 2)
  rates[[i]] = list(raw, scale, pic)
}

# evo_rates = cbind(
#   raw_beta = unlist(lapply(rates, function(x) {x[[1]]$opt$sigsq})),
#   scaled_beta = unlist(lapply(rates, function(x) {x[[2]]$opt$sigsq})),
#  pic_beta = unlist(lapply(rates, function(x) {x[[3]]})))
# evo_rates = as.data.frame(evo_rates)
# rownames(evo_rates) = keep1

evo_rates = physig2
beta = round(unlist(lapply(rates, function(x) {x[[1]]$opt$sigsq})), 2)
names(beta) = keep1
evo_rates$beta = beta[rownames(physig2)]
write.csv(evo_rates, "evolutionary_rates.csv", row.names = F)

###########################
# environmental analysis
###########################

d = read_csv("~/Dropbox/Encelia/analysis/spatial_analyses/encelia/all_points_thinned.csv")

climr = list.files("~/Dropbox/Encelia/spatial_rasters/wc2.0_30s_bio/", pattern=".tif", full.names = T)
soilr = list.files("~/Dropbox/Encelia/spatial_rasters/NACP_MSTMIP_UNIFIED_NA_SOIL_MA_1242/", pattern=".tif", full.names = T)
climr = stack(lapply(climr, raster))
soilr = stack(lapply(soilr, raster))

# need to make points unique
d1 = unique(d)
# for now, no artesmia or bailey
d1 = d1 %>% filter(!species %in% c("Artemisia tridentata", 
                                  "Baileya multiradiata"))
d1[(d1$species == "Encelia farinosa" & d1$Latitude < 31.6), "species"] = "Encelia farinosa farinosa"
d1[(d1$species == "Encelia farinosa" & d1$Latitude >= 31.6), "species"] = "Encelia farinosa phenicodonta"

pts = d1[, c("Longitude", "Latitude")]
pts = as.data.frame(pts)

clim = raster::extract(climr, pts)
soil = raster::extract(soilr, pts)

climpts = cbind(d1, clim)
climpts2 = climpts[complete.cases(climpts), ]
climpts2 = climpts2 %>% filter(Longitude != -74.37556)
colnames(climpts2) =  gsub("wc2.0_bio_30s_", "bc", colnames(climpts2))
colnames(climpts2) =  gsub("_Encelia", "", colnames(climpts2))
climpca = prcomp(climpts2 %>% 
                   dplyr::select(-species, -Longitude, -Latitude), 
                 center = T, scale. = T)
climpcax = as.data.frame(climpca[['x']])
names(climpcax) = paste("clim", names(climpcax), sep = "_")

soilpts = cbind(d1, soil)
soilpts2 = soilpts[complete.cases(soilpts), ]
colnames(soilpts2) =  gsub("Unified_NA_Soil_Map_", "", colnames(soilpts2))
soilpca = prcomp(soilpts2 %>% 
                   dplyr::select(-species, -Longitude, -Latitude), 
                 center = T, scale. = T)
soilpcax = as.data.frame(soilpca[['x']])
names(soilpcax) = paste("soil", names(soilpcax), sep = "_")

allpts1 = cbind(climpts2, climpcax)
allpts2 = cbind(soilpts2, soilpcax)
allpts3 = full_join(allpts1, allpts2)

allpts3$species = gsub(" ", "_", allpts3$species)
allpts3[allpts3$species == "Encelia_actonii", "species"] = "Encelia_actoni"
# all these spatial data appear to be frutescens frutescens
allpts3[allpts3$species == "Encelia_frutescens", "species"] = "Encelia_frutescens_frutescens"
span <- function(col) {
  return(max(col, na.rm = T) - min(col, na.rm = T))
}
median_na <- function(col) {
  return(median(col, na.rm = T))
}
mean_na <- function(col) {
  return(mean(col, na.rm = T))
}
max_na <- function(col) {
  return(max(col, na.rm = T))
}
min_na <- function(col) {
  return(min(col, na.rm = T))
}
allpts4 = allpts3 %>% group_by(species) %>% 
  dplyr::select(-Longitude, -Latitude) %>% 
  summarise_all(list(median = median_na, range = span,
                     min = min_na, max = max_na,
                     mean = mean_na)) %>% ungroup()

### climate bivariate plot ###
groups = list(c("Encelia_ravenii", "Encelia_frutescens_frutescens",
                "Encelia_actoni", "Encelia_virginensis",
                "Encelia_resinifera"),
              c("Encelia_farinosa", "Encelia_palmeri",
                "Encelia_canescens", "Encelia_californica",
                "Encelia_asperifolia"),
              c("Encelia_ventorum", "Encelia_densifolia"))
allpts3no = allpts3 %>% filter(species != "Xylorhiza_tortifolia",
                               species != "Enceliopsis_covillei")
biplots = vector("list", length(groups))
for (i in 1:length(groups)) {
  biplots[[i]] = ggplot(allpts3no, aes(clim_PC1, clim_PC2)) + 
    geom_point(col = "gray80") +
    geom_point(data = allpts3no %>% filter(species %in% groups[[i]]),
               aes(clim_PC1, clim_PC2), color = "black", size = 2) +
    geom_point(data = allpts3no %>% filter(species %in% groups[[i]]),
               aes(clim_PC1, clim_PC2, col = species)) +
    scale_color_manual(values = c(brewer.pal(length(groups[[i]]), "Set3")),
                       labels = gsub("_", " ", sort(groups[[i]]))) +
    xlab("clim. PC1") + ylab("clim. PC2") +
    theme(legend.text = element_text(face = "italic"))
}
abc = (biplots[[1]] / biplots [[2]] / biplots[[3]] )
save_plot("bivariate_plots.png", abc, nrow = 3, ncol = 1)

m = read_csv("~/Dropbox/Encelia/analysis/phylogeny/morphological_data/morphological_data.csv")
m[m$PUTATIVE_SPECIES == "Artemesia_tridentata", "PUTATIVE_SPECIES"] = "Artemisia_tridentata"
m[m$PUTATIVE_SPECIES == "Baiylea_multiradiata", "PUTATIVE_SPECIES"] = "Baileya_multiradiata"
m[grep("farinosa", m$PUTATIVE_SPECIES), "PUTATIVE_SPECIES"] = "Encelia_farinosa"
m2 = m %>% dplyr::select(-PLANT_ID) %>% 
  group_by(PUTATIVE_SPECIES, measurement) %>% 
  summarise(mean_estimate = mean(as.numeric(estimate))) %>% ungroup()
m3 = m2 %>% spread(measurement, mean_estimate)

m4 = full_join(m3, allpts4, by = c("PUTATIVE_SPECIES" = "species"))
m5 = full_join(m2, allpts4, by = c("PUTATIVE_SPECIES" = "species"))

###########################
# env mapping
###########################

envkeep = c("bc01_median", "bc12_median",
            "bc05_median", "bc06_median",
            "bc13_median", "bc14_median")
m4a = m4 %>% filter(!PUTATIVE_SPECIES %in% c("Xylorhiza_tortifolia",
                                             "Enceliopsis_covillei",
                                             "Artemisia_tridentata",
                                             "Baileya_multiradiata"))
td_nom2 = drop.tip(td_nom, "Encelia_farinosa_phenicodonta")
td_nom2$tip.label[ which(td_nom2$tip.label == "Encelia_farinosa_farinosa") ] = "Encelia_farinosa"

pdf("not_used/phenophylogenetic_env.pdf", width = 6, height = 8)
par(mar = c(3, 5, 0, 0), mfrow = c(3, 2))
envsig = vector("list", length(envkeep))
for (i in 1:length(envkeep)) {
  trait = pull(m4a[ , envkeep[i]])
  names(trait) = m4a$PUTATIVE_SPECIES
  trait = trait[complete.cases(trait)]
  keepsp = intersect(names(trait), td_nom2$tip.label)
  td_nomx = drop.tip(td_nom2, setdiff(td_nom2$tip.label, keepsp))
  
  names(trait) = shortnames[names(trait)]
  td_nomx$tip.label = shortnames[td_nomx$tip.label]
  
  envsig[[i]] = phylosig(td_nomx, trait,
           method = "lambda", test = T)
  phenogram(td_nomx, trait[td_nomx$tip.label], 
              spread.labels = T, ylab = envkeep[i])
}
dev.off()


# trait-o-gram showing correlations
# between traits
mt1 = m4a[, c("PUTATIVE_SPECIES", envkeep)]
mt1 = mt1[complete.cases(mt1), ]
td_nom3 = drop.tip(td_nom2, setdiff(td_nom2$tip.label, mt1$PUTATIVE_SPECIES))
mt1 = mt1[match(td_nom3$tip.label, mt1$PUTATIVE_SPECIES), ]
td_nom3$tip.label = shortnames[td_nom3$tip.label]
mt1 = data.frame(mt1)
rownames(mt1) = shortnames[mt1$PUTATIVE_SPECIES]
mt1$PUTATIVE_SPECIES = NULL
mt1 = as.matrix(mt1)
colnames(mt1) = envkeep
pdf("not_used/trees_trait_correlation_env.pdf", width = 10, height = 10)
obj<-fancyTree(td_nom3, type="scattergram", 
               X=mt1)
dev.off()


get_vals2 <- function(dd, var1, var2, var3, name) {
  
  par(mar=c(5, 1, 1, 1), tck=-0.02)
  xvals = pretty(range(dd[, c(var1, var2, var3)], na.rm = T), 5)
  plot.new()
  plot.window(xlim=range(xvals, na.rm = T), ylim=c(1, nrow(dd)))
  
  for (i in 1:nrow(dd)) {
    lines(x = c(dd[i, var1], dd[i, var2]), y = c(i, i), 
          lwd = 3, col = alpha("#377eb8", 0.4))
    points(dd[i, var3], i, pch=16, col="#377eb8")
  }
  
  axis(1, at=xvals, labels=NA, line=1)
  axis(1, at=xvals, labels=xvals, line=0.4, lwd = 0, cex.axis=1)
  mtext(name, side=1, line=3.3, cex=1)
}

m4c = m4[match(td_nom2$tip.label, m4$PUTATIVE_SPECIES), ]
m4c = m4c %>% mutate(month_avg_precip = bc12_median / 12)
pdf("Encelia_environment.pdf", width = 7, height = 3)
par(mfrow=c(1, 3))
par(mar=c(5,0,1,0))
plot(td_nom2, show.tip.label = F, edge.width=1.5, cex=1.2)
# PC1 is mainly capturing area
get_vals2(m4c, "bc05_median", "bc06_median", "bc01_median", "temperature")
get_vals2(m4c, "bc13_median", "bc14_median", "month_avg_precip", "precipitation")
dev.off()

pdf("Encelia_niche_breadth.pdf", width = 7, height = 3)
par(mfrow=c(1, 3))
par(mar=c(5,0,1,0))
plot(td_nom2, show.tip.label = F, edge.width=1.5, cex=1.2)
# PC1 is mainly capturing area
get_vals2(m4c, "bc01_min", "bc01_max", "bc01_median", "temperature")
get_vals2(m4c, "bc12_min", "bc12_max", "bc12_median", "precipitation")
dev.off()


###########################
# targeted trait correlations
###########################

m4o = m4 %>% slice(grep("Encelia_", PUTATIVE_SPECIES))

# bad science
keep1a = c("BTSA", "top_trichomes_per_mm2",
          "Kshoot_mmol_s_m2_MPa", "leaf_gray_mean",
          "leaf_Area", 
          "leaf_Round", "leaf_LMA",
          "wood_density", "vessel_diam_um")
m4o = m4o %>% rename(Kshoot_mmol_s_m2_MPa = 'Kshoot_mmol/s/m2/MPa')
extremeenv = c("bc05_mean", "bc06_mean",
                "bc13_mean", "bc14_mean",
               "soil_PC1_mean", "soil_PC2_mean")
corres = data.frame(env = rep(NA, 54),
                    trait = rep(NA, 54),
                    cor = rep(NA, 54),
                    pval = rep(NA, 54),
                    stringsAsFactors = F)
track = 1
for (i in 1:length(extremeenv)) {
  for (j in 1:length(keep1a)) {
    X = m4o %>% dplyr::select(PUTATIVE_SPECIES, extremeenv[i], keep1a[j]) %>%
      na.omit()
    keepsp = intersect(X$PUTATIVE_SPECIES, td_nom2$tip.label)
    td_nomx = drop.tip(td_nom2, setdiff(td_nom2$tip.label, keepsp))
    
    X[, 3] = log(X[, 3])
    X = X[match(td_nomx$tip.label, X$PUTATIVE_SPECIES), ]
    obj = phyl.vcv(as.matrix(X[, 2:3]), vcv(td_nomx), 1)
    
    corvar = cov2cor(obj$R)[1, 2]
    t.xy = corvar * sqrt((Ntip(td_nomx)-2)/(1-corvar^2))
    pval = 2*pt(abs(t.xy),df=Ntip(td_nomx)-2,lower.tail=F)
    
    corres[track, "env"] = extremeenv[i]
    corres[track, "trait"] = keep1a[j]
    corres[track, "cor"] = corvar
    corres[track, "pval"] = pval
    track = track + 1
  }
}

corres %>% filter(pval < 0.1) %>% arrange(trait)

heat = c("bc05_mean", "bc06_mean")
rain = c("bc13_mean", "bc14_mean")

# trichomes & dryness (fewer trichomes, more wet - neg)
corres %>% filter(env %in% rain, trait == "top_trichomes_per_mm2")
# in opposing directions, not sig

# trichomes & heat (more trichomes, more heat - pos)
corres %>% filter(env %in% heat, trait == "top_trichomes_per_mm2")
# in predicted direction, not sig

# leaf area & dryness (smaller leaves, more dry - pos)
corres %>% filter(env %in% rain, trait == "leaf_Area")
# in opposing direction, not sig

# leaf area & heat (more heat, smaller leaves - neg)
corres %>% filter(env %in% heat, trait == "leaf_Area")
# in opposing directions, not sig

# lma & dryness (more dry, more sla, less lma - pos)
corres %>% filter(env %in% rain, trait == "leaf_LMA")
#in opposing direction, not sig

# lma & heat (more heat, less sla, more LMA - pos)
corres %>% filter(env %in% heat, trait == "leaf_LMA")
# in predicted direction, sig

corres2 = filter(corres, grepl("bc", env))
corres2 = corres2 %>% arrange(trait) %>% dplyr::select(trait, env, cor, pval)
write.csv(corres2, "trait_env_correlations.csv", row.names = F, quote = F)

########################
# plot correlation data
#######################

traitlocs = seq(0, length.out= length(keep1a), by = 2)
envlocs = seq(0, length.out = length(extremeenv), by = 3)
ntraits = c("BTSA", "trichomes",
            "shoot hydr.",
            "leaf color", "leaf area",
            "leaf round", "LMA",
            "wood dens.", "vessel diam.")
nenv = c("warmest temp.", "coldest temp.",
         "wettest precip.", "driest precip.",
         "soil PC1", "soil PC2")

png("trait_environment_correlations.png", height = 8, 
    width = 5, units = "in", res = 200)
par(xpd = T, mar = c(0, 3, 0, 3))
plot.new()
plot.window(ylim = range(traitlocs), xlim=c(1, 10))
text(y = traitlocs, x = rep(1, length(traitlocs)),
     labels = ntraits)
text(y = envlocs, x = rep(10, length(envlocs)),
     labels = nenv)
for (i in 1:nrow(corres)) {
  if (corres[i, "pval"] < 0.05) {
    lcol = "red"
  } else {
    lcol = "gray60"
  }
  if (corres[i, "cor"] < 0) {
    llty = "dotted"
  } else {
    llty = "solid"
  }
  lcex = abs(corres[i, "cor"]) * 4
  ystart = traitlocs[which(keep1a == corres[i, "trait"])] 
  yend = envlocs[which(extremeenv == corres[i, "env"])] 
  lines(x = c(3, 8), y = c(ystart, yend), lwd = lcex, col = lcol, 
        lty = llty)
}
dev.off()

###########################
# extreme value - environmental analysis
###########################

d = read_csv("~/Dropbox/Encelia/analysis/spatial_analyses/encelia/all_points_thinned.csv")

climr = list.files("~/Dropbox/Encelia/spatial_rasters/wc2.0_30s_bio/", pattern=".tif", full.names = T)
# bio5: Max Temp  warmest month
# bio6: Min. Temp. coldest month
# bio14: Preciptn. driest month
# bio13: Preciptn. wettest month
climr_sub = climr[c(5, 6, 13, 14)]
climr = stack(lapply(climr_sub, raster))

# need to make points unique
d1 = unique(d)
# for now, no outgroups
d1 = d1 %>% filter(!species %in% c("Artemisia tridentata", 
                                   "Baileya multiradiata",
                                   "Xylorhiza tortifolia",
                                   "Enceliopsis covillei"))
pts = d1[, c("Longitude", "Latitude")]
pts = as.data.frame(pts)

clim = raster::extract(climr, pts)

climpts = cbind(d1, clim)
climpts2 = climpts[complete.cases(climpts), ]
colnames(climpts2) =  gsub("wc2.0_bio_30s_", "bc", colnames(climpts2))
colnames(climpts2) =  gsub("_Encelia", "", colnames(climpts2))

# remove outlier
climpts2 = climpts2 %>% filter(Longitude != -74.37556)

climpca = prcomp(climpts2 %>% 
                   dplyr::select(-species, -Longitude, -Latitude), 
                 center = T, scale. = T)
# pca1 is 47% and is more temp (but is a bit of both)
# pca2 is 80% and is more precip (but is a bit of both)
climpcax = as.data.frame(climpca[['x']])
names(climpcax) = paste("clim", names(climpcax), sep = "_")

allpts3 = cbind(climpts2, climpcax)

allpts3$species = gsub(" ", "_", allpts3$species)
allpts3[allpts3$species == "Encelia_actonii", "species"] = "Encelia_actoni"
# all these spatial data appear to be frutescens frutescens
allpts3[allpts3$species == "Encelia_frutescens", "species"] = "Encelia_frutescens_frutescens"

allpts4 = allpts3 %>% group_by(species) %>% 
  dplyr::select(-Longitude, -Latitude) %>% 
  summarise_all(median, na.rm=T) %>% ungroup()

### climate bivariate plot ###
clades = c("frutescens", "californica", "densifolia & ventorum")
groups = list(c("Encelia_ravenii", "Encelia_frutescens_frutescens",
                "Encelia_actoni", "Encelia_virginensis",
                "Encelia_resinifera"),
              c("Encelia_farinosa", "Encelia_palmeri",
                "Encelia_canescens", "Encelia_californica",
                "Encelia_asperifolia"),
              c("Encelia_ventorum", "Encelia_densifolia"))
biplots = vector("list", length(groups) * 2)
for (i in 1:length(groups)) {
  biplots[[i]] = ggplot(allpts3, aes(bc05, bc06)) + 
    geom_point(col = "gray80") +
    geom_point(data = allpts3 %>% filter(species %in% groups[[i]]),
               aes(bc05, bc06), color = "black", size = 2) +
    geom_point(data = allpts3 %>% filter(species %in% groups[[i]]),
               aes(bc05, bc06, col = species)) +
    scale_color_manual(values = c(brewer.pal(length(groups[[i]]), "Set3")),
                       labels = gsub("_", " ", sort(groups[[i]]))) +
    xlab("warmest temp.") + ylab("coldest temp.") +
    theme(legend.position = "none")
}
for (i in 1:length(groups)) {
  biplots[[i + 3]] = ggplot(allpts3, aes(bc13, bc14)) + 
    geom_point(col = "gray80") +
    geom_point(data = allpts3 %>% filter(species %in% groups[[i]]),
               aes(bc13, bc14), color = "black", size = 2) +
    geom_point(data = allpts3 %>% filter(species %in% groups[[i]]),
               aes(bc13, bc14, col = species)) +
    scale_color_manual(values = c(brewer.pal(length(groups[[i]]), "Set3")),
                       labels = gsub("_", " ", sort(groups[[i]]))) +
    xlab("wettest precip.") + ylab("driest precip.") +
    theme(legend.text = element_text(face = "italic")) + 
    labs(color = clades[i])
}
abc = ( (biplots[[1]] / biplots [[2]] / biplots[[3]]) |
        (biplots[[4]] / biplots [[5]] / biplots[[6]]) )
save_plot("bivariate_plots_extreme.png", abc, nrow = 3, ncol = 2)

### climate & traits analysis ###

m = read_csv("~/Dropbox/Encelia/analysis/phylogeny/morphological_data/morphological_data.csv")
m[m$PUTATIVE_SPECIES == "Artemesia_tridentata", "PUTATIVE_SPECIES"] = "Artemisia_tridentata"
m[m$PUTATIVE_SPECIES == "Baiylea_multiradiata", "PUTATIVE_SPECIES"] = "Baileya_multiradiata"
m[grep("farinosa", m$PUTATIVE_SPECIES), "PUTATIVE_SPECIES"] = "Encelia_farinosa"
m2 = m %>% dplyr::select(-PLANT_ID) %>% 
  group_by(PUTATIVE_SPECIES, measurement) %>% 
  summarise(mean_estimate = mean(as.numeric(estimate))) %>% ungroup()
m3 = m2 %>% spread(measurement, mean_estimate)

m4 = full_join(m3, allpts4, by = c("PUTATIVE_SPECIES" = "species"))
m5 = full_join(m2, allpts4, by = c("PUTATIVE_SPECIES" = "species"))

m5k = m5 %>% filter(measurement %in% keep1)
f1 = ggplot(m5k, aes(clim_PC1, mean_estimate)) + geom_point() + 
  facet_wrap(~ measurement, scales = "free") + 
  xlab("climate, PC1") + ylab("estimate")
save_plot("not_used/extreme_climate_PC1_traits.pdf", f1, 
          nrow = 3, ncol = 3, base_height = 2, base_width = 3)
f2 = ggplot(m5k, aes(clim_PC2, mean_estimate)) + geom_point() + 
  facet_wrap(~ measurement, scales = "free") + 
  xlab("climate, PC2") + ylab("estimate")
save_plot("not_used/extreme_climate_PC2_traits.pdf", f2, 
          nrow = 3, ncol = 3, base_height = 2, base_width = 3)

climvars = c("bc05", "bc06", "bc13", "bc14")
climnames= c("max. temp. of warmest month",
             "min. temp. of coldest month",
             "precip. of wettest month",
             "precip. of driest month")
for (i in 1:length(climvars)) {
  tmp = ggplot(m5k, aes_string(climvars[i], "mean_estimate")) + geom_point() + 
    facet_wrap(~ measurement, scales = "free") + 
    xlab(climnames[i]) + ylab("estimate")
  save_plot(paste0(climvars[i], "_traits.pdf"), tmp, 
            nrow = 3, ncol = 3, base_height = 2, base_width = 3)
}

##########################
# disparity
#########################

m = read_csv("~/Dropbox/Encelia/analysis/phylogeny/morphological_data/morphological_data_spread.csv")
me = m %>% slice(grep("Encelia", PUTATIVE_SPECIES))

dttr = vector("list", length(keep1))
for (i in 1:length(keep1)) {
  trait = pull(me %>% dplyr::select(keep1[i]))
  names(trait) = me$PUTATIVE_SPECIES
  trait = trait[complete.cases(trait)]
  dttr[[i]] = dtt(td_nom, trait, index = "avg.sq", 
             nsim = 100, calculateMDIp = T, plot = F)
}

envkeep = c("soil_PC1_mean",  "soil_PC2_mean", 
            "bc05_mean", "bc06_mean","bc13_mean", "bc14_mean")
dttr2 = vector("list", length(envkeep))
for (i in 1:length(envkeep)) {
  trait = pull(m4o %>% dplyr::select(envkeep[i]))
  names(trait) = m4o$PUTATIVE_SPECIES
  trait = trait[complete.cases(trait)]
  dttr2[[i]] = dtt(td_nom2, trait, index = "avg.sq", 
                  nsim = 100, calculateMDIp = T, plot = F)
}

dttres = data.frame(trait = c(keep1, envkeep),
                    MDI = c(unlist(lapply(dttr, function(x) {x$MDI})),
                            unlist(lapply(dttr2, function(x) {x$MDI}))),
                    MDIpval = c(unlist(lapply(dttr, function(x) {x$MDIpVal})),
                            unlist(lapply(dttr2, function(x) {x$MDIpVal})))
                    )


alldtt = c(dttr, dttr2)
ntraits = c("BTSA", "trichomes",
            "shoot hydr.", "leaf color",
            "leaf area", 
            "leaf roundness", "LMA",
            "wood density", "vessel diam.", 
            "soil PC1", "soil PC2", "warmest temp.", 
            "coldest temp.","wettest precip.", "driest precip.")
dttplts = vector("list", length(alldtt))
for (i in 1:length(alldtt)) {
  xx = alldtt[[i]]
  dttxx = data.frame(mean = apply(xx$sim, 1, mean),
                     actual = xx$dtt,
                     times = xx$times)
  
  dttp = data.frame(y = c(apply(xx$sim, 1, quantile, 0.025), 
                          apply(xx$sim, 1, quantile, 0.975)),
                    x = c(xx$times, rev(xx$times)))
  
  dttplts[[i]] = ggplot(dttxx) + geom_polygon(data = dttp, aes(x = x, y = y), 
                               fill = "gray80", alpha = 0.6) + 
    geom_line(aes(times, actual)) + 
    geom_line(aes(times, mean), linetype = "dashed") +
    xlab("time") + ylab(paste("disparity in", ntraits[i]))
}

png("disparity_through_time.png", height = 8, width = 14,
    res = 200, units = "in")
wrap_plots(dttplts, ncol = 5)
dev.off()


###########################
# CCA
###########################

envca = allpts4 %>% dplyr::select(species, bc05_mean, bc06_mean, 
                   bc13_mean, bc14_mean) %>%
                   filter(! species %in% c("Enceliopsis_covillei",
                                           "Xylorhiza_tortifolia",
                                           "Encelia_resinifera",
                                           "Encelia_ravenii"))


m = read_csv("~/Dropbox/Encelia/analysis/phylogeny/morphological_data/morphological_data.csv")
m[m$PUTATIVE_SPECIES == "Artemesia_tridentata", "PUTATIVE_SPECIES"] = "Artemisia_tridentata"
m[m$PUTATIVE_SPECIES == "Baiylea_multiradiata", "PUTATIVE_SPECIES"] = "Baileya_multiradiata"
m2 = m %>% dplyr::select(-PLANT_ID) %>% 
  group_by(PUTATIVE_SPECIES, measurement) %>% 
  summarise(mean_estimate = mean(as.numeric(estimate))) %>% ungroup()
m3 = m2 %>% spread(measurement, mean_estimate)

morphca = m3 %>% dplyr::select(PUTATIVE_SPECIES, leaf_Area,
                               leaf_Round, BTSA, leaf_gray_mean, 
                               top_trichomes_per_mm2, leaf_LMA) %>%
  filter(! PUTATIVE_SPECIES %in% c("Enceliopsis_covillei",
                          "Xylorhiza_tortifolia", "Baileya_multiradiata",
                          "Artemisia_tridentata", "Encelia_resinifera",
                          "Encelia_ravenii"))

td_nom2 = drop.tip(td_nom, c("Xylorhiza_tortifolia", "Enceliopsis_covillei",
                             "Encelia_resinifera",
                             "Encelia_ravenii"))

envca2 = data.frame(envca, stringsAsFactors = F)
rownames(envca2) = envca2$species
envca2$species = NULL
X = scale(as.matrix(envca2[td_nom2$tip.label, ]))

morphca2 = data.frame(morphca, stringsAsFactors = F)
rownames(morphca2) = morphca2$PUTATIVE_SPECIES
morphca2$PUTATIVE_SPECIES = NULL
Y = scale(as.matrix(morphca2[td_nom2$tip.label, ]))

phyl.cca(td_nom2, X, Y, fixed=F)
