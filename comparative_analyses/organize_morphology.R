library(dplyr)
library(readr)
library(readxl)
library(tidyr)

mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}

###########################################
# organize the morphology & physiology data
###########################################

shortsps = c("actoni", "asperifolia", "californica", "canescens",
             "densifolia", "enceliopsis", "far_far", "far_phen",
             "frutescens", "pal_ven", "ravenii", "resinifera",
             "ven_pal", "ventorum", "virginensis", "xylorhiza",
             "palmeri")
sps = c("Encelia_actoni", "Encelia_asperifolia", "Encelia_californica", 
        "Encelia_canescens", "Encelia_densifolia", "Enceliopsis_covillei", 
        "Encelia_farinosa_farinosa", "Encelia_farinosa_phenicodonta",
        "Encelia_frutescens", "palmeri_ventorum_hybrid", "Encelia_ravenii", 
        "Encelia_resinifera", "ventorum_palmeri_hybrid", "Encelia_ventorum", 
        "Encelia_virginensis", "Xylorhiza_tortifolia", "Encelia_palmeri")
names(sps) = shortsps

d = read_csv("~/Dropbox/Encelia/analysis/phylogeny/morphological_data/Encelia_Samples - GENERAL.csv")
cg = read_csv("~/Dropbox/Encelia/analysis/phylogeny/morphological_data/Encelia_Samples - COMMON_GARDEN.csv")

# trichome density
t = read_xlsx("~/Dropbox/Encelia/analysis/physiology/trichomes.xlsx")
t1 = t %>% group_by(taxon, scan, surface) %>% summarise_all(mean) %>% ungroup()
t2 = t1 %>% dplyr::select(taxon, surface, trichomes_per_mm2)

t2$PUTATIVE_SPECIES = sps[t2$taxon]
t3 = t2 %>% tidyr::gather("measurement", "estimate", trichomes_per_mm2)
t3$measurement = paste(t3$surface, t3$measurement, sep = "_")
t4 = t3 %>% dplyr::select(-surface, -taxon) %>% 
  tibble::add_column(PLANT_ID = NA, .before = 1)

# trichome type
tt = unique(t %>% dplyr::select(taxon, notes))
tt$PUTATIVE_SPECIES = sps[tt$taxon]
tt$measurement = "trichome_type"
tt$estimate = tt$notes
tt = tt %>% dplyr::select(PUTATIVE_SPECIES, measurement, estimate)  %>% 
  tibble::add_column(PLANT_ID = NA, .before = 1)

# shoot ramification
r = read_csv("~/Dropbox/Encelia/analysis/phylogeny/morphological_data/Encelia_Samples - SHOOT_RAMIFICATION.csv")
r = r[grep("ENCL0", r$PLANT_ID, invert = T), ]
r$avg_radius = (r$DIAMETER_CM1 + r$DIAMETER_CM2) / 2 / 2
r$area = 3.14159 * r$avg_radius ^ 2
r$BTSA = r$TIPS / r$area
r1 = left_join(r, d %>% dplyr::select(PLANT_ID, PUTATIVE_SPECIES))
r2 = r1 %>% dplyr::select(PLANT_ID, BTSA, PUTATIVE_SPECIES)
rr = r2 %>% tidyr::gather("measurement", "estimate", BTSA)

# leaf shape & size
l = read_csv("~/Dropbox/Encelia/analysis/phylogeny/morphological_data/LeafImages_ENCILIA_folder  - Sheet1.csv")
# get rid of experimental backcrosses
l1 = l[grep("^ENCL", l$SampleID, invert = T), ]
# get rid of note-taking rows
l2 = l1[complete.cases(l1$SampleID), ]
# get rid of incomplete data
l3 = l2[complete.cases(l2$Mean), ]

# need to get rid of individuals not measured for all three color channels!!
cts = l3 %>% group_by(SampleID, LeafID) %>% 
  summarise(cts = length(unique(ColorChannel)))
l4 = left_join(l3, cts) %>% filter(cts >= 3)

# color
vals = c("Area", "Mean", "Mode", "Perimeter", 
         "Major", "Minor", "Angle", "Circ.", "AR", 
         "Round", "Solidty")
cc = list.files("~/Dropbox/Encelia/analysis/physiology/leaf_color/", 
                full.names = T, pattern = "csv")
cc1 = lapply(cc, read_csv)
cc2 = do.call(rbind, cc1)
# get sample name
ln = read_csv("~/Dropbox/Encelia/leaf_images/2017-02-24_mapping.csv")
ln$file = gsub("2017-02-24 Leaf Photos/", "", ln$file)
cc2$SampleID = pull(ln[match(gsub("JPG-corrected.tif", "JPG", cc2$Document), ln$file), "match"])
cc3 = cc2[!is.na(cc2$SampleID), ]
cc4 = cc3 %>% group_by(SampleID) %>% 
  summarise(gray_mean = mean(`Gray Value (Mean)`)) %>% ungroup()

# sum across the four color channels
l4sum = l4 %>% dplyr::select(-FileName, -ColorChannel, -Number, -X1) %>% 
  group_by(SampleID, LeafID) %>% 
  summarize_all(sum) %>% ungroup()
l4sum2 = left_join(l4sum %>% dplyr::select(-cts), 
                    cc4)
# remove mean & mode because color based 
# and colors not accurate because we didn't correct brightness
# also remove colors bc that's a diff anlaysis
pcr = prcomp(l4sum2 %>% 
               dplyr::select(-SampleID, -LeafID, -Mean, -Mode, -gray_mean),
             scale. = T, center = T)
# look at PCA loadings XXX
l5 = cbind(l4sum2, pcr[["x"]][, 1:5])

# add in leaf mass
lm = read_xlsx("~/Dropbox/Encelia/analysis/physiology/encelia_leafdrymass_CGmarch2017.xlsx")
lm$LeafID = NA
lm = data.frame(lm, stringsAsFactors = F)
lm[lm$rep == 1, "LeafID"] = "A"
lm[lm$rep == 2, "LeafID"] = "B"
lm[lm$rep == 3, "LeafID"] = "C"
lm1 = lm %>% select("SampleID" = pos, LeafID, mass)
# lost some data along the way, but hard to figure out data issues
l5b = left_join(l5, lm1)
# leaf area & LMA are good, without leaf mass
l5c = l5b %>% mutate(LMA = mass / Area)

# average across all leaves for an individual
l6 = l5c %>% dplyr::select(-LeafID) %>% 
  group_by(SampleID) %>% summarise_all(mean, na.rm = T) %>% ungroup()
l7 = left_join(l6, cg %>% dplyr::select(COMBINED, PLANT_ID), 
               by = c("SampleID" = "COMBINED"))
l8 = left_join(l7, d %>% dplyr::select(PLANT_ID, PUTATIVE_SPECIES))
# who is missing
l8[!complete.cases(l8$PUTATIVE_SPECIES), ] 
# all these are species that died shortly after photographing
# so not in database! manually fix
l8[which(l8$SampleID == "13U"), "PUTATIVE_SPECIES"] = "Baiylea_multiradiata"
l8[which(l8$SampleID == "9V"), "PUTATIVE_SPECIES"] = "Encelia_farinosa_phenicodonta"
l8[which(l8$SampleID == "3G"), "PUTATIVE_SPECIES"] = "Baiylea_multiradiata"
ll = l8 %>% dplyr::select(-SampleID) %>% tidyr::gather(key = "measurement", 
                          value = "estimate", 
                          -PLANT_ID, -PUTATIVE_SPECIES)
ll$measurement = paste("leaf", ll$measurement, sep = "_")

# stem hydraulic
s = read_xlsx("~/Dropbox/Encelia/analysis/phylogeny/morphological_data/enceliaRiversideKstemLA.xlsx",
              sheet = 4)
# these individual ids map to a previous version of the garden, so not useful here
s$PUTATIVE_SPECIES = sps[s$species]
# all the farinosa here are farinosa farinosa
s[s$species == "farinosa", "PUTATIVE_SPECIES"] = "Encelia_farinosa_farinosa"
ss = s %>% dplyr::select(PUTATIVE_SPECIES, avg_leaf_size_m2, `Kshoot_mmol/s/m2/MPa`) %>% 
  tidyr::gather("measurement", "estimate", -PUTATIVE_SPECIES) %>% 
  tibble::add_column(PLANT_ID = NA, .before = 1)

# stem data
st = read_csv("~/Dropbox/Encelia/analysis/physiology/stem_vessel_diam.csv")
st$ind = gsub("rec\\d+_\\d+_", "", st$scan)
st$ind = gsub(".h5", "", st$ind)
st_inds = c("FAR-1", "FAR-2", "VEN-1", "ASP-1",
            "DEN-1", "ECal_Stem_Dry", "ECan_Stem_Dry",
            "Enceliopsis_1", "PAL-1", "EFru_Stem_Dry",
            "EAct_Stem_Dry", "EFARfar_1")
# not too sure about sp designations for FAR-1 & FAR-2
st_sps = c("Encelia_farinosa_phenicodonta", "Encelia_farinosa_phenicodonta", 
           "Encelia_ventorum", "Encelia_asperifolia",
            "Encelia_densifolia", "Encelia_californica", 
            "Encelia_canescens",
            "Enceliopsis_covillei", "Encelia_palmeri", 
            "Encelia_frutescens_frutescens",
            "Encelia_actoni", "Encelia_farinosa_farinosa")
names(st_sps) = st_inds
st$PUTATIVE_SPECIES = st_sps[st$ind]
st1 = st %>% group_by(PUTATIVE_SPECIES) %>% 
  dplyr::select(-scan, -ind, -vessel_rep) %>%
  summarize_all(mean) %>% 
  ungroup()
st2 = tidyr::gather(st1, "measurement", "estimate", -PUTATIVE_SPECIES) %>% 
  tibble::add_column(PLANT_ID = NA, .before = 1)

#####################
# woody data
#####################
w = read_csv("~/Dropbox/Encelia/analysis/physiology/wooddensity.csv")
sp3 = c("ASP", "DEN", "FAR", "HYB", "PAL", "VEN")
species3 = c("Encelia_asperifolia", "Encelia_densifolia", 
            "Encelia_farinosa_farinosa", 
            "palmeri_ventorum_hybrid", "Encelia_palmeri",
            "Encelia_ventorum")
names(species3) = sp3
w$PUTATIVE_SPECIES = species3[w$species]
w1 = w %>% mutate(measurement = "wood_density")
w2 = w1 %>% select(PUTATIVE_SPECIES, measurement, "estimate" = density) %>% 
  tibble::add_column(PLANT_ID = NA, .before = 1)

dd = do.call("rbind", list(t4, tt, rr, ll, ss, st2, w2))
dd[dd$PUTATIVE_SPECIES == "palmeri_ventorum_F1", "PUTATIVE_SPECIES"] = "palmeri_ventorum_hybrid"
dd[dd$PUTATIVE_SPECIES == "palmeri_ventorum_BC", "PUTATIVE_SPECIES"] = "palmeri_ventorum_hybrid"
dd[dd$PUTATIVE_SPECIES == "asperifolia_ventorum_BC", "PUTATIVE_SPECIES"] = "asperifolia_ventorum_hybrid"
dd[dd$PUTATIVE_SPECIES == "asperifolia_ventorum_F1", "PUTATIVE_SPECIES"] = "asperifolia_ventorum_hybrid"
dd[dd$PUTATIVE_SPECIES == "farinosa_frutescens_BC", "PUTATIVE_SPECIES"] = "farinosa_frutescens_hybrid"
dd[dd$PUTATIVE_SPECIES == "Encelia_farinosa", "PUTATIVE_SPECIES"] = "Encelia_farinosa_farinosa"
dd[dd$PUTATIVE_SPECIES == "Encelia_frutescens", "PUTATIVE_SPECIES"] = "Encelia_frutescens_frutescens"

dd1 = dd[grep("hybrid", dd$PUTATIVE_SPECIES, invert = T), ]

# do a species level summary
dd2 = dd1 %>% dplyr::select(-PLANT_ID) %>% 
  group_by(PUTATIVE_SPECIES, measurement) %>% 
  summarise(mean_estimate = mean(as.numeric(estimate), na.rm = T)) %>% ungroup()
dd3 = dd2 %>% spread(measurement, mean_estimate)

write_csv(dd1, "~/Dropbox/Encelia/analysis/phylogeny/morphological_data/morphological_data.csv")
write_csv(dd3, "~/Dropbox/Encelia/analysis/phylogeny/morphological_data/morphological_data_spread.csv")

# color is correlated with trichome density
cor.test(log(dd3$top_trichomes_per_mm2), dd3$leaf_gray_mean)

# how complete are the data
apply(dd3, 2, function(x) {sum(complete.cases(x))})