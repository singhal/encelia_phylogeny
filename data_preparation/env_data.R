library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggmap)
library(maps)
library(raster)
library(geosphere)
library(CoordinateCleaner)
library(rnaturalearth)
library(rnaturalearthdata)
library(mapdata)
library(spThin)

w2lr <- map_data("worldLores")
w2hr <- map_data("worldHires")
outdir = "/Users/sonal/Desktop/"

d = list.files("~/Desktop/activeWork/spatial_analyses/", pattern = "*csv",
               full.names = T)
d = lapply(d, function(x) {read.csv(x, sep="\t", 
                                    stringsAsFactors = F, na.string = c("NA", ""))})
d = do.call(rbind, d)
d1 = d[complete.cases(d$decimalLatitude), ]
d1 = d1[which(d1$basisOfRecord %in% c('PRESERVED_SPECIMEN', 'HUMAN_OBSERVATION')), ]

# remove imprecise points
hist(d1$coordinateUncertaintyInMeters / 1000, breaks = 20)
d1 = d1 %>% filter(coordinateUncertaintyInMeters / 1000 <= 100 | is.na(coordinateUncertaintyInMeters))

# keep only species of interest
sps = c("Encelia actonii", "Encelia asperifolia", "Encelia californica",
        "Encelia canescens", "Encelia densifolia",
        "Encelia farinosa", "Encelia frutescens", 
        "Encelia palmeri", "Encelia ravenii",
        "Encelia resinifera", "Encelia ventorum", "Encelia virginensis",
        "Artemisia tridentata", "Baileya multiradiata",
        "Enceliopsis covillei", "Xylorhiza tortifolia")
d2 = d1 %>% filter(species %in% sps)

rownames(d2)<-1:nrow(d2)
flags = clean_coordinates(x = d2, 
                          lon = "decimalLongitude", lat = "decimalLatitude",
                          countries = "countryCode", 
                          species = "species",
                          tests = c("capitals", "centroids",
                                    "equal","gbif", "institutions", "seas",
                                    "zeros"))
d2$flagged1 = !(flags$.summary)

sp_d = split(d2, d2$species)
sp_d_thin = vector('list', length(sp_d))
sensitivity = rep(2, length(sp_d))
names(sensitivity) = names(sp_d)
# need to manually adjust the sensitivity of some of these 
# species because the default doesn't work for all of them
sensitivity["Encelia asperifolia"] = 0.01
sensitivity["Encelia californica"] = 5
sensitivity["Encelia farinosa"] = 5
sensitivity["Encelia palmeri"] = 1
sensitivity["Encelia ventorum"] = 0.01
sensitivity["Encelia virginensis"] = 0.01

for (i in 1:length(sp_d)) {
  sp_pts = sp_d[[i]]
  rownames(sp_pts) = 1:nrow(sp_pts)
  outl = cc_outl(x = sp_pts, lon = "decimalLongitude", 
               lat = "decimalLatitude", mltpl = sensitivity[i],
               value = "flagged")
  sp_pts = data.frame(sp_pts, outlier =  as.factor(!outl))
  
  sp_ptsx = sp_pts[sp_pts$outlier == FALSE, ]
  sp_pts2 = distinct(sp_ptsx[, c("species", "decimalLatitude", "decimalLongitude")])
  sp_pts3 = thin(sp_pts2, "decimalLatitude", "decimalLongitude", "species", 
       thin.par = 1, reps = 1, 
       locs.thinned.list.return = T, write.files = F)[[1]]
  sp_d_thin[[i]] = sp_pts3
  sp_d[[i]] = sp_pts
}

d3 = do.call("rbind", sp_d)
rownames(d3) = 1:nrow(d3)
d3$flagged = rep(NA, nrow(d3))
d3[(d3$flagged1 == FALSE & d3$outlier == FALSE), "flagged"] = "GOOD"
d3[(d3$flagged1 == TRUE & d3$outlier == FALSE), "flagged"] = "BAD_CLEAN"
d3[(d3$flagged1 == TRUE & d3$outlier == TRUE), "flagged"] = "BAD_BOTH"
d3[(d3$flagged1 == FALSE & d3$outlier == TRUE), "flagged"] = "BAD_GEO"

plt = vector("list", length(sps))
names(plt) = sps

for (i in 1:length(sps)) {
  d4 = d3[which(d3$species == sps[i]), ]
  d4 = d4[d4$flagged %in% c("GOOD", "BAD_GEO"), ]
  d4$flagged = factor(d4$flagged, levels = c("GOOD", "BAD_GEO"))
  
  xlim1 = range(d4$decimalLongitude) * c(1.05, 0.95)
  ylim1 = ifelse(min(d4$decimalLatitude) > 0, 
                 min(d4$decimalLatitude) * 0.9, 
                 min(d4$decimalLatitude) * 1.1)
  ylim2 = ifelse(max(d4$decimalLatitude) > 0, 
                 max(d4$decimalLatitude) * 1.1, 
                 max(d4$decimalLatitude) * 0.9)
  
  a = ggplot() +
    geom_polygon(data = w2lr, fill = "grey90", 
                 aes(x = long, y= lat, group = group)) +
    theme_bw() +
    coord_fixed(xlim = xlim1,  ylim = c(ylim1, ylim2), ratio = 1.3)
  plt[[i]] = a + geom_point(data = d4, size = 0.4, 
                        aes(x = decimalLongitude, y = decimalLatitude, 
                        col = flagged)) +
      scale_colour_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")) + 
      theme(legend.position = "bottom")
}

for (i in 1:length(plt)) {
  save_plot(paste("~/Desktop/activeWork/spatial_analyses/encelia/", gsub(" ", "_", sps[i]), ".pdf", sep = ""),
            plt[[i]])
}
d4 = d3[which(d3$flagged == 'GOOD'), ]
write.csv(d4, "~/Desktop/activeWork/spatial_analyses/encelia/all_points.csv", row.names = F)
for (i in 1:length(sp_d)) {
  sp_d_thin[[i]]$species = names(sp_d)[i]
}
d4_thin = do.call("rbind", sp_d_thin)
write.csv(d4_thin, "~/Desktop/activeWork/spatial_analyses/encelia/all_points_thinned.csv", row.names = F)