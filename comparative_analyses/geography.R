###########################
# geographic range estimation
###########################

library(ggspatial)
library(alphahull)
library(rgeos)
library(readr)
library(tidyr)
library(dplyr)

source('/Users/sonal/Dropbox/scripts/eco_IBD/pascal_scripts/ah2sp.R')

make_polygon <- function(coords) {
  polygon = NA
  
  #if enough points for alpha hulls
  if (nrow(coords) > 3) {
    polygon <-try(getAlphaHullRange(coords, 
                                    percent=1, partCount= 1, 
                                    buff=5000,
                                    coordHeaders=c('Longitude','Latitude'), 
                                    verbose=TRUE), silent=TRUE)
    while ('try-error' %in% class(q)) {
      polygon <-try(getAlphaHullRange(coords, percent=1, 
                                      partCount= 1, buff=5000, 
                                      coordHeaders=c('Longitude','Latitude'), 
                                      verbose=TRUE), silent=TRUE)
    }
    
    #if 3 points, make minimum convex hull
  } else if (nrow(coords) == 3) {
    pt <- SpatialPoints(coords, proj4string=CRS('+proj=longlat +datum=WGS84'))
    pt <- spTransform(pt, CRS("+init=epsg:3395"))
    b <- gBuffer(gConvexHull(pt), width=5000)
    b <- spTransform(b, CRS('+proj=longlat +datum=WGS84'))
    polygon <- list(b, NA, 3)
    
    #if 2 points, make tube polygon
  } else if (nrow(coords) == 2){
    
    polygon <- list(linkTwoPoints(SpatialPoints(coords, 
                                                proj4string=CRS('+proj=longlat +datum=WGS84'))), NA, 2)
    
    # if 1 point, make buffered point
  } else if (nrow(coords) == 1) {
    
    pt <- SpatialPoints(coords, proj4string=CRS('+proj=longlat +datum=WGS84'))
    pt <- spTransform(pt, CRS("+init=epsg:3395"))
    b <- gBuffer(pt, width=5000)
    b <- spTransform(b, CRS('+proj=longlat +datum=WGS84'))
    polygon <- list(b, NA, 1)
  }
  return(polygon)	
}

write_range <- function(range, name, outdir) {
  # save our new ranges!
  filename = paste(outdir, name, sep="")
  if (class(range) == 'SpatialPolygons') {
    IDs <- sapply(slot(range, "polygons"), function(x) slot(x, "ID"))
    df <- data.frame(rep(0, length(IDs)), row.names=IDs)
    writeSpatialShape(SpatialPolygonsDataFrame(range, data=df), filename)
  } else {
    writeSpatialShape(range, filename)
  }
}

d = read_csv("~/Dropbox/Encelia/analysis/spatial_analyses/encelia/all_points_thinned.csv")

# need to make points unique
d1 = unique(d)
# for now, no outgroups
d1 = d1 %>% filter(!species %in% c("Artemisia tridentata", 
                                   "Baileya multiradiata",
                                   "Enceliopsis covillei",
                                   "Xylorhiza tortifolia"))
d1 = as.data.frame(d1) %>% filter(Longitude != -74.37556)
d1$species = gsub(" ", "_", d1$species)
d1[d1$species == "Encelia_actonii", "species"] = "Encelia_actoni"
d1[d1$species == "Encelia_frutescens", "species"] = "Encelia_frutescens_frutescens"

rdata = split(d1, d1$species)
poly = vector("list", length(rdata))
names(poly) = names(rdata)
for (i in 1:length(rdata)) {
  poly[[i]] = make_polygon(rdata[[i]] %>% dplyr::select(Longitude, Latitude)) 
}

areas = rep(NA, length(poly))
names(areas) = names(poly)
pdf("~/Desktop/Encelia_range_maps.pdf")
for (i in 1:length(poly)) {
  plot(poly[[i]][[1]])
  points(rdata[[i]] %>% dplyr::select(Longitude, Latitude))
  areas[i] = gArea(poly[[i]][[1]], byid = T)
  write_range(poly[[i]][[1]], names(areas)[i], "/Users/Sonal/Desktop/encelia_ranges/")
}
dev.off()

get_vals3 <- function(trait, name) {
  par(mar=c(5, 1, 1, 1), tck=-0.02)
  xvals = pretty(range(trait, na.rm = T), 2)
  plot.new()
  plot.window(xlim=range(xvals, na.rm = T), ylim=c(1, length(trait)))
  
  for (i in 1:length(trait)) {
    # lines(x = c(dd[i, varmin], dd[i, varmax]), y = c(i, i))
    points(trait[i], i, pch=21, bg="#377eb8", cex=1)
  }
  
  axis(1, at=xvals, labels=NA, line=1)
  axis(1, at=xvals, labels=xvals, line=0.4, lwd = 0, cex.axis=1)
  mtext(name, side=1, line=3.3, cex=1)
}

area2 = areas[td_nom2$tip.label]
pdf("~/Dropbox/Encelia/manuscripts/Encelia_Phylogeny/figures/Encelia_geography.pdf", width = 3, height = 3)
par(mfrow=c(1, 2))
par(mar=c(5,0,1,0))
plot(td_nom, show.tip.label = F, edge.width=1.5, cex=1.2)
# PC1 is mainly capturing area
get_vals3(log(area2), "log range size")
dev.off()