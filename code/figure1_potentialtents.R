library(tidyverse)
library(sf)
library(ggmap)
#devtools::install_github('oswaldosantos/ggsn')
library(ggsn)

setwd("~/Optimized tents paper/")

##-----------------------------------------------------------------------------
# 1. Load and clean data 
##-----------------------------------------------------------------------------

# Load ASA microred shapefile
asa <- st_read("data/MR_ALTO_SELVA_ALGRE.kml")
plot(st_geometry(asa))

# Load ASA district shapefile -- takeaway is that ASA 
asa_dist <- st_read("data/dist_AQPcity.shp") %>%
  filter(dist == "ALTO SELVA ALEGRE")
plot(st_geometry(asa_dist))
plot(st_geometry(asa), col = "red", add = T)


# Load supply points (potential vaccination sites)
supply<-read.csv(file='data/supplypoints.csv',sep=';') %>%
  st_as_sf(., coords = c("long", "lat"), crs = 4326)
plot(st_geometry(supply), add = T)

##-----------------------------------------------------------------------------
# 2. Create map of ASA with potential vaccination sites
##-----------------------------------------------------------------------------

# Create basemap for ASA 
bbox <- st_bbox(asa)
pad <- 0.005
pad2 <- 0.01
map_borders <- c(bottom = as.numeric(bbox$ymin) - pad, 
                 top = as.numeric(bbox$ymax) + pad, 
                 left = as.numeric(bbox$xmin) - pad2, 
                 right = as.numeric(bbox$xmax) + pad2)

# Get basemap -- Stamen terrain 
asa_terrain <- get_stamenmap(bbox = map_borders, zoom = 14, maptype = "terrain")
basemap <- ggmap(asa_terrain)

# Save different color versions

pdf("figures/R_output/figure1_potentialtents_white.pdf", height = 6, width = 5.5)
basemap +
  geom_sf(data = asa, inherit.aes = FALSE, fill = NA, color = "black", 
          lwd = 1) +
  geom_sf(data = supply, inherit.aes = FALSE, size = 2, color = "black", shape = 24,
          fill = "white") +
  scalebar(location = "bottomright", data = asa, dist_unit = "km",
           dist = 1, transform = TRUE, model = 'WGS84')
dev.off()

pdf("figures/R_output/figure1_potentialtents_lightred.pdf", height = 6, width = 5.5)
basemap +
  geom_sf(data = asa, inherit.aes = FALSE, fill = NA, color = "black", 
          lwd = 1) +
  geom_sf(data = supply, inherit.aes = FALSE, size = 2, color = "black", shape = 24,
          fill = "#f7cdcd") +
  scalebar(location = "bottomright", data = asa, dist_unit = "km",
              dist = 1, transform = TRUE, model = 'WGS84')
dev.off()

pdf("figures/R_output/figure1_potentialtents_darkred.pdf", height = 6, width = 5.5)
basemap +
  geom_sf(data = asa, inherit.aes = FALSE, fill = NA, color = "black", 
          lwd = 1) +
  geom_sf(data = supply, inherit.aes = FALSE, size = 2, color = "black", shape = 24,
          fill = "#bf4141") +
  scalebar(location = "bottomright", data = asa, dist_unit = "km",
           dist = 1, transform = TRUE, model = 'WGS84')
dev.off()


