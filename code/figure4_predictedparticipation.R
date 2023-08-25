library(tidyverse)
library(sf)
library(ggmap)
library(ggsn)
library(viridis)
setwd("~/Optimized tents paper/")

##-----------------------------------------------------------------------------
# 1. Load and clean data
##-----------------------------------------------------------------------------

# Load all possible supply points
supply<-read.csv(file='data/supplypoints.csv',sep=';') %>%
  st_as_sf(., coords = c("long", "lat"), crs = 4326)

# Load all demand points
houses <- read.csv(file='data/demandpoints.csv',sep=';') %>%
  st_as_sf(., coords = c("long", "lat"), crs = 4326)

# Load walking distance matrix
wd <- read.csv("data/distancematrixwalking.csv") %>%
  select(-X)

# Load ASA shapefile (needed to create map borders)
asa <- st_read("data/MR_ALTO_SELVA_ALGRE.kml")

# Indices of vax sites used in 2016 and the different optimized sets
v_2016    <- c(7, 8, 9, 10, 13, 15, 22, 24, 27, 31, 34, 38, 39, 42, 46, 54, 55, 
               63, 66, 67)
v_pcenter <- c(1, 2, 4, 9, 12, 15, 18, 19, 21, 22, 26, 29, 32, 35, 38, 48, 50, 
               51, 54, 61)
v_pmedian <- c(1, 2, 7, 9, 15, 18, 19, 20, 22, 26, 28, 29, 31, 37, 42, 45, 47, 
               57, 61, 63)
v_pprob   <- c(1, 2, 7, 9, 15, 18, 19, 20, 22, 26, 28, 29, 31, 37, 41, 45, 47,
               57, 61, 63)

# Use indices to obtain vaccination locations for the different sets
set_2016 <- supply[v_2016,]
set_center <- supply[v_pcenter,]
set_median <- supply[v_pmedian,]
set_prob <- supply[v_pprob,]

##-----------------------------------------------------------------------------
# 2. Find catchments and estimate probability of participation for each house
##-----------------------------------------------------------------------------

# Get the walking distance matrix wfor the different sets
wd_2016 <- wd[,v_2016]
wd_center <- wd[,v_pcenter]
wd_median <- wd[,v_pmedian]
wd_prob <- wd[,v_pprob]

# Get the minimum walking distance for each house
minwd_2016 <- apply(wd_2016, MARGIN = 1, min)
minwd_center <- apply(wd_center, MARGIN = 1, min)
minwd_median <- apply(wd_median, MARGIN = 1, min)
minwd_prob <- apply(wd_prob, MARGIN = 1, min)

# Rescale the minimum walking distance
minwd_2016_rescaled <- (minwd_2016 - 488.0174)/300.3352 
minwd_center_rescaled <- (minwd_center - 488.0174)/300.3352 
minwd_median_rescaled <- (minwd_median - 488.0174)/300.3352 
minwd_prob_rescaled <- (minwd_prob - 488.0174)/300.3352 

# Use the minimum walking distance and the mixed effects Poisson decay function 
# for 2016 to find prob vax for each house
probvax_2016 <- exp(-0.4909919 - 0.1261808*minwd_2016_rescaled)
probvax_center <- exp(-0.4909919 - 0.1261808*minwd_center_rescaled)
probvax_median <- exp(-0.4909919 - 0.1261808*minwd_median_rescaled)
probvax_prob <- exp(-0.4909919 - 0.1261808*minwd_prob_rescaled)
mean(probvax_2016)  # mean = 0.5910666
mean(probvax_center)  # mean = 0.6216073
mean(probvax_median)  # mean = 0.6321322
mean(probvax_prob)  # mean = 0.6321337

# Use the minimum walking distance and the mixed effects Poisson decay function 
# for 2016 to find prob vax for each house
houses$probvax_2016 <- probvax_2016
houses$probvax_center <- probvax_center
houses$probvax_median <- probvax_median
houses$probvax_prob <- probvax_prob

##-----------------------------------------------------------------------------
# 3. Create map of houses colored by probability of participation
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

# Find colors for mapping
precols <- plasma(50)
mycols <- precols[c(1:25, seq(26, 50, by = 2))]

# Find probvax limits for mapping
summary(houses$probvax_2016)
summary(houses$probvax_center)
summary(houses$probvax_median)
summary(houses$probvax_prob)

pdf("figures/R_output/figure4_predictedparticipation.pdf", height = 6, width = 7)

# Create map for 2016 sites
basemap +
  geom_sf(data = houses, aes(color = probvax_2016), alpha = 0.5, inherit.aes = FALSE) +
  scale_color_gradientn(colours = mycols, limits = c(0.174, 0.752)) +
  geom_sf(data = set_2016, inherit.aes = FALSE, size = 2, color = "black", shape = 24,
          fill = "white") +
  ggtitle("2016 sites") +
  scalebar(location = "bottomright", data = asa, dist_unit = "km",
           dist = 1, transform = TRUE, model = 'WGS84') 
  
# Create map for p-center sites
basemap +
  geom_sf(data = houses, aes(color = probvax_center), alpha = 0.5, inherit.aes = FALSE) +
  scale_color_gradientn(colours = mycols, limits = c(0.174, 0.752)) +
  geom_sf(data = set_center, inherit.aes = FALSE, size = 2, color = "black", shape = 24,
          fill = "white") +
  ggtitle("p-center sites") +
  scalebar(location = "bottomright", data = asa, dist_unit = "km",
           dist = 1, transform = TRUE, model = 'WGS84') 

# Create map for p-median sites
basemap +
  geom_sf(data = houses, aes(color = probvax_median), alpha = 0.5, inherit.aes = FALSE) +
  scale_color_gradientn(colours = mycols, limits = c(0.174, 0.752)) +
  geom_sf(data = set_median, inherit.aes = FALSE, size = 2, color = "black", shape = 24,
          fill = "white") +
  ggtitle("p-median sites") +
  scalebar(location = "bottomright", data = asa, dist_unit = "km",
           dist = 1, transform = TRUE, model = 'WGS84') 
  

# Create map for p-probability sites 
basemap +
  geom_sf(data = houses, aes(color = probvax_prob), alpha = 0.5, inherit.aes = FALSE) +
  scale_color_gradientn(colours = mycols, limits = c(0.174, 0.752)) +
  geom_sf(data = set_prob, inherit.aes = FALSE, size = 2, color = "black", shape = 24,
          fill = "white") +
  ggtitle("p-probability sites") +
  scalebar(location = "bottomright", data = asa, dist_unit = "km",
           dist = 1, transform = TRUE, model = 'WGS84') 

dev.off()
