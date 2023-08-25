library(tidyverse)
library(sf)
library(ggmap)
library(ggsn)
library(viridis)
library(scales)
setwd("~/Optimized tents paper/")

##-----------------------------------------------------------------------------
# 1. Load data
##-----------------------------------------------------------------------------

# Load demand points (houses) and convert to sf
houses <- read.csv(file='data/demandpoints.csv',sep=';') %>%
  st_as_sf(., coords = c("long", "lat"), crs = 4326)

# Load supply points (potential vax sites)
supply <- read.csv(file='data/supplypoints.csv',sep=';') %>%
  st_as_sf(., coords = c("long", "lat"), crs = 4326)

# Indices of vax sites used in 2016, p-center solution, p-median solution, and
# p-probability solution
v_16 <- c(7, 8, 9, 10, 13, 15, 22, 24, 27, 31, 34, 38, 39, 42, 46, 54, 55, 63,
         66, 67)
v_pcenter <- c(1, 2, 4, 9, 12, 15, 18, 19, 21, 22, 26, 29, 32, 35, 38, 48, 50, 
               51, 54, 61)
v_pmedian <- c(1, 2, 7, 9, 15, 18, 19, 20, 22, 26, 28, 29, 31, 37, 42, 45, 47, 
               57, 61, 63)
v_pprob   <- c(1, 2, 7, 9, 15, 18, 19, 20, 22, 26, 28, 29, 31, 37, 41, 45, 47,
               57, 61, 63)

# Create supply sets for comparison
s_2016 <- supply[v_16,]
s_center <- supply[v_pcenter,]
s_median <- supply[v_pmedian,]
s_prob   <- supply[v_pprob,]

# Load walking distance matrix for all supply points and make subsets for each
# comparison group
dm <- read.csv("data/distancematrixwalking.csv") %>%
  select(-X)
colnames(dm) <- supply$VaccPoint
rownames(dm) <- houses$UNICODE
dm_2016 <- dm[,v_16]
dm_center <- dm[,v_pcenter]
dm_median <- dm[,v_pmedian]
dm_pprob  <- dm[,v_pprob]

# Load ASA shapefile
asa <- st_read("data/MR_ALTO_SELVA_ALGRE.kml")

# Create basemap for ASA 
bbox <- st_bbox(houses)
pad <- 0.004
pad2 <- 0.0025
map_borders <- c(bottom = as.numeric(bbox$ymin) - pad, 
                 top = as.numeric(bbox$ymax) + pad, 
                 left = as.numeric(bbox$xmin) - pad2, 
                 right = as.numeric(bbox$xmax) + pad2)

# Get basemap -- Stamen terrain 
asa_terrain <- get_stamenmap(bbox = map_borders, zoom = 14, maptype = "terrain")
basemap <- ggmap(asa_terrain)

##-----------------------------------------------------------------------------
# 2. Make catchment map for 2016 vaccination sites
##-----------------------------------------------------------------------------

# Get vax site ID's for 2016 vax campaign
s_names <- colnames(dm_2016)

# Find closest vax site for each house in 2016 campaign and save info using vax
# site ID's
assignment_byhouse <- apply(dm_2016, MARGIN = 1, which.min)
s_names <- s_names[unique(assignment_byhouse)]  # Remove 1 site with no assignments
assignment_names <- factor(assignment_byhouse, labels = s_names, ordered = T)
table(assignment_names)  # Check labels assigned as expected

# Add variable to houses sf with site assignments
houses$assignment_2016 <- assignment_names

# Create color palette
cols_plasma <- rep(plasma(8), 3)[c(1:16, 1, 3, 5, 8)]
set.seed(1)
cols <- sample(cols_plasma, 20, replace = F)

# Tweak colors
cols[19] <- cols_plasma[7]
cols[3] <- cols_plasma[8]

# Create map of catchments in 2016
pdf("figures/R_output/figure3_2016catchments.pdf", height = 6, width = 6.5)
basemap +
  geom_sf(data = houses, aes(col = assignment_2016), inherit.aes = FALSE,
          size = 1, alpha = 0.3) +
  scale_color_manual(values= cols) +
  geom_sf(data = s_2016, fill = "white", shape = 24, size = 3,
          color = "black", inherit.aes = FALSE) +
  scalebar(location = "bottomright", data = houses, dist_unit = "km",
           dist = 1, transform = TRUE, model = 'WGS84')
dev.off()

##-----------------------------------------------------------------------------
# 3. Make catchment map for p-center solution
##-----------------------------------------------------------------------------

# Get vax site ID's for 2016 vax campaign
s_names <- colnames(dm_center)

# Find closest vax site for each house in 2016 campaign and save info using vax
# site ID's
assignment_byhouse <- apply(dm_center, MARGIN = 1, which.min)
assignment_names <- factor(assignment_byhouse, labels = s_names, ordered = T)
table(assignment_names)  # Check labels assigned as expected

# Add variable to houses sf with site assignments
houses$assignment_pcenter <- assignment_names

# Create color palette
cols_plasma <- rep(plasma(8), 3)[c(1:16, 1, 3, 5, 8)]
set.seed(17)  # 3 was good
cols <- sample(cols_plasma, 20, replace = F)

# Tweak colors
cols[13] <- cols_plasma[7]
cols[3] <- cols_plasma[4]
cols[6] <- cols_plasma[2]

# Create map of p-center solution
pdf("figures/R_output/figure3_pcentercatchments.pdf", height = 6, width = 6.5)
basemap +
  geom_sf(data = houses, aes(col = assignment_pcenter), inherit.aes = FALSE,
          size = 1, alpha = 0.3) +
  scale_color_manual(values= cols) +
  geom_sf(data = s_center, fill = "white", shape = 24, size = 3,
          color = "black", inherit.aes = FALSE) +
  scalebar(location = "bottomright", data = houses, dist_unit = "km",
           dist = 1, transform = TRUE, model = 'WGS84')
dev.off()


##-----------------------------------------------------------------------------
# 4. Make catchment map for p-median solution
##-----------------------------------------------------------------------------

# Get vax site ID's for p-median solution
s_names <- colnames(dm_median)

# Find closest vax site for each house in 2016 campaign and save info using vax
# site ID's
assignment_byhouse <- apply(dm_median, MARGIN = 1, which.min)
assignment_names <- factor(assignment_byhouse, labels = s_names, ordered = T)
table(assignment_names)  # Check labels assigned as expected

# Add variable to houses sf with site assignments
houses$assignment_pmedian <- assignment_names

# Create color palette
cols_plasma <- plasma(8)
cols_sample <- c(rep(cols_plasma, 2), cols_plasma[c(1,3,5,7)])
set.seed(5)
cols <- sample(cols_sample)

# Create map of catchments in 2016
pdf("figures/R_output/figure3_pmediancatchments.pdf", height = 6, width = 6.5)
basemap +
  geom_sf(data = houses, aes(col = assignment_pmedian), inherit.aes = FALSE,
          size = 1, alpha = 0.3) +
  scale_color_manual(values= cols) +
  geom_sf(data = s_median, fill = "white", shape = 24, size = 3,
          color = "black", inherit.aes = FALSE) +
  scalebar(location = "bottomright", data = houses, dist_unit = "km",
           dist = 1, transform = TRUE, model = 'WGS84')
dev.off()

##-----------------------------------------------------------------------------
# 5. Make catchment map for p-probability solution
##-----------------------------------------------------------------------------

# Get vax site ID's for p-probability solution
s_names <- colnames(dm_pprob)

# Find closest vax site for each house in 2016 campaign and save info using vax
# site ID's
assignment_byhouse <- apply(dm_pprob, MARGIN = 1, which.min)
assignment_names <- factor(assignment_byhouse, labels = s_names, ordered = T)
table(assignment_names)  # Check labels assigned as expected

# Add variable to houses sf with site assignments
houses$assignment_pprob <- assignment_names

# Create color palette
cols_plasma <- plasma(8)
cols_sample <- c(rep(cols_plasma, 2), cols_plasma[c(1,3,5,7)])
set.seed(5)
cols <- sample(cols_sample)

# Tweak colors
cols[15] <- cols_plasma[8]

# Create map of catchments in 2016
pdf("figures/R_output/figure3_pprobcatchments.pdf", height = 6, width = 6.5)
basemap +
  geom_sf(data = houses, aes(col = assignment_pprob), inherit.aes = FALSE,
          size = 1, alpha = 0.3) +
  scale_color_manual(values= cols) +
  geom_sf(data = s_pprob, fill = "white", shape = 24, size = 3,
          color = "black", inherit.aes = FALSE) +
  scalebar(location = "bottomright", data = houses, dist_unit = "km",
           dist = 1, transform = TRUE, model = 'WGS84')
dev.off()


##-----------------------------------------------------------------------------
# 6. Create histograms showing each house's distance to closest vax site for 
#    every set of sites (2016, p-center, p-median, and p-probability)
##-----------------------------------------------------------------------------

# Find distance to closest vax site for 2016 sites
distance_byhouse <- apply(dm_2016, MARGIN = 1, min)
houses$distance_2016 <- distance_byhouse

# Find distance to closest vax site for p-center solution
distance_byhouse <- apply(dm_center, MARGIN = 1, min)
houses$distance_pcenter<- distance_byhouse

# Find distance to closest vax site for p-median solution
distance_byhouse <- apply(dm_median, MARGIN = 1, min)
houses$distance_pmedian <- distance_byhouse

# Find distance to closest vax site for p-probability solution
distance_byhouse <- apply(dm_pprob, MARGIN = 1, min)
houses$distance_pprob <- distance_byhouse

pdf("figures/R_output/figure3_distancehistograms.pdf", height = 5, width = 7)
# Create histogram for 2016
ggplot(houses) +
  geom_histogram(aes(x = distance_2016), binwidth = 100, fill = "grey", col = "black", boundary = T) +
  xlab("Distance to closest site, m") +
  ylab("Number of houses") +
  ggtitle("2016 sites") +
  scale_x_continuous(breaks = pretty_breaks(), limits = c(0, 2100)) +
  scale_y_continuous(breaks = seq(0, 5000, by = 1000), 
                     minor_breaks = seq(0, 5500, by = 500),
                                        limits = c(0, 5500)) +
  theme_minimal() +
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=28, face = "bold")) 

ggplot(houses) +
  geom_histogram(aes(x = distance_pcenter), binwidth = 100, fill = "grey", col = "black", boundary = T) +
  xlab("Distance to closest site, m") +
  ylab("Number of houses") +
  ggtitle("P-center sites") +
  scale_x_continuous(breaks = pretty_breaks(), limits = c(0, 2100)) +
  scale_y_continuous(breaks = seq(0, 5000, by = 1000), 
                     minor_breaks = seq(0, 5500, by = 500),
                     limits = c(0, 5500)) +
  theme_minimal() +
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=28, face = "bold")) 

# Create histogram for p-median
ggplot(houses) +
  geom_histogram(aes(x = distance_pmedian), binwidth = 100, fill = "grey", col = "black", boundary = T) +
  xlab("Distance to closest site, m") +
  ylab("Number of houses") +
  ggtitle("P-median sites") +
  scale_x_continuous(breaks = pretty_breaks(), limits = c(0, 2100)) +
  scale_y_continuous(breaks = pretty_breaks(), limits = c(0, 5500)) +
  theme_minimal() +
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=28, face = "bold")) 

# Create histogram for p-probability
ggplot(houses) +
  geom_histogram(aes(x = distance_pprob), binwidth = 100, fill = "grey", col = "black", boundary = T) +
  xlab("Distance to closest site, m") +
  ylab("Number of houses") +
  ggtitle("P-probability sites") +
  scale_x_continuous(breaks = pretty_breaks(), limits = c(0, 2100)) +
  scale_y_continuous(breaks = pretty_breaks(), limits = c(0, 5500)) +
  theme_minimal() +
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=28, face = "bold")) 
dev.off()

##-----------------------------------------------------------------------------
# 7. Create histograms showing workload at each site
##-----------------------------------------------------------------------------

pdf("figures/R_output/figure3_workloadhistograms.pdf", height = 5, width = 7)

workload_2016 <- st_drop_geometry(houses) %>%
  group_by(assignment_2016) %>%
  summarise(count = n()) %>%
  arrange(desc(count))
workload_2016$label <- as.factor(1:19)

ggplot(workload_2016) +
  geom_bar(aes(x = label, y = count), stat = "identity", width = 1, 
           fill = "grey", col = "black") +
  xlab("Catchment") +
  ylab("Number of houses") +
  ggtitle("2016 sites") +
  scale_y_continuous(breaks = pretty_breaks(), limits = c(0, 3600)) +
  theme_minimal() +
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=28, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

workload_pcenter <- st_drop_geometry(houses) %>%
  group_by(assignment_pcenter) %>%
  summarise(count = n()) %>%
  arrange(desc(count))
workload_pcenter$label <- as.factor(1:20)

ggplot(workload_pcenter) +
  geom_bar(aes(x = label, y = count), stat = "identity", width = 1, 
           fill = "grey", col = "black") +  xlab("Catchment") +
  ylab("Number of houses") +
  ggtitle("P-center sites") +
  scale_y_continuous(breaks = pretty_breaks(), limits = c(0, 3600)) +
  theme_minimal() +
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=28, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

workload_pmedian <- st_drop_geometry(houses) %>%
  group_by(assignment_pmedian) %>%
  summarise(count = n()) %>%
  arrange(desc(count))
workload_pmedian$label <- as.factor(1:20)

ggplot(workload_pmedian) +
  geom_bar(aes(x = label, y = count), stat = "identity", width = 1, 
           fill = "grey", col = "black") +  xlab("Catchment") +
  ylab("Number of houses") +
  ggtitle("P-median sites") +
  scale_y_continuous(breaks = pretty_breaks(), limits = c(0, 3600)) +
  theme_minimal() +
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=28, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

workload_pprob <- st_drop_geometry(houses) %>%
  group_by(assignment_pprob) %>%
  summarise(count = n()) %>%
  arrange(desc(count))
workload_pprob$label <- as.factor(1:20)

ggplot(workload_pprob) +
  geom_bar(aes(x = label, y = count), stat = "identity", width = 1, 
           fill = "grey", col = "black") +  xlab("Catchment") +
  ylab("Number of houses") +
  ggtitle("P-probability sites") +
  scale_y_continuous(breaks = pretty_breaks(), limits = c(0, 3600)) +
  theme_minimal() +
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=28, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

dev.off()
