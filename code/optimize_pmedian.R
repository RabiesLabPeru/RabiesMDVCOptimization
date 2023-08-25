library(tbart)
library(sf)
library(tidyverse)

##-----------------------------------------------------------------------------
# 0. User-defined function
##-----------------------------------------------------------------------------

GetScores <- function(opt_vector, dm){
  
  # Update distance matrix keeping only the optimized supply points
  dm_opt <- dm[,opt_vector]
  
  # Find the walking distance to nearest supply point for each house
  row_mins <- rowMins(dm_opt)
  
  # Sum the total walking distance
  tot_dist <- sum(row_mins)
  
  return(tot_dist)
  
}

##-----------------------------------------------------------------------------
# 1. Load and clean data 
##-----------------------------------------------------------------------------

setwd("~/Optimized tents paper/")

# Demand points (houses)
# Load houses, convert to sf, and re-project to UTM coordinates
houses<-read.csv(file='data/demandpoints.csv',sep=';') %>%
  st_as_sf(., coords = c("long", "lat"), crs = 4326) %>%
  st_transform(29188)
demandcoord <- st_coordinates(houses)

# Supply points (potential vaccination sites)
supply<-read.csv(file='data/supplypoints.csv',sep=';') %>%
  st_as_sf(., coords = c("long", "lat"), crs = 4326) %>%
  st_transform(29188)
supplycoord <- st_coordinates(supply)

# Distance matrix
distancematrix <- read.csv("data/distancematrixwalking.csv") %>%
  select(-X)
distancematrix <- as.matrix(distancematrix)

# Load ASA microred shapefile for mapping
asa <- st_read("data/MR_ALTO_SELVA_ALGRE.kml") %>%
  st_transform(29188)
  
##-----------------------------------------------------------------------------
# 2. Run tb optimization algorithm
##-----------------------------------------------------------------------------

set.seed(123)

# Get 1000 sets of solutions using 1000 different seeds
rand_seed <- sample(999999, 1000)
solutions <- list()

for (i in 1:1000){
  set.seed(rand_seed[i])
  solutions[[i]] <- tb(demandcoord, supplycoord, p=20, metric=distancematrix,
                       verbose=TRUE)
}

scores <- numeric()
library(analyze.stuff)
for (i in 1:1000){
  set.seed(i)
  scores[i] <- GetScores(solutions[[i]], distancematrix)
}

# Save solutions
save(solutions, scores, avgwd, file = "pmedian_solutions.rda")

##-----------------------------------------------------------------------------
# 3. Obtain best solution
##-----------------------------------------------------------------------------

load("output/pmedian_solutions.rda")

length(unique(scores))  # 15
hist(scores)
which.min(scores)  # 4

solutions[[4]]  
sort(solutions[[4]])
#  [1] 1  2  7  9 15 18 19 20 22 26 28 29 31 37 42 45 47 57 61 63
supplysets_pmedian   <- c(1, 2, 7, 9, 15, 18, 19, 20, 22, 26, 28, 29, 31, 37, 
                          42, 45, 47, 57, 61, 63)
