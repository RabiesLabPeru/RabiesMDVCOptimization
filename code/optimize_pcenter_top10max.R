library(sf)
library(tidyverse)
library(analyze.stuff)
library(tictoc)

##-----------------------------------------------------------------------------
# 0. User-defined function
##-----------------------------------------------------------------------------

# Find the maximum walking distance between any house and its closest 
# vaccination site; also keep information about house and vaccination tent
# that is the largest distance apart

FindMaxDist <- function(dm){
  
  # Keep only shortest walking distance for each house
  row_mins <- apply(dm, 1, which.min) 
  dm_shortest <- dm
  for (i in 1:nrow(dm_shortest)){
    dm_shortest[i, -row_mins[i]] <- NA
    i <- i+1
  }
  # Find the maximum shortest walking distance for each vaccination tent and
  # sum the 10 greatest values
  max_wd <- sort(apply(dm_shortest, MARGIN = 2, max, na.rm = T))
  top10max <- sum(max_wd[11:20])
  
  return(top10max)
}

#Function solves the p-center problem 
SolvePCenter <- function(dmatrix, p) {
  n_supply <- ncol(dmatrix)
  
  # Generate random sample of p supply points 
  supplyset <- sample(c(1:n_supply), p)
  dm_set <- dmatrix[, supplyset]  # Distance matrix with only chosen supply points
  
  # Find maximum distance between any house and its closest vax site
  top10max <- FindMaxDist(dm_set)
  print(paste0(paste(c("Supply points:", supplyset), collapse = " "), 
               " Score=", top10max))
  
  # For loop across each selected site
  for (i in 1:length(supplyset)){
    othersites <- c(1:n_supply)[! c(1:n_supply) %in% supplyset]
    
    # For loop for the ith site to rotate through all remaining site
    for (j in othersites){
      supplyset_tmp <- supplyset
      supplyset_tmp[i] <- j  # Replace ith site with site j
      dm_tmp <- dmatrix[, supplyset_tmp]
      top10max_tmp <- FindMaxDist(dm_tmp)
      if(top10max_tmp < top10max){  # If replacing ith site with j improves score then replace
        top10max <- top10max_tmp
        supplyset <- supplyset_tmp
      }
    }
    print(paste0(paste(c("Supply points:", supplyset), collapse = " "), 
                 " Score=", top10max))
  }
  
  return(list(supply = supplyset, top10max = top10max))
}

##-----------------------------------------------------------------------------
# 1. Load and clean data 
##-----------------------------------------------------------------------------

setwd("~/Optimized tents paper/")

# Demand points (houses)
# Load houses, convert to sf, and re-project to UTM coordinates
# Load demand points (houses) and convert to sf
houses <- read.csv(file='data/demandpoints.csv',sep=';') %>%
  st_as_sf(., coords = c("long", "lat"), crs = 4326) %>%
  st_transform(29188)

# Supply points (potential vaccination sites)
supply<-read.csv(file='data/supplypoints.csv',sep=';') %>%
  st_as_sf(., coords = c("long", "lat"), crs = 4326) %>%
  st_transform(29188)

# Distance matrix
distancematrix <- read.csv("data/distancematrixwalking.csv") %>%
  select(-X)
distancematrix <- as.matrix(distancematrix)
colnames(distancematrix) <- supply$VaccPoint
rownames(distancematrix) <- houses$UNICODE

# Load ASA microred shapefile for mapping
asa <- st_read("data/MR_ALTO_SELVA_ALGRE.kml") %>%
  st_transform(29188)

##-----------------------------------------------------------------------------
# 2. Run tb optimization algorithm
##-----------------------------------------------------------------------------

set.seed(123)

# Get 100 sets of solutions using 100 different seeds
rand_seed <- sample(999999, 1000)
solutions <- list()

tic()
for (i in 1:1000){
  set.seed(rand_seed[i])
  solutions[[i]] <- SolvePCenter(distancematrix, p = 20)
}
toc()

scores <- numeric()
supplysets <- list()
for (i in 1:1000){
  scores[i] <- solutions[[i]]$top10max
  supplysets[[i]] <- solutions[[i]]$supply
}

# Save solutions
save(solutions, scores, supplysets, file = "pcenter_solutions.rda")

##-----------------------------------------------------------------------------
# 3. Obtain best solution
##-----------------------------------------------------------------------------

load("output/pcenter_solutions.rda")

length(unique(scores))  # 817
hist(scores)
which.min(scores)  # 430

c(1, 2, 4, 9, 12, 15, 18, 19, 21, 22, 26, 29, 32, 35, 38, 48, 50, 
  51, 54, 61)

sort(supplysets[[430]]) ## REPLACE WITH INDEX
#  [1]  1 2 4 9 12 15 18 19 21 22 26 29 32 35 38 48 50 51 54 61
supplysets_pcenter <- c(1, 2, 4, 9, 12, 15, 18, 19, 21, 22, 26, 29, 32, 35, 38, 
                        48, 50, 51, 54, 61)
