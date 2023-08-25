library(tidyverse)
library(sf)
library(tictoc)
library(analyze.stuff)

setwd("~/Optimized tents paper/")

##-----------------------------------------------------------------------------
# 0. User-defined functions
##-----------------------------------------------------------------------------

CalculateExpProb <- function(pm, supplyset = NULL){
  if (is.null(supplyset) == FALSE){
    pm <- pm[, supplyset]
  }
  # Find the maximum participation probability for each demand point
  row_maxs <- rowMaxs(pm)
  
  # Sum the probabilities and find average to get expected vax coverage
  tot_prob <- sum(row_maxs)
  exp_prob <- tot_prob/nrow(pm)
  return(exp_prob)
}


# Running optimizer with maximizing participation probability as the objective
# function
OptimizeParticipationProbability <- function(pmatrix, p, forced = NULL) {
  n_supply <- ncol(pmatrix)
  
  # Generate random sample of p supply points 
  supplyset <- sample(c(1:n_supply), p)
  pm_set <- pmatrix[, supplyset]  # Probability matrix with only chosen supply points
  
  exppar <- CalculateExpProb(pm_set)  # Find and print expected participation
  print(paste(paste(c("Supply points:", supplyset), collapse = " "), "Exp. coverage:", exppar))
  
  # For loop across each selected site
  for (i in (length(forced) + 1):length(supplyset)){
    othersites <- c(1:n_supply)[! c(1:n_supply) %in% supplyset]
    
    # For loop for the ith site to rotate through all remaining site
    for (j in othersites){
      supplyset_tmp <- supplyset
      supplyset_tmp[i] <- j  # Replace ith site with site j
      pm_tmp <- pmatrix[, supplyset_tmp]
      exppar_tmp <- CalculateExpProb(pm_tmp)
      if(exppar_tmp > exppar){  # If replacing ith site with j improves score then replace
        exppar <- exppar_tmp
        supplyset <- supplyset_tmp
      }
    }
    print(paste(paste(c("Supply points:", supplyset), collapse = " "), "Exp. coverage:", exppar))
  }
  
  supplynames <- colnames(pmatrix)[supplyset]
  
  return(list(supply = supplyset, supplynames = supplynames, score = exppar))
}

FindCatchment <- function(dm){
  size_houses <- apply(dm, MARGIN = 1, which.min)
  size_hhwdogs 
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

# Convert distances to probabilities using the mixed-effects Poisson regression
# function
# -- First we'll need to rescale the distances
rescaled_distancematrix <- (distancematrix-488.0174)/300.3352
# -- Next we'll use the Poisson function to transform rescaled distances to 
#    probabilities
probmatrix <- exp(-0.62392 - 0.12618*rescaled_distancematrix)

# Load ASA microred shapefile for mapping
asa <- st_read("data/MR_ALTO_SELVA_ALGRE.kml") %>%
  st_transform(29188)

##-----------------------------------------------------------------------------
# 2. Run algorithm over 1000 iterations
##-----------------------------------------------------------------------------

set.seed(123)

# Get 1000 sets of solutions using 1000 different seeds
rand_seed <- sample(999999, 1000)
solutions <- list()

tic()
for (i in 1:1000){
  set.seed(rand_seed[i])
  solutions[[i]] <- OptimizeParticipationProbability(probmatrix, p = 20)
}
toc()

# Parse solutions and save
scores <- numeric()
supplysets <- list()
for (i in 1:1000){
  scores[i] <- solutions[[i]]$score
  supplysets[[i]] <- solutions[[i]]$supply
}

save(solutions, scores, supplysets, file = "pprob_solutions.rda")

##-----------------------------------------------------------------------------
# 3. Obtain best solution
##-----------------------------------------------------------------------------

load("output/pprob_solutions.rda")

length(unique(scores))  # 647
hist(scores)
which.max(scores)  # 14

solutions[[14]]  
sort(supplysets[[14]]) ## REPLACE WITH INDEX
#  [1]  1  2  7  9 15 18 19 20 22 26 28 29 31 37 41 45 47 57 61 63
supplysets_pprob   <- c(1, 2, 7, 9, 15, 18, 19, 20, 22, 26, 28, 29, 31, 37, 41,
                        45, 47, 57, 61, 63)

