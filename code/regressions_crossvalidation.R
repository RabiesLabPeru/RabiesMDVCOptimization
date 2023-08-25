library(tidyverse)
library(gridExtra)
library(AICcmodavg)
library(h2o)
library(boot)
library(lme4)

setwd("~/Optimized tents paper/")

##-----------------------------------------------------------------------------
# 0. User-defined functions
##-----------------------------------------------------------------------------

# Bin data by 30 m
BinData <- function(df){
  
  # Make 30m bin ref table
  BinWidth <- 30 # set bin width
  MaxSeq <- ceiling(max(df$distance)/BinWidth) # Number of bins
  bins <- data.frame(Bin = seq(1:MaxSeq),  # Make ref table
                     StartDist= BinWidth*(seq(1:MaxSeq) - 1),
                     EndDist = BinWidth*seq(1:MaxSeq)) #ref table
  
  # Assign a bin number based on ref table
  df$Bin <- NA
  for(i in 1:nrow(bins)){
    row = bins[i,]
    df <- df %>%
      dplyr::mutate(Bin = ifelse(distance >= row$StartDist[[1]] & distance < row$EndDist[[1]],
                                 row$Bin[[1]], Bin))
  }
  
  # Find average distance for each bin and year 
  vac_binavg <- df %>%
    group_by(Bin, year) %>%
    summarise(mean_distance = mean(distance))
  
  # Find the number of houses in each bin by vaccine status and year and merge
  # other bin variables
  df_out <- data.frame(plyr::count(df%>%select(vac_status, Bin, year)))%>%
    tidyr::pivot_wider(names_from = vac_status, values_from = freq)%>%
    rename(num_not_vac = '0',
           num_vac ='1')%>%
    mutate(num_vac = tidyr::replace_na(num_vac, 0),
           num_not_vac = tidyr::replace_na(num_not_vac, 0),
           num_of_houses = num_not_vac+num_vac,
           VaccFreq = num_vac/num_of_houses) %>%
    arrange(year, Bin) %>%
    left_join(., vac_binavg, by = c("year", "Bin")) %>%
    left_join(., bins)
  
  return(df_out)
}

BinPredictions <- function(pred_vec, df){
  
  # Make 30m bin ref table
  BinWidth <- 30 # set bin width
  MaxSeq <- ceiling(max(df$distance)/BinWidth) # Number of bins
  bins <- data.frame(Bin = seq(1:MaxSeq),  # Make ref table
                     StartDist= BinWidth*(seq(1:MaxSeq) - 1),
                     EndDist = BinWidth*seq(1:MaxSeq)) #ref table
  
  # Assign a bin number based on ref table
  df$Bin <- NA
  for(i in 1:nrow(bins)){
    row = bins[i,]
    df <- df %>%
      dplyr::mutate(Bin = ifelse(distance >= row$StartDist[[1]] & distance < row$EndDist[[1]],
                                 row$Bin[[1]], Bin))
  }
  
  # Bin predictions by summing probabilities for each distance bin and year
  df$year_bin <- paste0(df$year, "_", df$Bin)
  year_bin <- sort(unique(df$year_bin))
  pred_bin <- numeric()
  j <- 1
  for(i in year_bin){
    index <- df$year_bin == i
    pred_bin[j] <- sum(pred_vec[index])
    j <- j+1
  }
  
  return(pred_bin)
}

##-----------------------------------------------------------------------------
# 1. Load and clean data
##-----------------------------------------------------------------------------

set.seed(123)

# Load raw vaccination data
vac_individual <- readRDS("data/cleanedhouses_3.23.23.rds")
vac_binned <- BinData(vac_individual)

# Load regression models
load(file = "models_standard.rda")

##-----------------------------------------------------------------------------
# 2. Find prediction error for Poisson regression model using 5-fold CV
##-----------------------------------------------------------------------------

# K-Fold Cross Validation

# Divide data into K = 5 groups
N <- nrow(vac_individual)
K <- 5
set.seed(123)
s <- sample(rep(1:K, ceiling(N/K))[1:N], size=N)

# Create empty vectors to store CV results
prederror_poisson <- numeric()
prederror_negbin <- numeric()
prederror_binom <- numeric()
prederror_binquad <- numeric()
prederror_poissonme <- numeric()
prederror_negbinme <- numeric()
prederror_binomme <- numeric()
prederror_binquadme <- numeric()
rmse_poisson <- numeric()
rmse_negbin <- numeric()
rmse_binom <- numeric()
rmse_binquad <- numeric()
rmse_poissonme <- numeric()
rmse_negbinme <- numeric()
rmse_binomme <- numeric()
rmse_binquadme <- numeric()

for(i in 1:K){
  
  # Divide into train and test sets
  train <- filter(vac_individual, s != i)
  test <- filter(vac_individual, s == i)
  obs_outputs[1:length(s[s==i]) + offset] <- test$vac_status
  
  # Bin data for Poisson and negative binomial models
  train_binned <- BinData(train)
  test_binned <- BinData(test)
  
  # Poisson train/test
  poisreg <- glm(num_vac ~ mean_distance + offset(log(num_of_houses)),
               family = poisson(link=log), data=train_binned)
  predicted_vac <- predict(poisreg, newdata = test_binned, type = "response")
  prederror_poisson[i] <- sum(test_binned$num_vac-predicted_vac)/sum(test_binned$num_vac)*100 
  rmse_poisson[i] <- sqrt(sum((test_binned$num_vac-predicted_vac)^2)/sum(test_binned$num_vac))  

  # Negative binomial train/test
  negbinreg<-MASS::glm.nb(num_vac~mean_distance+offset(log(num_of_houses)), 
                          data=train_binned)
  predicted_vac <- predict(negbinreg, newdata = test_binned, type = "response")
  prederror_negbin[i] <- sum(test_binned$num_vac-predicted_vac)/sum(test_binned$num_vac)*100 
  rmse_negbin[i] <- sqrt(sum((test_binned$num_vac-predicted_vac)^2)/sum(test_binned$num_vac))
  
  # Binomial with linear terms train/test
  # **Note we will bin predictions to make output comparable to those for the 
  #   Poisson and negative binomial models
  binomreg <-glm(formula = vac_status ~ distance, family = "binomial", 
                 data = train)
  predicted_ind <- predict(binomreg, newdata = test, type = "response")
  predicted_vac <- BinPredictions(predicted_ind, test)
  prederror_binom[i] <- sum(test_binned$num_vac-predicted_vac)/sum(test_binned$num_vac)*100
  rmse_binom[i] <- sqrt(sum((test_binned$num_vac-predicted_vac)^2)/sum(test_binned$num_vac)) 
  
  # Binomial with linear + quadratic terms train/test
  # **As above, we'll bin predictions 
  binquadreg <- glm(formula = vac_status ~ distance + I(distance^2), 
                    family = "binomial", data = train)
  predicted_ind <- predict(binquadreg, newdata = test, type = "response")
  predicted_vac <- BinPredictions(predicted_ind, test)
  prederror_binquad[i] <- sum(test_binned$num_vac-predicted_vac)/sum(test_binned$num_vac)*100  
  rmse_binquad[i] <- sqrt(sum((test_binned$num_vac-predicted_vac)^2)/sum(test_binned$num_vac))
  
  # Poisson mixed effects train/test
  # **We need to rescale mean_distance for model to converge
  poisme <- glmer(num_vac ~ scale(mean_distance) + offset(log(num_of_houses)) + (1|year),
                  data = train_binned, family = poisson)
  predicted_vac <- predict(poisme, newdata = test_binned, type = "response")
  prederror_poissonme[i] <- sum(test_binned$num_vac-predicted_vac)/sum(test_binned$num_vac)*100 
  rmse_poissonme[i] <- sqrt(sum((test_binned$num_vac-predicted_vac)^2)/sum(test_binned$num_vac))  
  
  # Negative binomial mixed effects train/test
  negbinme <- glmer.nb(num_vac ~ scale(mean_distance) + offset(log(num_of_houses)) + (1|year),
                       data = train_binned) 
  predicted_vac <- predict(negbinme, newdata = test_binned, type = "response")
  prederror_negbinme[i] <- sum(test_binned$num_vac-predicted_vac)/sum(test_binned$num_vac)*100 
  rmse_negbinme[i] <- sqrt(sum((test_binned$num_vac-predicted_vac)^2)/sum(test_binned$num_vac))  
  
  # Binomial - linear terms mixed effects train/test
  binomme <- glmer(formula = vac_status ~ log(distance) + (1|year), family = "binomial", 
                   data = train)
  predicted_ind <- predict(binomme, newdata = test, type = "response")
  predicted_vac <- BinPredictions(predicted_ind, test)
  prederror_binomme[i] <- sum(test_binned$num_vac-predicted_vac)/sum(test_binned$num_vac)*100
  rmse_binomme[i] <- sqrt(sum((test_binned$num_vac-predicted_vac)^2)/sum(test_binned$num_vac)) 
  
  # Binomial - linear + quadratic terms mixed effects train/test
  binquadme <- glmer(formula = vac_status ~ log(distance) + I(log(distance)^2)+ (1|year), 
                     family = "binomial", data = vac_individual)
  predicted_ind <- predict(binquadme, newdata = test, type = "response")
  predicted_vac <- BinPredictions(predicted_ind, test)
  prederror_binquadme[i] <- sum(test_binned$num_vac-predicted_vac)/sum(test_binned$num_vac)*100
  rmse_binquadme[i] <- sqrt(sum((test_binned$num_vac-predicted_vac)^2)/sum(test_binned$num_vac)) 
  
}

# Save cv results in a data frame
models <- c("Poisson", "Neg binomial", "Binomial - linear", 
            "Binomial - quadratic", "Poisson ME", "Neg binomial ME", 
            "Binomial - linear ME", "Binomial - quadratic ME")
results <- data.frame(model = models, prederror = NA, RMSE = NA)

results$prederror <- c(mean(prederror_poisson), mean(prederror_negbin), 
                       mean(prederror_binom), mean(prederror_binquad),
                       mean(prederror_poissonme), mean(prederror_negbinme),
                       mean(prederror_binomme), mean(prederror_binquadme))

results$RMSE <- c(mean(rmse_poisson), mean(rmse_negbin), 
                       mean(rmse_binom), mean(rmse_binquad),
                       mean(rmse_poissonme), mean(rmse_negbinme),
                       mean(rmse_binomme), mean(rmse_binquadme))

results

##-----------------------------------------------------------------------------
# 3. Find prediction error, RMSE and chi-square results for NEGATIVE BINOMIAL
#    regression model
##-----------------------------------------------------------------------------

cv.err <- cv.glm(vac_binned, negbinreg, K = 10)
summary(cv.err)
cv.err$delta

# Prediction error of the model = 6.384055

##-----------------------------------------------------------------------------
# 4. Find prediction error, RMSE and chi-square results for BINOMIAL-LINEAR
#    regression model
##-----------------------------------------------------------------------------

cv.err <- cv.glm(vac_raw, binomreg, K = 10)
summary(cv.err)
cv.err$delta  

##-----------------------------------------------------------------------------
# 5. Find prediction error, RMSE and chi-square results for BINOMIAL-QUADRATIC
#    regression model
##-----------------------------------------------------------------------------

# Fit model using train data 
binquadreg <- glm(formula = vac_status ~ distance + I(distance^2), 
                  family = "binomial", data = train_individual)
summary(binquadreg)
# Intercept (beta0)  = 1.311; distance (beta1) = -3.677e-03; 
# distance^2 (beta2) = 2.306e-06

# Use predict to test data and obtain prediction error
# **Note we are binning predictions by mean distance to make output comparable
# to those for the Poisson and negative binomial models
predicted_vac<-plogis((1.311 - 3.677e-03*test_binned$mean_distance+2.306e-06*(test_binned$mean_distance)^2))*test_binned$num_of_houses
sum(test_binned$num_vac-predicted_vac)/sum(test_binned$num_vac)*100  # -5.697%

# RMSE
sqrt(sum((test_binned$num_vac-predicted_vac)^2)/sum(test_binned$num_vac))  # 3.149063

# Chi-squared
chisq.test(cbind(test_binned$num_vac, predicted_vac))

