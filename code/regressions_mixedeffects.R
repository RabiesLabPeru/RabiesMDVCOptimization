library(tidyverse)
library(gridExtra)
library(AICcmodavg)
library(h2o)
library(caret)
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

##-----------------------------------------------------------------------------
# 1. Load and clean data
##-----------------------------------------------------------------------------

# Bin raw data
vac_individual <- readRDS("data/cleanedhouses_3.23.23.rds")
vac_binned <- BinData(vac_individual)

##-----------------------------------------------------------------------------
# 2. Poisson mixed effects regression
##-----------------------------------------------------------------------------

# Fit model - to avoid error, we need to rescale the mean_distance variable
# Note the scale function does the following transformation: 
# x_rescaled = (x-mean(x))/sd(x)
mean(vac_binned$mean_distance)  # 488.0174
sd(vac_binned$mean_distance)  # 300.3352
poisme <- glmer(num_vac ~ scale(mean_distance) + offset(log(num_of_houses)) + (1|year),
                data = vac_binned, family = poisson)
summary(poisme)
# Fixed effects:
#                      Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -0.62392    0.07297  -8.551   <2e-16 ***
# scale(mean_distance) -0.12618    0.05035  -2.506   0.0122 *  

coef(poisme)
# $year
#      (Intercept) scale(mean_distance)
# 2016  -0.4909919           -0.1261808
# 2017  -0.6244793           -0.1261808
# 2018  -0.5991846           -0.1261808
# 2019  -0.7757199           -0.1261808

# Create regression curve - note we'll have to rescale the distance (xvals)
xvals <- seq(13, 1068, by = 1)
y_shift <- mean(vac_binned$mean_distance)  # 488.0174 
m_shift <- sd(vac_binned$mean_distance)  # 300.3352
xvals_rescaled <- (xvals-y_shift)/m_shift
yvals_pois <- exp(-0.62392 - 0.12618*xvals_rescaled)
xy_pois <- bind_rows(x = xvals, y = yvals_pois, type = "Poisson")

##-----------------------------------------------------------------------------
# 5. Negative binomial
##-----------------------------------------------------------------------------

# Fit model
negbinme <- glmer.nb(num_vac ~ scale(mean_distance) + offset(log(num_of_houses)) + (1|year),
                     data = vac_binned) 
summary(negbinme)
# Fixed effects:
#                      Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -0.62393    0.07308  -8.538   <2e-16 ***
# scale(mean_distance) -0.12618    0.05038  -2.505   0.0123 *  

coef(negbinme)
# $year
#      (Intercept) scale(mean_distance)
# 2016  -0.4909941           -0.1261835
# 2017  -0.6244792           -0.1261835
# 2018  -0.5991892           -0.1261835
# 2019  -0.7757116           -0.1261835

# Create regression curve
yvals_negbin <- exp(-0.62393 - 0.12618 *xvals_rescaled)
xy_negbin <- bind_rows(x = xvals, y = yvals_negbin, type = "Negative binomial")

# Use likelihood ratio test to compare Poisson vs. negative binomial model
# Source: https://stats.oarc.ucla.edu/r/dae/negative-binomial-regression/
pchisq(2 * (logLik(poisme) - logLik(negbinme)), df = 1, lower.tail = FALSE)
# 'log Lik.' 0.9657839 (df=3) 
2 * (logLik(poisme) - logLik(negbinme))  # test statistic

##-----------------------------------------------------------------------------
# 6. Binomial regression with linear term
##-----------------------------------------------------------------------------

# Fit model
binomme <- glmer(formula = vac_status ~ log(distance) + (1|year), 
                 family = "binomial", data = vac_individual)
summary(binomme)
# Fixed effects:
#               Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    2.30334    0.49910   4.615 3.93e-06 ***
# log(distance) -0.35225    0.08268  -4.261 2.04e-05 ***

coef(binomme)
# $year
#      (Intercept) log(distance)
# 2016    2.653014    -0.3522544
# 2017    2.284441    -0.3522544
# 2018    2.328484    -0.3522544
# 2019    1.944850    -0.3522544

# Create regression curve
xvals_log <- log(xvals)  # Because model was fit using log-transformed predictor variables
yvals_binom <- 1/(1+exp(-2.30334 + 0.35225*xvals_log))
xy_binom <- bind_rows(x = xvals, y = yvals_binom, type = "Binomial - linear")

##-----------------------------------------------------------------------------
# 7. Binomial regression with QUADRATIC term
##-----------------------------------------------------------------------------

# Fit model - note distance needed to be log-transformed for the model to run
binquadme <- glmer(formula = vac_status ~ log(distance) + I(log(distance)^2)+ (1|year), 
                   family = "binomial", data = vac_individual)
summary(binquadme)
# Fixed effects:
#                    Estimate Std. Error z value Pr(>|z|)
# (Intercept)         3.04995    1.98984   1.533    0.125
# log(distance)      -0.63483    0.73108  -0.868    0.385
# I(log(distance)^2)  0.02615    0.06704   0.390    0.697
# Create regression curve 

coef(binquadme)
# $year
# (Intercept) log(distance) I(log(distance)^2)
# 2016    3.405377    -0.6348331         0.02614824
# 2017    3.027436    -0.6348331         0.02614824
# 2018    3.078289    -0.6348331         0.02614824
# 2019    2.686062    -0.6348331         0.02614824

# Create regression curve
yvals_binom <- 1/(1+exp(-3.04995 + 0.63483*xvals_log - 0.02615*xvals_log^2))
xy_binquad <- bind_rows(x = xvals, y = yvals_binom, type = "Binomial - quadratic")

##-----------------------------------------------------------------------------
# 8. Save regression curves 
##-----------------------------------------------------------------------------

# Combine x/y vals for the 4 regression models
reg_curves <- bind_rows(xy_pois, xy_negbin, xy_binom, xy_binquad)

# Convert curve type to ordered factor
reg_curves$type <- factor(reg_curves$type, ordered = T, 
                          levels = c("Poisson", "Negative binomial", 
                                     "Binomial - linear", "Binomial - quadratic"))

saveRDS(reg_curves, "data/regressioncurves_mixedeffects.rds")

