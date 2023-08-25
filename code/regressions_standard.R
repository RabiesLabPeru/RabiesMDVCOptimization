library(tidyverse)
library(gridExtra)
library(AICcmodavg)
library(h2o)
library(caret)

setwd("~/Optimization")

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
vac_raw <- readRDS("data_cleaned/cleanedhouses_3.23.23.rds")
vac_binned <- BinData(vac_raw)
max_d <- max(vac_binned$mean_distance)  # 1068.38

# Check that house totals is the same for raw and binned data
nrow(vac_raw)  # 1855
sum(vac_binned$num_of_houses)  # 1855

# Check Poisson assumption -- does mean ~ var?
plot(vac_binned$num_of_houses, vac_binned$num_vac, xlab = "Number of houses",
     ylab = "Number vaccinated")
plot(vac_binned$mean_distance, vac_binned$VaccFreq, xlab = "Mean distance",
     ylab = "Proportion vaccinated")

##-----------------------------------------------------------------------------
# 2. Make scaled scatter plots of vaccination coverage vs. mean walking dist.
##-----------------------------------------------------------------------------

# Pick color scale (qualitative scale from https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=4)
cols <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")
#cols <- c("#bf4141", "#377eb8", "#4daf4a", "#984ea3")  # red, blue, green, purple

# Make plot of ALL data
vac_binned %>%
  ggplot(., aes(x=mean_distance, y=VaccFreq, size=num_of_houses, color=as.factor(year))) +
  scale_color_manual(values = cols) +
  geom_point()+
  theme_classic()+
  ggtitle("ALL YEARS")+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
  xlab("Bin mean walking distance (m)") + ylab("Vaccination coverage (%)")
p1 <- vac_binned %>%
  filter(year == 2016) %>%
  ggplot(., aes(x=mean_distance, y=VaccFreq, size=num_of_houses)) +
  geom_point(color=cols[1])+
  theme_classic()+
  ggtitle("2016")+
  xlim(0, max_d)+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
  xlab("Bin mean walking distance (m)") + ylab("Vaccination coverage (%)")
p2 <- vac_binned %>%
  filter(year == 2017) %>%
  ggplot(., aes(x=mean_distance, y=VaccFreq, size=num_of_houses)) +
  geom_point(color=cols[2])+
  theme_classic()+
  ggtitle("2017")+
  xlim(0, max_d)+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
  xlab("Bin mean walking distance (m)") + ylab("Vaccination coverage (%)")
p3 <- vac_binned %>%
  filter(year == 2018) %>%
  ggplot(., aes(x=mean_distance, y=VaccFreq, size=num_of_houses)) +
  geom_point(color=cols[3])+
  theme_classic()+
  ggtitle("2018")+
  xlim(0, max_d)+
  ylim(0, 1)+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
  xlab("Bin mean walking distance (m)") + ylab("Vaccination coverage (%)")
p4 <- vac_binned %>%
  filter(year == 2019) %>%
  ggplot(., aes(x=mean_distance, y=VaccFreq, size=num_of_houses)) +
  geom_point(color=cols[4])+
  theme_classic()+
  ggtitle("2019")+
  xlim(0, max_d)+
  ylim(0, 1)+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
  xlab("Bin mean walking distance (m)") + ylab("Vaccination coverage (%)")
grid.arrange(p1, p2, p3, p4, ncol = 2)

##-----------------------------------------------------------------------------
# 4. Poisson regression
##-----------------------------------------------------------------------------

# Fit model
poisreg<-glm(num_vac~mean_distance+offset(log(num_of_houses)),
             family=poisson(link=log), data=vac_binned)
summary(poisreg)

# Calculate regression curve
coef(poisreg)
#   (Intercept) mean_distance 
# -0.2890246742 -0.0006483417 
summary(vac_binned$mean_distance)  # range from 13-1068 m
xvals <- seq(13, 1068, by = 1)
yvals_pois <- exp(-0.2890246742 - 0.0006483417*xvals)
xy_pois <- bind_rows(x = xvals, y = yvals_pois, type = "Poisson")

# Plot scatter plot with regression curve
#pdf("output_pdf/Poissonregression.pdf", width = 6, height = 4)
ggplot()+
  theme_classic()+
  geom_point(data=vac_binned, aes(x=mean_distance, y=VaccFreq, size=num_of_houses,
                                    color=as.factor(year)))+
  scale_color_manual(values = cols) +
  geom_line(data = xy_pois,
            aes(x = x, y = y, linetype = type)) +
  xlab("Bin mean walking distance (m)") + ylab("Vaccination coverage (%)")
#dev.off()

# View model residuals 
p_res <- resid(poisreg)
plot(fitted(poisreg), p_res, col='steelblue', pch=16,
     xlab='Predicted Offers', ylab='Standardized Residuals', main='Poisson')
abline(0,0)

##-----------------------------------------------------------------------------
# 5. Negative binomial
##-----------------------------------------------------------------------------

# Fit model
negbinreg<-MASS::glm.nb(num_vac~mean_distance+offset(log(num_of_houses)), 
                        data=vac_binned)  # xx: warning iteration reached
summary(negbinreg)
# Theta estimated to be very large (762980)  
# This implies that the variance is equal to the mean and the model is 
# equivalent to Poisson

# Calculate regression line
coef(negbinreg)
#   (Intercept) mean_distance 
# -0.2890296865 -0.0006483349  
yvals_negbin <- exp(-0.2890296865 - 0.0006483349*xvals)
xy_negbin <- bind_rows(x = xvals, y = yvals_negbin, type = "Negative binomial")

# Use likelihood ratio test to compare Poisson vs. negative binomial model
# Source: https://stats.oarc.ucla.edu/r/dae/negative-binomial-regression/
pchisq(2 * (logLik(poisreg) - logLik(negbinreg)), df = 1, lower.tail = FALSE)
# 'log Lik.' 0.9458546 (df=2)  # More clear in the code below which is stat
# and which is the pvalue

# p-value 
# Source: https://stats.stackexchange.com/questions/127505/compare-poisson-and-negative-binomial-regression-with-lr-test
stat <- as.numeric(2 * (logLik(poisreg) - logLik(negbinreg)))
stat  # 0.004612225
pchisq(stat, df = 1, lower.tail = FALSE)
# pvalue = 0.9458546

##-----------------------------------------------------------------------------
# 6. Binomial regression with linear term
##-----------------------------------------------------------------------------

# Fit model
binomreg <-glm(formula = vac_status ~ distance, family = "binomial", 
                 data = vac_raw)
summary(binomreg)

# Calculate regression curve
coef(binomreg)
# (Intercept)     distance 
# 0.961606574 -0.001537142 
summary(vac_raw$distance)  # range from 2-1079 m
xvals <- seq(2, 1079, by = 1)
yvals_binom <- 1/(1+exp(-0.961606574 + 0.001537142*xvals))
xy_binom <- bind_rows(x = xvals, y = yvals_binom, type = "Binomial - linear")

# Plot scatter plot with regression curve
ggplot()+
  theme_classic()+
  geom_point(data=vac_binned, aes(x=mean_distance, y=VaccFreq, size=num_of_houses,
                                  color=as.factor(year)))+
  scale_color_manual(values = cols) +
  geom_line(data = xy_binom,
            aes(x = x, y = y, linetype = type)) +
  xlab("Bin mean walking distance (m)") + ylab("Vaccination coverage (%)")


##-----------------------------------------------------------------------------
# 7. Binomial regression with QUADRATIC term
##-----------------------------------------------------------------------------

# Fit model - note commented out version is equivalent to the below
#model_binom_quad <- glm(formula = vac_status ~ poly(distance, 2, raw = T), 
#                        family = "binomial", data = vac_raw)
binquadreg <- glm(formula = vac_status ~ distance + I(distance^2), 
                        family = "binomial", data = vac_raw)
summary(binquadreg)

# Calculate regression curve
coef(binquadreg)
# (Intercept)      distance I(distance^2) 
# 1.309577e+00 -3.564847e-03  2.214972e-06 
yvals_binquadreg <- 1/(1+exp(-1.309577 + 3.564847e-03*xvals - 2.214972e-06*xvals^2))
xy_binquadreg <- bind_rows(x = xvals, y = yvals_binquadreg, type = "Binomial - quadratic")

# Plot scatter plot with regression curve
ggplot()+
  theme_classic()+
  geom_point(data=vac_binned, aes(x=mean_distance, y=VaccFreq, size=num_of_houses,
                                       color=as.factor(year)))+
  scale_color_manual(values = cols) +
  geom_line(data = xy_binquadreg,
            aes(x = x, y = y, linetype = type)) +
  xlab("Bin mean walking distance (m)") + ylab("Vaccination coverage (%)")

# Out of curiosity let's see what the binomial-quadratic function looks like at
# greater distances
xvals2 <- 1:2000
yvals_binquadreg2 <- 1/(1+exp(-1.309577 + 3.564847e-03*xvals2 - 2.214972e-06*xvals2^2))
xy_binquadreg2 <- bind_rows(x = xvals2, y = yvals_binquadreg2, type = "Binomial - quadratic")

ggplot()+
  theme_classic()+
  geom_point(data=vac_binned, aes(x=mean_distance, y=VaccFreq, size=num_of_houses,
                                  color=as.factor(year)))+
  scale_color_manual(values = cols) +
  geom_line(data = xy_binquadreg2,
            aes(x = x, y = y, linetype = type)) +
  xlim(c(0, 2000)) +
  xlab("Bin mean walking distance (m)") + ylab("Vaccination coverage (%)")


# Save output from models
setwd("~/Optimized tents paper/")
save(poisreg, negbinreg, binomreg, binquadreg, file = "models_standard.rda")

##-----------------------------------------------------------------------------
# 8. Plot regression lines for all models on same graph
##-----------------------------------------------------------------------------

# Combine x/y vals for the 4 regression models
reg_curves <- bind_rows(xy_pois, xy_negbin, xy_binom, xy_binquadreg)

# Convert curve type to ordered factor
reg_curves$type <- factor(reg_curves$type, ordered = T, 
                          levels = c("Poisson", "Negative binomial", 
                                     "Binomial - linear", "Binomial - quadratic"))

saveRDS(reg_curves, "regressioncurves_standard.rds")
