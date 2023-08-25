library(tidyverse)

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

# Load regression curves for standard and mixed effects models
curves_standard <- readRDS("data/regressioncurves_standard.rds")
curves_mixed <- readRDS("data/regressioncurves_mixedeffects.rds")

##-----------------------------------------------------------------------------
# 2. Plot regression curves for STANDARD models
##-----------------------------------------------------------------------------

# Plot scatter plot with regression lines
cols <- c("#db91bd", "#82AFD3", "#DC9147", "#828F71")

pdf("figures/R_output/figure2_regressionSTANDARD_4.2x5.5.pdf", height = 4.2, width = 5.5)
ggplot()+
  theme_classic()+
  geom_point(data=vac_binned, aes(x=mean_distance, y=VaccFreq*100, size=num_of_houses,
                                  color=as.factor(year)))+
  scale_color_manual(values = cols) +
  geom_line(data = curves_standard,
            aes(x = x, y = y*100, linetype = type)) +
  scale_linetype_manual(values = c("solid", "solid", "dashed", "dotted")) + 
  xlab("Bin mean walking distance, m") + ylab("Vaccination coverage, %")
dev.off()

##-----------------------------------------------------------------------------
# 2. Plot regression curves for MIXED EFFECTS models
##-----------------------------------------------------------------------------

pdf("figures/R_output/figure2_regressionMIXEDEFFECTS_4.2x5.5.pdf", height = 4.2, width = 5.5)
ggplot()+
  theme_classic()+
  geom_point(data=vac_binned, aes(x=mean_distance, y=VaccFreq*100, size=num_of_houses,
                                  color=as.factor(year)))+
  scale_color_manual(values = cols) +
  geom_line(data = curves_mixed,
            aes(x = x, y = y*100, linetype = type)) +
  scale_linetype_manual(values = c("solid", "solid", "dashed", "dotted")) + 
  xlab("Bin mean walking distance, m") + ylab("Vaccination coverage, %")
dev.off()
