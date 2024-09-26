install.packages("ggplot2")
library(ggplot2)
install.packages("ggjoy")
library(ggjoy)
library(dplyr)
library(ggpubr)
library(tidyr)
library(RColorBrewer)

data <- read.csv(file="monthly_deforestation_data.csv")

mapsFG <- data %>%
  group_by(mun_exp) %>%
  summarise(accumulatedDF = sum(annual_deforestation_km2),
            meanDF = mean(annual_deforestation_km2),
            accumulatedcases = sum(cases),
            meancases = mean(cases),
            meanindigenouscases = mean(indigenous_cases),
            accumulatedindigenouscases = sum(indigenous_cases),
            accumulatedgarimpocases = sum(garimpo_cases))

mapsC <- data %>%
  group_by(mun_exp) %>%
  summarise(maxdf = max(annual_deforestation_km2),
            annual_df = annual_deforestation_km2,
            year = year)

mapsC <- unique(mapsC)
mapsC$max_year <- 0
mapsC$max_year[mapsC$maxdf==mapsC$annual_df] <- 1

mapsC_final <- mapsC[mapsC$max_year==1,]
mapsC_final <- mapsC_final[mapsC_final$maxdf>0,]

write.csv(mapsC_final,"C:/Users/nicho/OneDrive/Dissertation/Chapter II/Current Analyses/Revised Figure 1/Data/mapCdata.csv")
write.csv(mapsFG,"C:/Users/nicho/OneDrive/Dissertation/Chapter II/Current Analyses/Revised Figure 1/Data/mapFGdata.csv")
