library(sp)
library(spdep)
library(sf)
library(data.table)
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)
library(dplyr)
library(splines)

#Set your working directory
### setwd()

d1_full <- fread(file="monthly_deforestation_data.csv")
muns = st_read("municipality_shapefile.shp")

muns$mun <- substr(muns$CD_GEOCMU,1,6)
muns$UF <- substr(muns$CD_GEOCMU,1,2)
muns$UF <- as.numeric(muns$UF)
muns <- muns[muns$UF%in%c(11:17,21,51)&!(muns$mun%in%150475),]
muns <- muns[!is.na(muns$UF),]
muns2 <- subset(muns,select="mun")


# Create the spatial neighborhood matrix
temp <- poly2nb(muns2)
listmuns <- as.data.frame(cbind(1:807,muns2$mun2))
nb2INLA("mun.graph", temp)
mun.adj <- paste(getwd(),"/mun.graph",sep="")


# Define the cosine formula with a phase shift
cosine_formula <- function(x, phase_shift) {
  amplitude <- (1.193519 - 0.8383214) / 2
  vertical_shift <- 1.02
  frequency <- pi / 5
  
  y <- amplitude * cos(frequency * (x + phase_shift)) + vertical_shift
  return(y)
}

# Calculate the phase shift based on the desired point (x = 1, y = 0.98)
desired_x <- 1
desired_y <- 0.98
amplitude <- (1.193519 - 0.8383214) / 2
frequency <- pi / 5
vertical_shift <- 1.02
phase_shift <- (asin((desired_y - vertical_shift) / amplitude) / frequency - desired_x)*5.9
phase_shift
# Generate x values
x_values <- seq(1, 12, length.out = 100)  # Adjust the length.out for smoother plot

# Calculate corresponding y values using the formula with the calculated phase shift
y_values <- cosine_formula(x_values, phase_shift)

d1_full$V1 <- NULL
d1_full$month_control_cosine <- cosine_formula(d1_full$month, phase_shift)
muns_unique <- as.data.frame(unique(d1_full$mun_exp))
muns_unique$munID_correct <- 1:807
colnames(muns_unique)[1] <- "mun_exp"
d1_full <- merge(d1_full,muns_unique,by="mun_exp",all.x=T)
d1_full$mun.1 <- d1_full$munID_correct
d1_full$mun.2 <- d1_full$munID_correct

set.seed(324)

# 
# Model
m21 <- cases ~ 1 +
  scale(dflag1m_perc_area) +
  scale(area) +
  scale(perc_area_difference_deforestation_current_to_prev_year) +
  scale(forested_area_total/area) +
  scale(population) +
  scaled_imp_local_ratio +
  scale(conservation_area_km2/area) +
  scale(settlement_cases_perc) +
  scale(rural_cases_perc) +
  scale(indigenous_cases_perc) +
  scale(Unpaved) +
  scale(Paved) +
  month_control_cosine +
  as.factor(year) + ## Add year fixed effect
  scale(dflag1m_perc_area)*scale(forested_area_total/area) +
  f(inla.group(tmax_max), model="ar1") +
  f(inla.group(p_tot), model="ar1") +
  f(mun.1, model="bym", graph=mun.adj, scale.model=T, param=c(0.001,0.001)) +
  f(int_mun_time, model="iid") +
  f(uf_exp, model="iid", param=c(0.001,0.001)) +
  f(mun.2, model="iid", param=c(0.001,0.001)) +
  f(ONI,model="ar1")

m21_results <- inla(m21, family="zeroinflatedpoisson0", data=d1_full,
                    # control.predictor=list(compute=TRUE,link=NA),
                    control.compute=list(dic=TRUE, waic = TRUE, cpo = TRUE, config=TRUE),
                    # control.inla=list(tolerance=0.01),
                    verbose = T, safe = T)

m21_results_pois <- inla(m21, family="poisson", data=d1_full,
                         # control.predictor=list(compute=TRUE,link=NA),
                         control.compute=list(dic=TRUE, waic = TRUE, cpo = TRUE, config=TRUE),
                         # control.inla=list(tolerance=0.01),
                         verbose = T, safe = T)

m21_results_nb <- inla(m21, family="zeroinflatednbinomial0", data=d1_full,
                       # control.predictor=list(compute=TRUE,link=NA),
                       control.compute=list(dic=TRUE, waic = TRUE, cpo = TRUE, config=TRUE),
                       # control.inla=list(tolerance=0.01),
                       verbose = T, safe = T)

m21_summary <- exp(m21_results$summary.fixed)
m21_summary_pois <- exp(m21_results_pois$summary.fixed)
m21_results_nb <- exp(m21_results_nb$summary.fixed)

write.csv(m21_summary,file="global_model_output.csv")
