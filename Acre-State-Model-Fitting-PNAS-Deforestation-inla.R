library(rgdal)
library(sp)
library(rgeos)
library(spdep)
library(sf)
install.packages("maptools")
library(maptools)
library(rgeos)
library(raster)
install.packages("spacetime")
library(spacetime)
library(data.table)
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)
library(dplyr)
library(splines)

d1_full <- fread(file="monthly_deforestation_data.csv")
d1_full <- d1_full[d1_full$uf_exp==12]

muns = st_read("municipality_shapefile.shp")

muns$mun <- substr(muns$CD_GEOCMU,1,6)
muns$UF <- substr(muns$CD_GEOCMU,1,2)
muns$UF <- as.numeric(muns$UF)
muns <- muns[!is.na(muns$UF),]
muns2 <-  muns[muns$UF==12,]
muns2 <- subset(muns2,select="mun")

# Create the spatial neighborhood matrix
temp <- poly2nb(muns2,row.names = muns2)
listmuns <- as.data.frame(cbind(1:22,muns2$mun2))
nb2INLA("mun.graph", temp)
mun.adj <- paste(getwd(),"/mun.graph",sep="")
H <- inla.read.graph(filename="mun.graph")
image(inla.graph2matrix(H),xlab="",ylab="")

cosine_formula <- function(t) {
  y <- -28.94 * sin(2*pi*(t+207.86)/213.79) + 29.78
  return(y)
}

d1_full$V1 <- NULL
d1_full$month_control_cosine <- cosine_formula(d1_full$month)
muns_unique <- as.data.frame(unique(d1_full$mun_exp))
muns_unique$munID_correct <- 1:22
colnames(muns_unique)[1] <- "mun_exp"
d1_full <- merge(d1_full,muns_unique,by="mun_exp",all.x=T)
d1_full$mun.1 <- d1_full$munID_correct
d1_full$mun.2 <- d1_full$munID_correct

set.seed(324)
# 
# #Model 21 AC Only
m21_AC <- cases ~ 1 +
  scale(dflag1m_perc_area) +
  scale(area) +
  scale(perc_area_difference_deforestation_current_to_prev_year) +
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
  f(mun.2, model="iid", param=c(0.001,0.001)) +
  f(ONI,model="ar1")

m21_AC_results <- inla(m21_AC, family="zeroinflatedpoisson0", data=d1_full[d1_full$uf_exp==12],
                       # control.predictor=list(compute=TRUE,link=NA),
                       control.compute=list(dic=TRUE, waic = TRUE, cpo = TRUE, config=TRUE),
                       # control.inla=list(tolerance=0.01),
                       verbose = T, safe = T)

m21_AC_summary <- exp(m21_AC_results$summary.fixed)

write.csv(m21_AC_summary,file="AC_model_output.csv")
