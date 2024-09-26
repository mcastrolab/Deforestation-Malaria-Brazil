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
library(dlnm)

d1_full <- fread(file="/n/home04/nicholasarisco/PNAS -- Deforestation Paper/Updating Data 2003 to 2022/final_model_data_2003to2022_Mar1.csv")
muns = st_read("/n/home04/nicholasarisco/PNAS -- Deforestation Paper/Model data -- Ch2/BRMUE250GC_SIR.shp")

d1_full %>%
  group_by(year) %>%
  summarise(cases = sum(cases))

muns$mun <- substr(muns$CD_GEOCMU,1,6)
muns$UF <- substr(muns$CD_GEOCMU,1,2)
muns$UF <- as.numeric(muns$UF)
muns <- muns[!is.na(muns$UF),]
muns2 <- subset(muns,select="mun")

# Create the spatial neighborhood matrix
temp <- poly2nb(muns2,row.names = muns2)
listmuns <- as.data.frame(cbind(1:807,muns2$mun2))
nb2INLA("mun.graph", temp)
mun.adj <- paste(getwd(),"/mun.graph",sep="")
H <- inla.read.graph(filename="mun.graph")
image(inla.graph2matrix(H),xlab="",ylab="")

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

dfknots <- as.vector(c(quantile(d1_full$dflag1m_perc_area,0.10),
                                      quantile(d1_full$dflag1m_perc_area,0.75),
                                      quantile(d1_full$dflag1m_perc_area,0.95)))

df_CB <- crossbasis(d1_full$dflag1m_perc_area, as.factor(d1_full$mun_exp), argvar=list(fun="bs", knots=dfknots), lag=6,
                                 arglag=list(fun="bs", df=4))#choose two sets of basis functions (variable + lag dimension), i.e. set cross basis for exposure and lag dimension 

# 
# Model
m21 <- cases ~ 1 +
  bs(dflag1m_perc_area, knots=dfknots, df=3) +
  scale(area) +
  scale(perc_area_difference_deforestation_current_to_prev_year) +
  scale(forested_area_total/area) +
  # scale(cons_area_perc) +
  scale(population) +
  scaled_imp_local_ratio +
  scale(conservation_area_km2/area) +
  # scale(homologized_area_indigenous_km2/area) +
  # scale(nonhomologized_area_indigenous_km2/area) +
  scale(garimpo_cases_perc) +
  scale(settlement_cases_perc) +
  scale(rural_cases_perc) +
  scale(indigenous_cases_perc) +
  scale(Unpaved) +
  scale(Paved) +
  month_control_cosine +
  # as.factor(uf_exp) +
  as.factor(year) + ## Add year fixed effect
  scale(dflag1m_perc_area)*scale(forested_area_total/area) +
  # as.factor(globalfund) +
  # as.factor(elimination_plan) +
  # as.factor(Bolsonaro) +
  # as.factor(bednets2007) +
  # as.factor(ACT) +
  # as.factor(national_policy) + ## Redoing policy information
  # f(year_order,model="ar1", hyper=list(prec=list(prior="loggamma",param=c(1,5e-5)), rho=list(prior="normal",param=c(0,0.4)))) +
  f(inla.group(tmax_max), model="ar1") +
  f(inla.group(p_tot), model="ar1") +
  f(mun.1, model="bym", graph=mun.adj, scale.model=T, param=c(0.001,0.001)) +
  # f(month, model="ar1", hyper=list(prec=list(prior="loggamma",param=c(1,5e-5)), rho=list(prior="normal",param=c(0,0.4)))) +
  f(int_mun_time, model="iid") +
  f(uf_exp, model="iid", param=c(0.001,0.001)) +
  f(mun.2, model="iid", param=c(0.001,0.001)) +
  f(ONI,model="ar1")

m21_results_nonlinear <- inla(m21, family="zeroinflatedpoisson0", data=d1_full,
                    # control.predictor=list(compute=TRUE,link=NA),
                    control.compute=list(dic=TRUE, waic = TRUE, cpo = TRUE, config=TRUE),
                    # control.inla=list(tolerance=0.01),
                    verbose = T, safe = T)



pred.precipitation1week <- as.data.frame(with(pred.precipitation1, t(rbind(allRRfit,allRRlow,allRRhigh))))
pred.precipitation1week$prec <- rownames(pred.precipitation1week)
pred.precipitation1week$lag <- "1 weeks"


m21_results_nonlinear$summary.fixed


# m21_summary <- exp(m21_results_nonlinear$summary.fixed)

# write.csv(m21_summary,file="/n/home04/nicholasarisco/PNAS -- Deforestation Paper/GitHub/Primary Analysis/INLA Models/Global Model/global_model_output.csv")
# 
# m21_results$waic$waic
# sum(m21_results$cpo$failure)
# 
# summary(m21_results)
# m21_results$mlik
