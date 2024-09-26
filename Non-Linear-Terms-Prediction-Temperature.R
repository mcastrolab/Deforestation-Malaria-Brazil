library(ggplot2)

## Load Model ##
# m21_results <- readRDS("/n/holyscratch01/mcastro_lab/Users/nicholasarisco/PNAS -- Deforestation Paper/Sensitivity Analysis Placeholder/M21 Spec/m21_model.rds")

# Create an empty data frame to store the results
median_table <- data.frame(Index = numeric(25), Value = numeric(25))
LCI_table <- data.frame(Index = numeric(25), Value = numeric(25))
UCI_table <- data.frame(Index = numeric(25), Value = numeric(25))

# Loop through the values of i
for (i in 1:25) {
  # Replace the number after "index." with i
  index_name <- paste("index.", i, sep = "")
  
  # Extract the median value for the specified index
  median_value <- median(as.data.frame(m21_results$marginals.random$`inla.group(tmax_max)`[[index_name]])$x)
  
  # Store the results in the data frame
  median_table[i, ] <- c(i, exp(median_value))
}

# Print or use the results_table as needed
plot(median_table)

# Loop through the values of i
for (i in 1:25) {
  # Replace the number after "index." with i
  index_name <- paste("index.", i, sep = "")
  
  # Extract the median value for the specified index
  LCI_value <- quantile(as.data.frame(m21_results$marginals.random$`inla.group(tmax_max)`[[index_name]])$x,0.05)
  
  # Store the results in the data frame
  LCI_table[i, ] <- c(i, exp(LCI_value))
}

# Print or use the results_table as needed
plot(LCI_table)


# Loop through the values of i
for (i in 1:25) {
  # Replace the number after "index." with i
  index_name <- paste("index.", i, sep = "")
  
  # Extract the median value for the specified index
  UCI_value <- quantile(as.data.frame(m21_results$marginals.random$`inla.group(tmax_max)`[[index_name]])$x,0.95)
  
  # Store the results in the data frame
  UCI_table[i, ] <- c(i, exp(UCI_value))
}

# Print or use the results_table as needed
plot(UCI_table)

colnames(median_table) <- c("id","median")
colnames(LCI_table) <- c("id","LCI")
colnames(UCI_table) <- c("id","UCI")

tmax_CrI <- merge(LCI_table,UCI_table,by="id")
tmax_results <- merge(median_table,tmax_CrI,by="id")


p <- ggplot(tmax_results) + 
  geom_line(aes(y=median, x=id, colour = "sin")) +
  geom_ribbon(aes(ymin=LCI, ymax=UCI, x=id, fill = "band"), alpha = 0.3) +
  scale_colour_manual("",values="blue") +
  scale_fill_manual("",values="grey12")
p

(max(na.omit(d1_full$tmax_max))-min(na.omit(d1_full$tmax_max)))/24

tmax_results$tmax <- seq(min(na.omit(d1_full$tmax_max)),max(na.omit(d1_full$tmax_max)),0.6356964)

write.csv(tmax_results,file="/n/home04/nicholasarisco/PNAS -- Deforestation Paper/Predicting nonlinear terms/temp_results.csv")

(max(d1_full$tmax)-min(d1_full$tmax))/24
(max(d1_full$ONI)-min(d1_full$ONI))/36
(max(d1_full$p_tot)-min(d1_full$p_tot))/16


length(unique(d1_full$ONI))
