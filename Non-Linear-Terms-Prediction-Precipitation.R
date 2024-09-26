library(ggplot2)

## Load Model ##
# m21_results <- readRDS("m21_model.rds")

# Create an empty data frame to store the results
median_table <- data.frame(Index = numeric(20), Value = numeric(20))
LCI_table <- data.frame(Index = numeric(20), Value = numeric(20))
UCI_table <- data.frame(Index = numeric(20), Value = numeric(20))

# Loop through the values of i
for (i in 1:20) {
  index_name <- paste("index.", i, sep = "")
  
  # Extract the median value for the specified index
  median_value <- median(as.data.frame(m21_RR_results$marginals.random$`inla.group(p_tot)`[[i]])$x)
  
  # Store the results in the data frame
  median_table[i, ] <- c(i, exp(median_value))
}

# Print or use the results_table as needed
plot(median_table)

# Loop through the values of i
for (i in 1:20) {
  # Replace the number after "index." with i
  index_name <- paste("index.", i, sep = "")
  
  # Extract the median value for the specified index
  LCI_value <- quantile(as.data.frame(m21_RR_results$marginals.random$`inla.group(p_tot)`[[index_name]])$x,0.05)
  
  # Store the results in the data frame
  LCI_table[i, ] <- c(i, exp(LCI_value))
}

# Print or use the results_table as needed
plot(LCI_table)


# Loop through the values of i
for (i in 1:20) {
  # Replace the number after "index." with i
  index_name <- paste("index.", i, sep = "")
  
  # Extract the median value for the specified index
  UCI_value <- quantile(as.data.frame(m21_RR_results$marginals.random$`inla.group(p_tot)`[[index_name]])$x,0.95)
  
  # Store the results in the data frame
  UCI_table[i, ] <- c(i, exp(UCI_value))
}

# Print or use the results_table as needed
plot(UCI_table)

colnames(median_table) <- c("id","median")
colnames(LCI_table) <- c("id","LCI")
colnames(UCI_table) <- c("id","UCI")

prec_CrI <- merge(LCI_table,UCI_table,by="id")
prec_results <- merge(median_table,prec_CrI,by="id")


p <- ggplot(prec_results) + 
  geom_line(aes(y=median, x=id, colour = "sin")) +
  geom_ribbon(aes(ymin=LCI, ymax=UCI, x=id, fill = "band"), alpha = 0.3) +
  scale_colour_manual("",values="blue") +
  scale_fill_manual("",values="grey12")
p

(max(na.omit(d1_full$p_tot))-min(na.omit(d1_full$p_tot)))/24
prec_results$p_tot <- seq(0,max(na.omit(d1_full$p_tot)),5.094178)

write.csv(prec_results,file="prec_results.csv")
