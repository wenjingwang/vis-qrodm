library(ggplot2)
library(ALDqr)
library(quantreg)
library(quokar)
library(tidyr)
library(dplyr)
library(rggobi)

###########################################################
data <- ais
ais_female <- subset(ais, Sex == 1)
actual_data <- ais_female[, c('BMI', 'LBM', 'Bfat')]
actual_data$flag <- "observations"
taus <- c(0.1, 0.5, 0.9)
model <- BMI ~ LBM + Bfat
for(tau in taus){
  
  br <- rq(model, tau = tau, data = ais_female, method = "br")
  
  br_low <- generate_data(ais_female[,c("LBM","Bfat")], n = 10000, method = "random")
  
  right_side <- as.matrix(cbind(1, br_low)) %*% br$coef
  
  sim_data <- cbind(right_side,br_low)
 
  sim_data$flag <- paste0("plane",tau)
  colnames(sim_data) <- colnames(actual_data)
  
  actual_data <- plyr::rbind.fill(actual_data, sim_data)
}


g <- ggobi(actual_data)
d <- g[1]

gc <- ifelse(actual_data$flag =="observations",2,
                         ifelse(actual_data$flag=="plane0.1",5,
                                ifelse(actual_data$flag=="plane0.5",4,
                                       ifelse(actual_data$flag=="plane0.9",9,12))))
gt <- c(rep(6,nrow(ais_female)),rep(7,nrow(actual_data)-nrow(ais_female)))
gs <- c(rep(4,nrow(ais_female)),rep(1,nrow(actual_data)-nrow(ais_female)))

glyph_color(d) <- gc
glyph_type(d) <- gt
glyph_size(d) <- gs
