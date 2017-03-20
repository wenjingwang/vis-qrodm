library(ggplot2)
library(ALDqr)
library(quantreg)
library(quokar)
library(tidyr)
library(dplyr)
library(rggobi)

br_ggobi <- function(data, model, tau){
  n <- nrow(data)
  idx_y = which(colnames(data) == all.vars(model)[1])
  idx_x = which(colnames(data) %in% all.vars(model)[-1])
  y <- as.matrix(data[idx_y])
  x <- as.matrix(data[idx_x])
  object <- rq(model, tau, method = 'br', data = data)
  ntau <- length(tau)
  br_flag <- rep("non-use", n)
  i=1
  for(i in 1:ntau){
    points <- frame_br(object, tau[i])$fitting_point$index
    br_flag[points] <- paste("use", tau[i], sep="")
  }
  br_data <- cbind(data, br_flag)
}

###########################################################
data <- ais
ais_female <- subset(ais, Sex == 1)
actual_data <- ais_female[, c('BMI', 'LBM', 'Bfat')]
actual_data$flag <- "observations"
taus <- c(0.1, 0.5, 0.9)
for(tau in taus){
  br_flag <- paste0("use",tau)
  br <- rq(BMI ~ LBM + Bfat, tau = tau, data = ais_female, method = "br")
  model <- BMI ~ LBM + Bfat
  br_data <- br_ggobi(data, model, tau)
  equation_data <- br_data[which(br_data$br_flag == br_flag), ]
  left_side <- data.frame(LBM = equation_data$LBM, Bfat = equation_data$Bfat)
  br_low <- generate_data(left_side, n = 10000, method = "random")
  right_side <- as.matrix(cbind(1, br_low)) %*% br$coef
  sim_data <- cbind(right_side,br_low)
  sim_data$flag <- paste0("plane",tau)
  colnames(sim_data) <- colnames(actual_data)
  actual_data <- plyr::rbind.fill(actual_data, sim_data)
  actual_data$flag[which(br_data$br_flag == br_flag)] <- paste0("plane",tau)
}


g <- ggobi(actual_data)
d <- g[1]
glyph_color(d) <- ifelse(actual_data$flag =="observations",2,
                         ifelse(actual_data$flag=="plane0.1",5,
                                ifelse(actual_data$flag=="plane0.5",4,
                                       ifelse(actual_data$flag=="plane0.9",9,12))))

