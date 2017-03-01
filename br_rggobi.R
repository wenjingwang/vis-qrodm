library(ggplot2)
library(ALDqr)
library(quantreg)
library(quokar)
library(tidyr)
library(dplyr)
library(rggobi)
###################
br_ggobi <- function(data, model, tau){
  n <- nrow(data)
  idx_y = which(colnames(data) == all.vars(model)[1])
  idx_x = which(colnames(data) %in% all.vars(model)[-1])
  y <- as.matrix(data[idx_y])
  x <- as.matrix(data[idx_x])
  object <- rq(model, tau, method = 'br', data = data)
  ntau <- length(tau)
  br_flag <- rep("non-use", n)
  for(i in 1:ntau){
    points <- frame_br(object, tau[i])$choose$index
    br_flag[points] <- paste("use", tau[i], sep="")
  }
  br_data <- cbind(data, br_flag)
}


br <- rq(BMI ~ LBM + Ht, tau = 0.5, data = ais, method = "br")
model <- BMI ~ LBM + Ht
data <- ais
tau <- 0.5
br_data <- br_ggobi(data, model, tau)
equation_data <- br_data[which(br_data$br_flag == 'use0.5'), ]
left_side <- data.frame(LBM = equation_data$LBM, Ht = equation_data$Ht)
x1_sim <- seq(min(left_side$LBM), max(left_side$LBM), by = 0.01)
x2_sim <- seq(min(left_side$Ht), max(left_side$Ht), by = 0.01)
#br_low <- generate_data(left_side, n = 1000, method = "random")
br_low <- data.frame(LBM = x1_sim, Ht = x2_sim)
right_side <- as.matrix(cbind(1, br_low)) %*% br$coef
sim_data <- cbind(br_low, right_side)
sim_data$flag <- "plane"
actual_data <- ais[, c('BMI', 'LBM', 'Ht')]
actual_data$flag <- "observations"
colnames(sim_data) <- colnames(actual_data)
br_plane <- plyr::rbind.fill(actual_data, sim_data)
br_plane$flag[which(br_data$br_flag == 'use0.5')] <- "plane"



g <- ggobi(br_plane)
d <- g[1]
glyph_color(d) <- as.numeric(as.factor(br_plane$flag))

glyph_type(d) <- as.numeric(br_data$br_flag)


install.packages("rgl")
library(rgl)

# Simulate some data
br2 <- rq(BMI ~ LBM + Ht, tau = 0.5, data = ais, method = "br")
br3 <- rq(BMI ~ LBM + Ht, tau = 0.9, data = ais, method = "br")

x1 <- runif(50)
x2 <- runif(50)
x3 <- rep(1, 50)


y <- matrix(c(x1, x2, x3), ncol = 3) %*% br$coef
y2 <- matrix(c(x1, x2, x3), ncol = 3) %*% br2$coef
y3 <- matrix(c(x1, x2, x3), ncol = 3) %*% br3$coef

