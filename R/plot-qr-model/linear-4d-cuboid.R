data <- ais
ais_female <- subset(ais, Sex == 1)
actual_data <- ais_female[, c('BMI', 'LBM', 'Bfat', 'Ferr')]
actual_data$flag <- "observations"
taus <- c(0.1, 0.5, 0.9)
for(tau in taus){
  br_flag <- paste0("use",tau)
  br <- rq(BMI ~ LBM + Bfat + Ferr, tau = tau, data = ais_female)
  model <- BMI ~ LBM + Bfat + Ferr
  br_data <- br_ggobi(data, model, tau)
  equation_data <- br_data[which(br_data$br_flag == br_flag), ]
  left_side <- data.frame(LBM = equation_data$LBM, Bfat = equation_data$Bfat,
                          Ferr = equation_data$Ferr)
  br_low <- generate_data(left_side, n = 10000, method = "random")
  right_side <- as.matrix(cbind(1, br_low)) %*% br$coef
  sim_data <- cbind(right_side,br_low)
  sim_data$flag <- paste0("cuboid",tau)
  colnames(sim_data) <- colnames(actual_data)
  actual_data <- plyr::rbind.fill(actual_data, sim_data)
  actual_data$flag[which(br_data$br_flag == br_flag)] <- paste0("cuboid",tau)
}


g <- ggobi(actual_data)
d <- g[1]
glyph_color(d) <- ifelse(actual_data$flag =="observations",2,
                         ifelse(actual_data$flag=="cuboid0.1",5,
                                ifelse(actual_data$flag=="cuboid0.5",4,
                                       ifelse(actual_data$flag=="cuboid0.9",9,12))))