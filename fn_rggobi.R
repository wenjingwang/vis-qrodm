library(ggplot2)
library(ALDqr)
library(quantreg)
library(quokar)
library(tidyr)
library(dplyr)
library(rggobi)
####fn algorithm
fn_ggobi <- function(data, model, tau){
  p <- length(all.vars(model)[-1]) + 1
  data$weight <-  "small_w"
  object <- rq(model, tau, method = 'fn', data = data)
  fn <- frame_fn_obs(object, tau)
  for(i in seq(tau)){
  active <- fn[[i]]
  sep <- sort(active,TRUE)[p]
  idx <- which(active>=sep)
  data$weight[idx] <-  paste("large_w", tau[i], sep="")
  }
  data$size <- "small_w"
  data$size[which(data$weight != "small_w")] <- "large"
  return(data)
}


data <- ais
model <- BMI ~ LBM + Ht
tau <- c(0.1, 0.5, 0.9)
fn_data <- fn_ggobi(data, model, tau)

g <- ggobi(fn_data)
d <- g[1]
glyph_color(d) <- -1* as.numeric(as.factor(fn_data$weight)) +5
glyph_size(d) <-  (-1*as.numeric(as.factor(fn_data$size))+3)*3







