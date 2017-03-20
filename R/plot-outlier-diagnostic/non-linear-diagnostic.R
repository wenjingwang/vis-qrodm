library(quantreg)
library(rggobi)

cook.nlrq <- function(formula, data,tau){
  n <- nrow(data)
  object <- nlrq(formula, data = data, tau = tau)
  y_i<-sapply(1:n, function(z){
    nlrq2 <- nlrq(formula, data=data[-z,],tau=tau)
    coefNlrq2 <- nlrq2$m$getAllPars()
    y_ii <- t((coefNlrq2-coef(object))) %*% (coefNlrq2-coef(object))
    return(y_ii)
  })
  dev <- deviance(object)
  d_i <- y_i / ((p + 1) * dev^2)
  return(d_i)
}


cook.nlrq_i <- function(formula, data, tau){
  n <- nrow(data)
  object <- nlrq(formula, data = data, tau = tau)
  nlrq2 <- nlrq(formula, data = data[-n,], tau = tau)
  coefNlrq2 <-nlrq2$m$getAllPars()
  y_ii <- t((coefNlrq2-coef(object))) %*% (coefNlrq2-coef(object))
  dev <- deviance(nlrq2)
  d_i <- y_ii / ((p + 1) * dev^2)
  return(d_i)
}

origin_model_nl <- function(formula, data, tau){
  nl <- cook.nlrq(formula, data, tau)
  upper_bound <- mean(nl) + sd(nl)
  lower_bound <- mean(nl) - sd(nl)
  n <- nrow(data)
  outlier_flag <- rep(0, n)
  for(i in 1:n){
    if(nl[i] >= upper_bound | nl[i] <= lower_bound){
      outlier_flag[i] <- "outlier"
    }else{
      outlier_flag[i] <- "normal"
    }
  }
  actual_data <- cbind(data, nl, outlier_flag)
  return(list(actual_data = actual_data, upper_bound = upper_bound,
              lower_bound = lower_bound))
}


sim_model_nl <- function(formula, data, tau, n = 1000){
  new_data <- generate_data(data, n = n, method = "random")
  colnames(new_data) <- colnames(data)
  m <- nrow(new_data)
  nl <- rep(0, m)
  outlier_flag <- rep(0, m)
  upper_bound <- origin_model_nl(formula, data, tau)$upper_bound
  lower_bound <- origin_model_nl(formula, data, tau)$lower_bound
  for(i in 1:m){
    data_n <- plyr::rbind.fill(data, new_data[i, ])
    nl[i] <- cook.nlrq_i(formula, data_n, tau)
    if(nl[i] >= upper_bound | nl[i] <= lower_bound){
      outlier_flag[i] = "outlier"
    }else{
      outlier_flag[i] = "normal"
    }
  }
  sim_data <- cbind(new_data, nl, outlier_flag)
  return(sim_data)
}

Dat <- NULL
Dat$x <- rep(1:25, 20)
set.seed(1)
Dat$y <- SSlogis(Dat$x, 10, 12, 2)*rnorm(500, 1, 0.1)
Dat <- data.frame(x = Dat$x, y = Dat$y)

formula <- y ~ SSlogis(x, Asym, mid, scal)
tau <- 0.1
data <- Dat
p <- 3


actual_data <- origin_model_nl(formula, data, tau)$actual_data
sim_data <- sim_model_nl(formula, data, tau, n = 10000)

actual_data$.Type <- factor("actual")
sim_data$.Type <- factor("simulation")
dat <- plyr::rbind.fill(actual_data, sim_data)
g <- ggobi(dat)
d <- g[1]
gcolor <- ifelse(dat$outlier_flag == "outlier", 9,4)
glyph_color(d) <- gcolor
glyph_type(d) <- ifelse(dat$.Type == "simulated", 1, 6)

