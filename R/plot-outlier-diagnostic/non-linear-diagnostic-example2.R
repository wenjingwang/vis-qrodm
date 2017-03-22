#library(quantreg)

library(rggobi)

## Another non-linear example: \beta_{1} + \beta_{2}*exp(\beta_{3}x)

## Simulate original data


cook.nlrq <- function(formula, data, tau, start){
  n <- nrow(data)
  object <- nlrq(formula, data = data, tau = tau, start = start)
  y_i<-sapply(1:n, function(z){
    nlrq2 <- nlrq(formula, data=data[-z,],tau=tau, start = start)
    coefNlrq2 <- nlrq2$m$getAllPars()
    y_ii <- t((coefNlrq2-coef(object))) %*% (coefNlrq2-coef(object))
    return(y_ii)
  })
  dev <- deviance(object)
  d_i <- y_i / ((p + 1) * dev^2)
  return(d_i)
}


cook.nlrq_i <- function(formula, data, tau, start){
  n <- nrow(data)
  object <- nlrq(formula, data = data, tau = tau, start = start)
  nlrq2 <- nlrq(formula, data = data[-n,], tau = tau, start = start)
  coefNlrq2 <-nlrq2$m$getAllPars()
  y_ii <- t((coefNlrq2-coef(object))) %*% (coefNlrq2-coef(object))
  dev <- deviance(nlrq2)
  d_i <- y_ii / ((p + 1) * dev^2)
  return(d_i)
}

origin_model_nl <- function(formula, data, tau, start){
  nl <- cook.nlrq(formula, data, tau, start)
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


sim_model_nl <- function(formula, data, tau, start, n = 1000){
  new_data <- generate_data(data, n = n, method = "random")
  colnames(new_data) <- colnames(data)
  m <- nrow(new_data)
  nl <- rep(0, m)
  outlier_flag <- rep(0, m)
  upper_bound <- origin_model_nl(formula, data, tau, start)$upper_bound
  lower_bound <- origin_model_nl(formula, data, tau, start)$lower_bound
  for(i in 1:m){
    data_n <- plyr::rbind.fill(data, new_data[i, ])
    nl[i] <- cook.nlrq_i(formula, data_n, tau, start)
    if(nl[i] >= upper_bound | nl[i] <= lower_bound){
      outlier_flag[i] = "outlier"
    }else{
      outlier_flag[i] = "normal"
    }
  }
  sim_data <- cbind(new_data, nl, outlier_flag)
  return(sim_data)
}


x <- rep(1:25, 20)
y <- 15 + 1* exp(0.2*x)*rnorm(500, 0, 10)
Dat <- data.frame(x, y)
#ggplot(Dat, aes(x = x, y = y)) +
#  geom_point() 

formula <- y ~ a + b*exp(c*x)
data <- Dat
tau <-  0.9
p <- 3
start <- list(a = 1, b = 1, c = 1)


actual_data <- origin_model_nl(formula, data, tau, start)$actual_data
sim_data <- sim_model_nl(formula, data, tau, start, n = 10000)

actual_data$.Type <- factor("actual")
sim_data$.Type <- factor("simulation")
dat <- plyr::rbind.fill(actual_data, sim_data)
g <- ggobi(dat)
d <- g[1]
gcolor <- ifelse(dat$outlier_flag == "outlier", 9,4)
glyph_color(d) <- gcolor
glyph_type(d) <- ifelse(dat$outlier_flag == "oultier", 7, 6)


