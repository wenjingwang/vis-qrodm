library(purrr)
library(quokar)
library(rggobi)

origin_model_gcd <- function(data, model, tau){
  idx_y = which(colnames(data) == all.vars(model)[1])
  idx_x = which(colnames(data) %in% all.vars(model)[-1])
  y <- as.matrix(data[idx_y])
  x <- as.matrix(data.frame(intercept = 1, data[idx_x]))
  gcd <- ALDqr_GCD(y, x, tau, iter = 1000, error = 1e-06)
  upper_bound <- mean(gcd) + 3*sd(gcd)
  lower_bound <- mean(gcd) - 3*sd(gcd)
  n <- length(y)
  outlier_flag <- rep(0, n)
  for(i in 1:n){
    if(gcd[i] >= upper_bound | gcd[i] <= lower_bound){
      outlier_flag[i] <- "outlier"
    }else{
      outlier_flag[i] <- "normal"
    }
  }
  actual_data <- cbind(data, gcd, outlier_flag)
  return(list(actual_data = actual_data, upper_bound = upper_bound,
              lower_bound = lower_bound))
}
##generate new obs#########
simvar <- function(x, n = 10, method = "grid") UseMethod("simvar")
simvar.factor <- function(x, n = 10, method = "grid"){
  switch(method,
         random = x[sample(length(x), n, replace = TRUE)],
         factor(levels(x), levels = levels(x))
  )
}

simvar.numeric <- function(x, n = 10, method = "grid"){
  rng <- range(x)
  switch(method,
         random = runif(n, rng[1], rng[2]),
         seq(rng[1], rng[2], length = n))
}

generate_data <- function(data, n = 1000, method = "grid"){
  if(method != "random"){
    n <- floor(n ^ (1/ncol(data)))
    df <- data.frame(expand.grid(lapply(data, simvar, n = n,
                                        method = "grid")))
    if(method == "nonaligned"){
      cont <- !sapply(df, is.factor)
      ranges <- lapply(df[, cont], function(x) diff(range(x)))
      df[,cont] <- df[,cont] +
        do.call(cbind, lapply(ranges, function(rng)
          runif(-rng/(2*n), rng/(2*n), n=nrow(df))))
    }
    df
  }else{
    data.frame(sapply(data, simvar, n=n, method=method))
  }
}

##calculate gcd and outlier_flag##################
sim_model_gcd <- function(data, model, tau, n = 1000){
  new_data <- generate_data(data, n = n, method = "random")
  colnames(new_data) <- colnames(data)
  m <- nrow(new_data)
  gcd <- rep(0, m)
  outlier_flag <- rep(0, m)
  upper_bound <- origin_model_gcd(data, model, tau)$upper_bound
  lower_bound <- origin_model_gcd(data, model, tau)$lower_bound
  for(i in 1:m){
    data_n <- plyr::rbind.fill(data, new_data[i, ])
    last_o <- nrow(data_n)
    idx_y = which(colnames(data_n) == all.vars(model)[1])
    idx_x = which(colnames(data_n) %in% all.vars(model)[-1])
    y <- as.matrix(data_n[idx_y])
    x <- as.matrix(data.frame(intercept = 1, data_n[idx_x]))
    gcd[i] <- ALDqr_GCD_i(y, x, tau, iter = 1000, error = 1e-06)
    if(gcd[i] >= upper_bound | gcd[i] <= lower_bound){
      outlier_flag[i] = "normal"
      }else{
      outlier_flag[i] = "outlier"
      }
  }
  sim_data <- cbind(new_data, gcd, outlier_flag)
  return(sim_data)
}
### Single variable
data <-  ais
model <- BMI ~ LBM
tau <-  0.5
### Multivariables
data <-  ais
model <- BMI ~ LBM + Ht
tau <-  0.5

actual_data <- origin_model_gcd(data, model, tau)$actual_data
sim_data <- sim_model_gcd(data, model, tau, n = 1000)

actual_data$.Type <- factor("actual")
sim_data$.Type <- factor("simulation")
dat <- plyr::rbind.fill(actual_data, sim_data)
g <- ggobi(dat)
d <- g[1]
glyph_color(d) <- as.numeric(dat$outlier_flag) + 1
glyph_type(d) <- ifelse(dat$.Type == "simulated", 1, 6)

###################

























