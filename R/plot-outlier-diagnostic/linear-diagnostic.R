##library
library(MASS)
library(quantreg)
library(ggplot2)
library(gridExtra)
library(purrr)
library(tidyr)
library(dplyr)
library(rggobi)

## generate data
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

######functions used for modeling
ALDqr_GCD_i <- function(y, x, tau, error, iter){
  n <- length(y)
  p <- ncol(x)
  theta <- EM.qr(y, x, tau, error, iter)$theta
  beta_qr <- theta[1:p, ]
  sigma_qr <- theta[p+1]
  taup2 <- (2/(tau * (1 - tau)))
  thep <- (1 - 2 * tau) / (tau * (1 - tau))
  delta2 <- (y - x %*% beta_qr)^2/(taup2 * sigma_qr)
  gamma2 <- (2 + thep^2/taup2)/sigma_qr
  muc <- y - x %*% beta_qr
  vchpN <- besselK(sqrt(delta2 * gamma2), 0.5 - 1)/
    (besselK(sqrt(delta2* gamma2), 0.5))*(sqrt(delta2 / gamma2))^(-1)
  vchp1 <- besselK(sqrt(delta2 * gamma2), 0.5 + 1)/
    (besselK(sqrt(delta2*gamma2), 0.5))*(sqrt(delta2 / gamma2))
  suma2 <- x[-n,] * c(vchpN[-n] * (y[-n] - x[-n,] %*% beta_qr) - thep)
  E1 <- apply(suma2, 2, sum)/(taup2)
  muc_i <- y[-n] - x[-n, ]%*%beta_qr
  E2 <- sum(3*sigma_qr - (vchpN[-n] * muc_i^2 -
                            2 * muc_i * thep + vchp1[-n] *(thep^2 + 2 * taup2))/taup2)
  Q1_beta <- E1/sigma_qr
  Q1_sigma <- -E2/(2*sigma_qr^2)
  xM <- c(sqrt(vchpN)) * x
  suma1 <- t(xM) %*% (xM)
  Q2_beta <- -(suma1)/(sigma_qr * taup2)
  Q2_sigma <- 3/(2*sigma_qr^2) - sum((vchpN*muc^2-2*muc*thep +
                                        vchp1*(thep^2 + 2*taup2)))/(sigma_qr^3*taup2)
  GCD_beta <- c(Q1_beta) %*% ginv(-Q2_beta) %*% matrix(Q1_beta, ncol = 1)
  GCD_sigma <- Q1_sigma*ginv(-Q2_sigma)*Q1_sigma
  GCD <- as.vector(GCD_beta + GCD_sigma)
  return(GCD)
}

ALDqr_case_deletion_i <- function(y, x, tau, error, iter)
{
  n <- length(y)
  p <- ncol(x)
  qr <- EM.qr(y,x,tau,error,iter)
  beta_qr <- qr$theta[1:p,]
  sigma_qr <- qr$theta[p+1]
  taup2 <- (2/(tau * (1 - tau)))
  thep <- (1 - 2 * tau) / (tau * (1 - tau))
  delta2 <- (y - x %*% beta_qr)^2/(taup2 * sigma_qr)
  gamma2 <- (2 + thep^2/taup2)/sigma_qr
  muc <- y - x %*% beta_qr
  vchp1 <- besselK(sqrt(delta2 * gamma2), 0.5 + 1)/(besselK(sqrt(delta2 *
                                                                   gamma2), 0.5)) * (sqrt(delta2 / gamma2))
  vchpN <- besselK(sqrt(delta2 * gamma2), 0.5 - 1)/(besselK(sqrt(delta2 *
                                                                   gamma2), 0.5)) * (sqrt(delta2 / gamma2))^(-1)
  suma2 <- x[-n,] * c(vchpN[-n] * (y[-n] - x[-n,] %*% beta_qr) - thep)
  E1 <- apply(suma2, 2, sum)/(taup2)
  muc_i <- y[-n] - x[-n, ]%*%beta_qr
  E2 <- sum(3*sigma_qr - (vchpN[-n] * muc_i^2 -
                            2 * muc_i * thep + vchp1[-n] *(thep^2 + 2 * taup2))/taup2)
  Q1_beta <- E1/sigma_qr
  Q1_sigma <- -E2/(2*sigma_qr^2)
  xM <- c(sqrt(vchpN)) * x
  suma1 <- t(xM) %*% (xM)
  Q2_beta <- -(suma1)/(sigma_qr * taup2)
  Q2_sigma <- 3/(2*sigma_qr^2) - sum((vchpN*muc^2-2*muc*thep +
                                        vchp1*(thep^2 + 2*taup2)))/(sigma_qr^3*taup2)
  beta_i <- beta_qr + taup2*ginv(suma1)%*% E1
  sigma_i2 <-  sigma_qr^2 - ginv(Q2_sigma)*E2/(2*sigma_qr^2)
  sigma_i <- sqrt(simplify2array(sigma_i2))
  theta_i <- list(beta_i = beta_i, sigma_i = sigma_i)
  return(theta_i)
}

ALDqr_QD_i <- function(y, x, tau, error, iter){
  p <- ncol(x)
  n <- length(y)
  theta_all <- EM.qr(y, x, tau, error, iter)$theta
  beta_all <- theta_all[1:p, ]
  sigma_all <- theta_all[p+1]
  beta_i <- as.vector(ALDqr_case_deletion_i(y, x, tau, error, iter)$beta_i)
  sigma_i <-  as.vector(ALDqr_case_deletion_i(y, x, tau, error, iter)$sigma_i)
  Q_function <- function(beta_qr, sigma_qr, tau){
    taup2 <- (2/(tau * (1 - tau)))
    thep <- (1 - 2 * tau) / (tau * (1 - tau))
    delta2 <- (y - x %*% beta_qr)^2/(taup2 * sigma_qr)
    gamma2 <- (2 + thep^2/taup2)/sigma_qr
    muc <- y - x %*% beta_qr
    vchpN <- besselK(sqrt(delta2 * gamma2), 0.5 - 1)/(besselK(sqrt(delta2 *
                                                                     gamma2), 0.5)) * (sqrt(delta2 / gamma2))^(-1)
    vchp1 <- besselK(sqrt(delta2 * gamma2), 0.5 + 1)/(besselK(sqrt(delta2 *
                                                                     gamma2), 0.5)) * (sqrt(delta2 / gamma2))
    Q <- (-3*log(sigma_qr)/2)*n - sum(vchpN * muc^2 - 2 * muc * thep +
                                        vchp1 *(thep^2 + 2 * taup2))/(2 * sigma_qr * taup2)
    return(Q)
  }
  Q_all <- Q_function(beta_all, sigma_all, tau)
  Q_i <- Q_function(beta_i, sigma_i, tau)
  QD <- 2*(Q_all - Q_i)
  return(QD)
}

ALDqr_QD <- function(y, x, tau, error, iter)
{
  p <- ncol(x)
  n <- length(y)
  theta_all <- EM.qr(y, x, tau, error, iter)$theta
  beta_all <- theta_all[1:p, ]
  sigma_all <- theta_all[p+1]
  beta_i <- ALDqr_case_deletion(y, x, tau, error, iter)$beta_i
  sigma_i <-  ALDqr_case_deletion(y, x, tau, error, iter)$sigma_i
  Q_function <- function(beta_qr, sigma_qr, tau){
    n <- length(y)
    taup2 <- (2/(tau * (1 - tau)))
    thep <- (1 - 2 * tau) / (tau * (1 - tau))
    delta2 <- (y - x %*% beta_qr)^2/(taup2 * sigma_qr)
    gamma2 <- (2 + thep^2/taup2)/sigma_qr
    muc <- y - x %*% beta_qr
    vchpN <- besselK(sqrt(delta2 * gamma2), 0.5 - 1)/(besselK(sqrt(delta2 *
                                                                     gamma2), 0.5)) * (sqrt(delta2 / gamma2))^(-1)
    vchp1 <- besselK(sqrt(delta2 * gamma2), 0.5 + 1)/(besselK(sqrt(delta2 *
                                                                     gamma2), 0.5)) * (sqrt(delta2 / gamma2))
    Q <- (-3*log(sigma_qr)/2)*n - sum(vchpN * muc^2 - 2 * muc * thep +
                                        vchp1 *(thep^2 + 2 * taup2))/(2 * sigma_qr * taup2)
    return(Q)
  }
  Q_all <- Q_function(beta_all, sigma_all, tau)
  Q_i <- rep(0, n)
  for(i in 1:n){
    Q_i[i] <- Q_function(beta_i[,i], sigma_i[i], tau)
  }
  QD <- rep(0, n)
  for(i in 1:n){
    QD[i] <- 2*(Q_all - Q_i[i])
  }
  return(QD)
}

ALDqr_GCD <- function(y, x, tau, error, iter)
{
  n <- length(y)
  p <- ncol(x)
  theta <- EM.qr(y, x, tau, error, iter)$theta
  beta_qr <- theta[1:p, ]
  sigma_qr <- theta[p+1]
  taup2 <- (2/(tau * (1 - tau)))
  thep <- (1 - 2 * tau) / (tau * (1 - tau))
  delta2 <- (y - x %*% beta_qr)^2/(taup2 * sigma_qr)
  gamma2 <- (2 + thep^2/taup2)/sigma_qr
  muc <- y - x %*% beta_qr
  vchpN <- besselK(sqrt(delta2 * gamma2), 0.5 - 1)/
    (besselK(sqrt(delta2* gamma2), 0.5))*(sqrt(delta2 / gamma2))^(-1)
  
  vchp1 <- besselK(sqrt(delta2 * gamma2), 0.5 + 1)/(besselK(sqrt(delta2              *gamma2), 0.5)) * (sqrt(delta2 / gamma2))
  E1 <- matrix(0, nrow = p, ncol = n)
  for(i in 1:n){
    suma2 <- x[-i,] * c(vchpN[-i] * (y[-i] - x[-i,] %*% beta_qr) - thep)
    E1[,i] <- apply(suma2, 2, sum)/(taup2)
  }
  E2 <- 1: n %>%
    map(function(i) {
      muc_i <- y[-i] - x[-i, ]%*%beta_qr
      sum(3*sigma_qr - (vchpN[-i] * muc_i^2 -
                          2 * muc_i * thep + vchp1[-i] *(thep^2 + 2 * taup2))/taup2)
    })
  E2 <- simplify2array(E2)
  Q1_beta <- E1/sigma_qr
  Q1_sigma <- -E2/(2*sigma_qr^2)
  xM <- c(sqrt(vchpN)) * x
  suma1 <- t(xM) %*% (xM)
  Q2_beta <- -(suma1)/(sigma_qr * taup2)
  Q2_sigma <- 3/(2*sigma_qr^2) - sum((vchpN*muc^2-2*muc*thep +
                                        vchp1*(thep^2 + 2*taup2)))/(sigma_qr^3*taup2)
  GCD_beta <- 1:n %>%
    map(function(i){
      c(Q1_beta[,i]) %*% ginv(-Q2_beta) %*% matrix(Q1_beta[,i], ncol = 1)
    })
  GCD_beta <- simplify2array(GCD_beta)
  GCD_sigma <- 1:n %>%
    map(function(i) {
      Q1_sigma[i]*ginv(-Q2_sigma)*Q1_sigma[i]
    })
  GCD_sigma <- simplify2array(GCD_sigma)
  GCD <- GCD_beta + GCD_sigma
  return(GCD)
}

ALDqr_case_deletion <- function(y, x, tau, error, iter)
{
  n <- length(y)
  p <- ncol(x)
  qr <- EM.qr(y,x,tau,error,iter)
  beta_qr <- qr$theta[1:p,]
  sigma_qr <- qr$theta[p+1]
  taup2 <- (2/(tau * (1 - tau)))
  thep <- (1 - 2 * tau) / (tau * (1 - tau))
  delta2 <- (y - x %*% beta_qr)^2/(taup2 * sigma_qr)
  gamma2 <- (2 + thep^2/taup2)/sigma_qr
  muc <- y - x %*% beta_qr
  vchp1 <- besselK(sqrt(delta2 * gamma2), 0.5 + 1)/(besselK(sqrt(delta2 *
                                                                   gamma2), 0.5)) * (sqrt(delta2 / gamma2))
  vchpN <- besselK(sqrt(delta2 * gamma2), 0.5 - 1)/(besselK(sqrt(delta2 *
                                                                   gamma2), 0.5)) * (sqrt(delta2 / gamma2))^(-1)
  E1 <- matrix(0, nrow = p, ncol = n)
  for(i in 1:n){
    suma2 <- x[-i,] * c(vchpN[-i] * (y[-i] - x[-i,] %*% beta_qr) - thep)
    E1[,i] <- apply(suma2, 2, sum)/(taup2)
  }
  E2 <- 1: n %>%
    map(function(i) {
      muc_i <- y[-i] - x[-i, ]%*%beta_qr
      sum(3*sigma_qr - (vchpN[-i] * muc_i^2 -
                          2 * muc_i * thep + vchp1[-i] *(thep^2 + 2 * taup2))/taup2)
    })
  E2 <- simplify2array(E2)
  Q1_beta <- E1/sigma_qr
  Q1_sigma <- -E2/(2*sigma_qr^2)
  xM <- c(sqrt(vchpN)) * x
  suma1 <- t(xM) %*% (xM)
  Q2_beta <- -(suma1)/(sigma_qr * taup2)
  Q2_sigma <- 3/(2*sigma_qr^2) - sum((vchpN*muc^2-2*muc*thep +
                                        vchp1*(thep^2 + 2*taup2)))/(sigma_qr^3*taup2)
  beta_i <- matrix(0, nrow=p, ncol = n)
  for(i in 1:n){
    beta_i[,i] <- beta_qr + taup2*ginv(suma1)%*% E1[,i]
  }
  sigma_i2 <- 1:n %>%
    map(function(i) sigma_qr^2 - ginv(Q2_sigma)*E2[i]/(2*sigma_qr^2))
  sigma_i <- sqrt(simplify2array(sigma_i2))
  theta_i <- list(beta_i = beta_i, sigma_i = sigma_i)
  return(theta_i)
}

###rggobi_result
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
      outlier_flag[i] = "outlier"
    }else{
      outlier_flag[i] = "normal"
    }
  }
  sim_data <- cbind(new_data, gcd, outlier_flag)
  return(sim_data)
}

origin_model_qd <- function(data, model, tau){
  idx_y = which(colnames(data) == all.vars(model)[1])
  idx_x = which(colnames(data) %in% all.vars(model)[-1])
  y <- as.matrix(data[idx_y])
  x <- as.matrix(data.frame(intercept = 1, data[idx_x]))
  qd <- ALDqr_QD(y, x, tau, iter = 1000, error = 1e-06)
  upper_bound <- mean(qd) + 3*sd(qd)
  lower_bound <- mean(qd) - 3*sd(qd)
  n <- length(y)
  outlier_flag <- rep(0, n)
  for(i in 1:n){
    if(qd[i] >= upper_bound | qd[i] <= lower_bound){
      outlier_flag[i] <- "outlier"
    }else{
      outlier_flag[i] <- "normal"
    }
  }
  actual_data <- cbind(data, qd, outlier_flag)
  return(list(actual_data = actual_data, upper_bound = upper_bound,
              lower_bound = lower_bound))
}
sim_model_qd <- function(data, model, tau, n = 1000){
  new_data <- generate_data(data, n = n, method = "random")
  colnames(new_data) <- colnames(data)
  m <- nrow(new_data)
  qd <- rep(0, m)
  outlier_flag <- rep(0, m)
  upper_bound <- origin_model_qd(data, model, tau)$upper_bound
  lower_bound <- origin_model_qd(data, model, tau)$lower_bound
  for(i in 1:m){
    data_n <- plyr::rbind.fill(data, new_data[i, ])
    last_o <- nrow(data_n)
    idx_y = which(colnames(data_n) == all.vars(model)[1])
    idx_x = which(colnames(data_n) %in% all.vars(model)[-1])
    y <- as.matrix(data_n[idx_y])
    x <- as.matrix(data.frame(intercept = 1, data_n[idx_x]))
    qd[i] <- ALDqr_QD_i(y, x, tau, iter = 1000, error = 1e-06)
    if(qd[i] >= upper_bound | qd[i] <= lower_bound){
      outlier_flag[i] = "outlier"
    }else{
      outlier_flag[i] = "normal"
    }
  }
  sim_data <- cbind(new_data, qd, outlier_flag)
  return(sim_data)
}

## EM.qr
EM.qr<-function(y,x=NULL,tau=NULL, error = 0.000001 ,iter=2000, 
                envelope=FALSE)
{
  
  if(envelope==TRUE){
    n <-length(y)
    
    #### Regressao Quantilica: Envelope   \rho_p(y-mu)/sigma^2 \sim exp(1)
    
    rq      <- EM.qr(y,x,tau)
    columas <-  ncol(x)
    muc     <- (y-x%*%rq$theta[1:columas])
    Ind     <- (muc<0)+0  
    d2s     <- muc*(tau-Ind)  ### Distancia de mahalobonisb
    d2s     <- sort(d2s)
    
    xq2  <- qexp(ppoints(n), 1/(rq$theta[4]))
    
    Xsim <- matrix(0,100,n)
    for(i in 1:100){
      Xsim[i,] <- rexp(n, 1/(rq$theta[4]))
    }
    
    Xsim2 <- apply(Xsim,1,sort)
    d21   <- matrix(0,n,1)
    d22   <- matrix(0,n,1)
    for(i in 1:n){
      d21[i]  <- quantile(Xsim2[i,],0.05)
      d22[i]  <- quantile(Xsim2[i,],0.95)
    }
    
    d2med <-apply(Xsim2,1,mean)
    
    fy <- range(d2s,d21,d22)
    plot(xq2,d2s,xlab = expression(paste("Theoretical ",exp(1), " quantiles")), 
         ylab="Sample values and simulated envelope",pch=20,ylim=fy)
    par(new=T)
    plot(xq2,d21,type="l",ylim=fy,xlab="",ylab="")
    par(new=T)
    plot(xq2,d2med,type="l",ylim=fy,xlab="",ylab="")
    par(new=T)
    plot(xq2,d22,type="l",ylim=fy,xlab="",ylab="")
    
  }
  
  
  
  
  ################################################################################
  ###                    MI Empirica: Veja givens
  ################################################################################
  
  MI_empirica<-function(y,x,tau,theta){
    
    p      <- ncol(x)
    n      <- nrow(x)
    taup2  <- (2/(tau*(1-tau)))
    thep   <- (1-2*tau)/(tau*(1-tau))
    
    beta   <- theta[1:p]
    sigma  <- theta[p+1]
    mu     <- x%*%beta
    muc    <- y-mu
    
    delta2 <- (y-x%*%beta)^2/(taup2*sigma)
    gamma2 <- (2+thep^2/taup2)/sigma
    K05P   <- 2*besselK(sqrt(delta2*gamma2), 0.5)*(sqrt(delta2/gamma2)^0.5)
    K05N   <- 2*besselK(sqrt(delta2*gamma2), -0.5)*(sqrt(delta2/gamma2)^(-0.5))
    K15P   <- 2*besselK(sqrt(delta2*gamma2), 1.5)*(sqrt(delta2/gamma2)^(1.5))  
    
    DerG   <- matrix(0,nrow=(p+1),ncol=(p+1))
    
    for (i in 1:n)
    {
      dkibeta   <- (muc[i]/(taup2*sigma))*(K05N[i])*x[i,]
      dkisigma  <- sqrt(delta2[i])/(2*sigma)*K05N[i]+sqrt(gamma2)/(2*sigma)*K15P[i]
      GradBeta  <- -thep/(taup2*sigma)*x[i,]+(K05P[i])^(-1)*dkibeta
      Gradsigma <- -1.5/sigma-thep*muc[i]/(taup2*sigma^2)+ (K05P[i])^(-1)*dkisigma      
      GradAux   <- as.matrix(c(GradBeta,Gradsigma),p+1,1)
      DerG      <- DerG+GradAux%*%t(GradAux)
    }
    
    EP         <- sqrt(diag(ginv(DerG)))
    
    obj.out <- list(EP = as.vector(EP))
    
    return(obj.out)    
    
  }
  ################################################################################
  
  ################################################################################
  ## Verossimilhanca da RQ: Usando a funcao de perda e usando a Bessel funcion
  ################################################################################
  
  logVQR   <- function(y,x,tau,theta)
  {
    p      <- ncol(x)
    n      <- nrow(x)
    
    beta   <- theta[1:p]
    sigma  <- theta[p+1]
    mu     <- x%*%beta
    muc    <- (y-mu)/sigma
    Ind    <- (muc<0)+0  
    logver <- sum(-log(sigma)+log(tau*(1-tau))-muc*(tau-Ind))
    return(logver)
  }
  
  p     <- ncol(x)
  n     <- nrow(x)
  reg   <- lm(y ~ x[,2:p])
  taup2 <- (2/(tau*(1-tau)))
  thep  <- (1-2*tau)/(tau*(1-tau))
  #Inicializa beta e sigma2 com os estimadores de minimos quadrados
  beta  <- as.vector(coefficients(reg),mode="numeric")
  sigma <- sqrt(sum((y-x%*%beta)^2)/(n-p))
  
  
  lk = lk1 = lk2   <- logVQR(y,x,tau,c(beta,sigma))## log-likelihood
  
  teta_velho <- matrix(c(beta,sigma),ncol=1)
  cont       <- 0
  criterio   <- 1
  
  while( criterio> error)
  { print(criterio)
    cont        <- (cont+1)
    
    muc         <- (y-x%*%beta)
    delta2      <- (y-x%*%beta)^2/(taup2*sigma)
    gamma2      <- (2+thep^2/taup2)/sigma
    
    vchpN       <- besselK(sqrt(delta2*gamma2), 0.5-1)/(besselK(sqrt(delta2*gamma2), 0.5))*(sqrt(delta2/gamma2))^(-1)
    vchp1       <- besselK(sqrt(delta2*gamma2), 0.5+1)/(besselK(sqrt(delta2*gamma2), 0.5))*(sqrt(delta2/gamma2))
    
    xM          <- c(sqrt(vchpN))*x  	
    suma1       <- t(xM)%*%(xM)
    suma2       <- x*c(vchpN*y-thep)		
    
    sigma       <- sum(vchpN*muc^2-2*muc*thep+vchp1*(thep^2+2*taup2))/(3*n*taup2)
    beta        <- ginv(suma1)%*%apply(suma2,2,sum)
    
    teta_novo   <- matrix(c(beta,sigma),ncol=1)
    criterio    <- sqrt(sum((teta_velho-teta_novo)^2)) 
    
    lk3         <- logVQR(y,x,tau,c(beta,sigma))
    if(cont<2) criterio <- abs(lk2 - lk3)/abs(lk3)
    else {
      tmp       <- (lk3 - lk2)/(lk2 - lk1)
      tmp2      <- lk2 + (lk3 - lk2)/(1-tmp)
      criterio  <- abs(tmp2 - lk3)
    }
    
    lk2 <- lk3
    
    if (cont==iter)
    { 
      break
    }
    
    
    teta_velho <- teta_novo
  }
  
  Weights <- vchpN*vchp1
  
  
  EP      <- MI_empirica(y,x,tau,teta_novo)$EP
  logver  <- logVQR(y,x,tau,teta_novo)
  return(list(theta=teta_novo,EP=EP,logver=logver,iter=cont,Weights=Weights,di=abs(muc)/sigma))
}

## Examples for plotting for 4D model in low quantile
data <-  subset(ais, Sex == 1)
model <- BMI ~ LBM + Bfat + Ferr
tau <-  0.1

actual_data <- origin_model_gcd(data, model, tau)$actual_data
sim_data <- sim_model_gcd(data, model, tau, n = 10000)
save(sim_data, file = "data/sim_gcd_4d_low.rdata")

actual_data$.Type <- factor("actual")
sim_data$.Type <- factor("simulation")
dat <- plyr::rbind.fill(actual_data, sim_data)
g <- ggobi(dat)
d <- g[1]
gcolor <- ifelse(dat$outlier_flag == "outlier", 9,4) # 9: purple, 4:green
glyph_color(d) <- gcolor
glyph_type(d) <- ifelse(dat$.Type == "simulated", 1, 6)

