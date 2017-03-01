ALDqr_GCD_i <- function(y, x, tau, error, iter){
  n <- length(y)
  p <- ncol(x)
  theta <- ALDqr::EM.qr(y, x, tau, error, iter)$theta
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
  GCD_beta <- c(Q1_beta) %*% solve(-Q2_beta) %*% matrix(Q1_beta, ncol = 1)
  GCD_sigma <- Q1_sigma*solve(-Q2_sigma)*Q1_sigma
  GCD <- as.vector(GCD_beta + GCD_sigma)
  return(GCD)
}

########################
ALDqr_case_deletion_i <- function(y, x, tau, error, iter)
{
  n <- length(y)
  p <- ncol(x)
  qr <- ALDqr::EM.qr(y,x,tau,error,iter)
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
  beta_i <- beta_qr + taup2*solve(suma1)%*% E1
  sigma_i2 <-  sigma_qr^2 - solve(Q2_sigma)*E2/(2*sigma_qr^2)
  sigma_i <- sqrt(simplify2array(sigma_i2))
  theta_i <- list(beta_i = beta_i, sigma_i = sigma_i)
  return(theta_i)
}
#######################
ALDqr_QD_i <- function(y, x, tau, error, iter){
  p <- ncol(x)
  n <- length(y)
  theta_all <- ALDqr::EM.qr(y, x, tau, error, iter)$theta
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

