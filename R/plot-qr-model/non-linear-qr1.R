library(quantreg)
library(rggobi)


simulateSSlogis <- function(param, x){
  x$x1^2 * param[1] + x$x2^2 * param[2]
}


x1 <- runif(500, -3, 3)
x2 <- runif(500, -3, 3)

y = (x1^2 + x2^2 ) + rnorm(500, 0, 9)
realDat <- data.frame(x1, x2, y, flag = "actual")
resDat <- realDat

taus <- c(0.1,0.5,0.9)

for(tau in taus){
  realDat.nlrq <- nlrq(y ~ x1^2*a + x2^2*b, data = realDat, tau = tau,
                       start = list(a = 1, b = 1))
  param <- realDat.nlrq$m$getAllPars()
  simX <- generate_data(data.frame(x1, x2),n = 1000)
  simY <- simulateSSlogis(param, simX)
  simDat <- data.frame(simX, y = simY)
  simDat$flag = paste0("nlrq",tau)
  resDat <- rbind(resDat,simDat)
}

g <- ggobi(resDat)
d <- g[1]
glyph_color(d) <- ifelse(resDat$flag =="actual", 2,
                         ifelse(resDat$flag=="nlrq0.1", 5,
                                ifelse(resDat$flag=="nlrq0.5", 4,
                                       ifelse(resDat$flag=="nlrq0.9", 9, 12))))

