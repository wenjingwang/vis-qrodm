## ---- load
library(quokar)
library(quantreg)
library(tidyverse)
library(gridExtra)
library(purrr)

## ---- vis-if
F <- function(x){
  1-exp(-x)
}
f <- function(x){
  exp(-x)
}
Q <- function(p){
  -log(1-p)
}
x <- seq(0, 10, by = 0.001)
y <- sort(1-exp(-x))
inf_quantile <- function(y, tau){
  n <- length(y)
  n1 <- sum(y <= Q(tau))
  inf_q <- rep(tau/f(Q(tau)), n)
  inf_q[1:n1] <- rep((tau-1)/f(Q(tau)), n1)
  return(inf_q)
}
tau <- 0.1
data_q <- data.frame(inf_q = inf_quantile(y, 0.1), y = y)
p1 <- ggplot(data_q, aes(y, inf_q))+
  geom_line() + 
  geom_point(data = data.frame(a = Q(tau), b = (tau-1)/f(Q(tau))), aes(a, b), shape = 1) +
  geom_point(data = data.frame(a = Q(tau), b = tau/f(Q(tau))), aes(a, b), shape = 16) +
  geom_vline(xintercept = Q(tau), colour = "red", linetype = "longdash") +
  geom_hline(yintercept = 0, colour = "red", linetype = "longdash") +
  xlab("y") +
  ylab("Influence function for quantile (tau = 0.1)")
tau <- 0.5
data_q <- data.frame(inf_q = inf_quantile(y, 0.5), y = y)
p2 <- ggplot(data_q, aes(y, inf_q))+
  geom_line() + 
  geom_point(data = data.frame(a = Q(tau), b = (tau-1)/f(Q(tau))), aes(a, b), shape = 1) +
  geom_point(data = data.frame(a = Q(tau), b = tau/f(Q(tau))), aes(a, b), shape = 16) +
  geom_hline(yintercept = 0, colour = "red", linetype = "longdash") +
  geom_vline(xintercept = Q(tau), colour = "red", linetype = "longdash") +
  xlab("y") +
  ylab("Influence function for quantile (tau = 0.5)")
tau <- 0.9
data_q <- data.frame(inf_q = inf_quantile(y, 0.9), y = y)
y_outlier1 <- seq(0, Q(tau), by = 0.01)
y_outlier2 <- seq(Q(tau), 3, by = 0.01)
ny1 <- length(y_outlier1)
ny2 <- length(y_outlier2)
data_q_outlier1 <- data.frame(h = rep((tau-1)/f(Q(tau)), ny1), y = y_outlier1)
data_q_outlier2 <- data.frame(h = rep(tau/f(Q(tau)), ny2), y = y_outlier2)
data_q_outlier <- rbind(data_q_outlier1, data_q_outlier2)
p3 <- ggplot(data_q_outlier, aes(y, h))+
  geom_line() + 
  geom_point(data = data.frame(a = Q(tau), b = (tau-1)/f(Q(tau))), aes(a, b), shape = 1) +
  geom_point(data = data.frame(a = Q(tau), b = tau/f(Q(tau))), aes(a, b), shape = 16) +
  geom_hline(yintercept = 0, colour = "red", linetype = "longdash") +
  geom_vline(xintercept = Q(tau), colour = "red", linetype = "longdash") +
  xlab("y") +
  ylab("Influence function for quantile (tau = 0.9)")
grid.arrange(p1, p2, p3, nrow = 1)

## ---- qr-outlier
x <- sort(runif(100))
y1 <- 40*x + x*rnorm(100, 0, 10)
df <- data.frame(y1, x)
add_outlier <- data.frame(y1 = c(60,61,62), x = c(0.71, 0.73,0.75))
df_o <- rbind(df, add_outlier)

coef1 <- rq(y1 ~ x, tau = c(0.1, 0.5, 0.9), 
            data = df, method = "br")$coef
rq_coef1 <- data.frame(intercept = coef1[1, ], coef = coef1[2, ],
                       tau_flag =colnames(coef1))

coef2 <- rq(y1 ~ x, tau = c(0.1, 0.5, 0.9),
            data = df_o, method = "br")$coef
rq_coef2 <- data.frame(intercept = coef2[1, ], coef = coef2[2, ],
                       tau_flag =colnames(coef2))
p1 <- ggplot(df_o) +
  geom_point(aes(x = x, y = y1), alpha = 0.1) +
  geom_abline(data = rq_coef1, aes(intercept = intercept,
                                   slope = coef, colour = tau_flag))+
  geom_abline(data = rq_coef2, aes(intercept = intercept,
                                   slope = coef, colour = tau_flag))

#####
x <- sort(runif(100))
y2 <- 40*x + x*rnorm(100, 0, 10)
df <- data.frame(y2, x)
add_outlier <- data.frame(y2 = c(1,2,3), x = c(0.71, 0.73,0.75))
df_o <- rbind(df, add_outlier)

coef1 <- rq(y2 ~ x, tau = c(0.1, 0.5, 0.9), 
            data = df, method = "br")$coef
rq_coef1 <- data.frame(intercept = coef1[1, ], coef = coef1[2, ], tau_flag = colnames(coef1))

coef2 <- rq(y2 ~ x, tau = c(0.1, 0.5, 0.9), 
            data = df_o, method = "br")$coef
rq_coef2 <- data.frame(intercept = coef2[1, ], coef = coef2[2, ], tau_flag = colnames(coef2))

p2 <- ggplot(df_o) +
  geom_point(aes(x = x, y = y2), alpha = 0.1) +
  geom_abline(data = rq_coef1, aes(intercept = intercept,
                                   slope = coef, colour = tau_flag))+
  geom_abline(data = rq_coef2, aes(intercept = intercept,
                                   slope = coef, colour = tau_flag))

grid.arrange(p1, p2, nrow = 1)

## ---- move-y1
x <- sort(runif(100))
y <- 40*x + x*rnorm(100, 0, 10)
selectedX <- sample(50:100,5)
y2<- y
y2[selectedX] <- x[1:5]*rnorm(5, 0, 10)
y3 <- y2
y3[selectedX] <- y2[selectedX] - 10
y4 <- y3
y4[selectedX] <- y3[selectedX] - 10
df <- data.frame(x, y, y2, y3, y4)
df_m <- df %>% gather(variable, value, -x)
coefs <- 2:5 %>%
  map(~ rq(df[, .] ~ x, data = df, seq(0.1, 0.9, 0.1))) %>%
  map_df(~ as.data.frame(t(as.matrix(coef(.)))))
colnames(coefs) <- c("intercept", "slope")
variable <- rep(c("y", "y2", "y3", "y4"), each = 9)
tau_flag <- paste("tau=", rep(seq(0.1, 0.9, by = 0.1), 4), sep = "")
qr_lines <- data.frame(coefs, variable, tau_flag)
ggplot(df_m, aes(x = x, y=value)) +
  geom_point(alpha = 0.3) +
  geom_abline(data = qr_lines, aes(intercept = intercept, 
                                   slope =  slope, colour = tau_flag), 
              size = 0.8)+
  xlab("x") +
  ylab("y") +
  facet_wrap(~variable, ncol=2, scale = "free_y") +
  scale_colour_brewer(palette="YlOrRd")+
  theme_grey()

## ---- move-y2
tau <- rep(seq(0.1, 0.9, by = 0.1), 4)
model <- paste('rq', rep(1:4, each = 9), sep="")
df_m1 <- data.frame(model, tau, coefs)
df_mf <- df_m1 %>% gather(variable, value, -c(model, tau))
ggplot(df_mf, aes(x = tau, y = value, colour = model)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ variable, scale = "free_y") +
  xlab('quantiles') +
  ylab('coefficients')

## ---- move-y-multi1
n <- 100
set.seed(101)
x1 <- sort(rnorm(n, 0, 1))
x2 <- sort(rnorm(n, 1, 2))
y <- 40*(x1 + x2) + x1*rnorm(100, 0, 10) + x2*rnorm(100, 0, 10)
selectedX <- sample(50:100,5)
y2<- y
y2[selectedX] <- x1[1:5]*rnorm(5, 0, 10) + x2[1:5]*rnorm(5, 0, 10)
y3 <- y2
y3[selectedX] <- y2[selectedX] - 100
y4 <- y3
y4[selectedX] <- y3[selectedX] - 100
df <- data.frame(y, y2, y3, y4, x1, x2)
coefs <- 1:4 %>%
  map(~ rq(df[, .] ~ x1 + x2, data = df, seq(0.1, 0.9, 0.1))) %>%
  map_df(~ as.data.frame(t(as.matrix(coef(.)))))
colnames(coefs) <- c("intercept", "slope_x1", "slope_x2")
tau <- rep(seq(0.1, 0.9, by = 0.1), 4)
model <- paste('rq', rep(1:4, each = 9), sep="")
df_m1 <- data.frame(model, tau, coefs)
df_mf <- df_m1 %>% gather(variable, value, -c(model, tau))
ggplot(df_mf, aes(x = tau, y = value, colour = model)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ variable, scale = "free_y") +
  xlab('quantiles') +
  ylab('coefficients')

## ----  move-x1

x <- sort(runif(100))
y <- 40*x + x*rnorm(100, 0, 10)
selectedIdx <- sample(50:100,5)
df <- data.frame(y)
df$y2 <- y
df$x <- x
df$y2[selectedIdx] <- df$x[1:5]*rnorm(5, 0, 10)
df$x2 <- x
df$x2[selectedIdx] <- df$x2[selectedIdx] + 0.2
df$x3 <- df$x2
df$x3[selectedIdx] <- df$x3[selectedIdx] + 0.2
df$x4 <- df$x3
df$x4[selectedIdx] <- df$x4[selectedIdx] + 0.2
df_m <- df %>% gather(variable, value, -y, -y2)
coefs <- 3:6 %>%
  map(~ rq(df$y2 ~ df[, .], data = df, seq(0.1, 0.9, 0.1))) %>%
  map_df(~ as.data.frame(t(as.matrix(coef(.)))))
colnames(coefs) <- c("intercept", "slope")
variable <- rep(c("x", "x2", "x3", "x4"), each = 9)
tau_flag <- paste("tau=", rep(seq(0.1, 0.9, by = 0.1), 4), sep = "")
qr_lines <- data.frame(coefs, variable, tau_flag)
ggplot(df_m, aes(x = value, y=y2)) +
  geom_point(alpha = 0.3) +
  geom_abline(data = qr_lines, aes(intercept = intercept, 
                                   slope = slope, colour = tau_flag),
              size = 0.8) +
  xlab("x") +
  ylab("y") +
  facet_wrap(~variable, ncol=2, scale = "free") +
  scale_colour_brewer(palette="YlOrRd")+
  theme_grey()

## ---- move-x2
tau <- rep(seq(0.1, 0.9, by = 0.1), 4)
model <- paste('rq', rep(1:4, each = 9), sep="")
df_m1 <- data.frame(model, tau, coefs)
df_mf <- df_m1 %>% gather(variable, value, -c(model, tau))
ggplot(df_mf, aes(x = tau, y = value, colour = model)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ variable, scale = "free_y") +
  xlab('quantiles') +
  ylab('coefficients')

## ---- single-case
data(ais)
ais_female <- filter(ais, Sex == 1)
case <- 1 : nrow(ais_female)
ais_female <- cbind(case, ais_female)
coef_rq <- coef(rq(BMI ~ LBM, tau = c(0.1, 0.5, 0.9),
                   data = ais_female, method = "br"))

br_coef <- data.frame(intercept = coef_rq[1, ],
                      coef = coef_rq[2, ],
                      tau_flag = colnames(coef_rq))
ggplot(ais_female)+
  geom_point(aes(x = LBM, y = BMI)) +
  geom_abline(data = br_coef, aes(intercept = intercept,
                                  slope = coef,
                                  colour = tau_flag), size = 1) +
  geom_text(data = subset(ais_female, case %in% c(1, 75)),
            aes(x = LBM, y = BMI, label = case), 
            colour = "red",hjust = 0, vjust = 0) +
  scale_colour_brewer(palette="YlOrRd")+
  theme_grey()

## ---- multi-case
ais_female_f <- dplyr::select(ais_female, c(case, BMI, LBM, Bfat))
ais_female_f_long <- tidyr::gather(ais_female_f, variable, value, -case, -BMI)
ggplot(ais_female_f_long, aes(x = value, y = BMI))+
  geom_point(alpha = 0.5) +
  geom_text(data = subset(ais_female_f_long, case %in% c(56, 75)),
            aes(x = value, y = BMI, label = case), 
            colour = "red", vjust = 0, hjust = 0) +
  facet_wrap(~variable, scales = "free_x") +
  scale_colour_brewer(palette="YlOrRd")+
  theme_grey()

## ---- robust-distance
n <- nrow(object$model)
case <- rep(1:n, length(tau))
distance <- cbind(case, distance)
distance$residuals <- abs(distance$residuals)
tau_f <- paste("tau", tau, sep="")
text_flag <- 1:length(cutoff_h) %>%
  map(function(i){
    distance %>% 
      filter((residuals > cutoff_h[i] |rd > cutoff_v)
             & tau_flag == tau_f[i])})

text_flag_d <- rbind(text_flag[[1]], text_flag[[2]], text_flag[[3]])
ggplot(distance, aes(x = rd, y = residuals)) +
  geom_point() +
  geom_hline(data = data.frame(tau_flag = paste("tau", tau, sep=""), 
                               cutoff_h = cutoff_h),   
             aes(yintercept = cutoff_h), colour = "red") +
  geom_vline(xintercept = cutoff_v, colour = "red") +
  geom_text(data = text_flag_d, aes(label = case), hjust = 0, vjust = 0) +
  facet_wrap(~ tau_flag, scales = 'free_y') +
  xlab("Robust Distance") +
  ylab("|Residuals|")

## ---- generalized-cook
y <- ais_female$BMI
x <- cbind(1, ais_female$LBM, ais_female$Bfat)
case <- rep(1:length(y), length(tau))
GCD <- frame_mle(y, x, tau, error = 1e-06, iter = 10000,
                 method = 'cook.distance')
GCD_m <- cbind(case, GCD)
ggplot(GCD_m, aes(x = case, y = value )) +
  geom_point() +
  facet_wrap(~variable, scale = 'free_y') +
  geom_text(data = subset(GCD_m, value > mean(value) + 2*sd(value)),
            aes(label = case), hjust = 0, vjust = 0) +
  xlab("case number") +
  ylab("Generalized Cook Distance")

## ---- q-function
QD <- frame_mle(y, x, tau, error = 1e-06, iter = 10000,
                method = 'qfunction')
QD_m <- cbind(case, QD)
ggplot(QD_m, aes(x = case, y = value)) +
  geom_point() +
  facet_wrap(~variable, scale = 'free_y')+
  geom_text(data = subset(QD_m, value > mean(value) + sd(value)),
            aes(label = case), hjust = 0, vjust = 0) +
  xlab('case number') +
  ylab('Qfunction Distance')

## ---- bayes-prob
prob_m <- cbind(case, prob)
ggplot(prob_m, aes(x = case, y = value )) +
  geom_point() +
  facet_wrap(~variable, scale = 'free') +
  geom_text(data = subset(prob_m, value > mean(value) + 2*sd(value)),
            aes(label = case), hjust = 0, vjust = 0) +
  xlab("case number") +
  ylab("Mean probability of posterior distribution")

## ---- bayes-kl


kl_m <- cbind(case, kl)
ggplot(kl_m, aes(x = case, y = value)) +
  geom_point() +
  facet_wrap(~variable, scale = 'free')+
  geom_text(data = subset(kl_m, value > mean(value) + sd(value)),
            aes(label = case), hjust = 0, vjust = 0) +
  xlab('case number') +
  ylab('Kullback-Leibler')


## ---- engel-plot
data(engel)
engel$case <- 1:nrow(engel)
influential_points <- filter(engel, income > 2500)
coef_rq <- coef(rq(foodexp ~ income, 
                   tau = c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95),
                   data = engel, method = "br"))
br_coef <- data.frame(intercept = coef_rq[1, ],
                      coef = coef_rq[2, ],
                      tau_flag = colnames(coef_rq))
ggplot(engel)+
  geom_point(aes(x = income, y = foodexp)) +
  geom_abline(data = br_coef, aes(intercept = intercept,
                                  slope = coef,
                                  colour = tau_flag), size = 1) +
  geom_text(data = influential_points, aes(x = income, y = foodexp,
                                           label = case),hjust = 1, 
            vjust = 0) +
  scale_colour_brewer(palette="YlOrRd")+
  theme_grey()

## --- engel-diagnosing-robust
tau <- c(0.05, 0.1, 0.25,0.75, 0.9, 0.95)
object <- rq(foodexp ~ income, data = engel, tau = tau)
plot_distance <- frame_distance(object, tau = tau)
distance <- plot_distance[[1]]
cutoff_v <- plot_distance[[2]]; cutoff_v
cutoff_h <- plot_distance[[3]]; cutoff_h
distance$residuals <- abs(distance$residuals)
n <- nrow(object$model)
case <- rep(1:n, length(tau))
distance <- cbind(case, distance)
tau_f <- paste("tau", tau, sep="")
text_flag <- 1:length(cutoff_h) %>%
  map(function(i){
    distance %>%
      filter(tau_flag == tau_f[i]) %>%
      filter(residuals < cutoff_h[i] & rd > cutoff_v)
  })
text_flag_leverage <- rbind(text_flag[[1]], text_flag[[2]], text_flag[[3]],
                     text_flag[[4]], text_flag[[5]], text_flag[[6]])
text_flag2 <- 1:length(cutoff_h) %>%
  map(function(i){
    distance %>%
      filter(tau_flag == tau_f[i]) %>%
      filter(residuals > cutoff_h[i] & rd > cutoff_v)
  })

text_flag_influential <- rbind(text_flag2[[1]], text_flag2[[2]], 
                               text_flag2[[3]],
                            text_flag2[[4]], 
                            text_flag2[[5]], 
                            text_flag2[[6]])
ggplot(distance, aes(x = rd, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(data = data.frame(tau_flag = paste("tau", tau, sep=""), 
                               cutoff_h = cutoff_h),   
             aes(yintercept = cutoff_h), colour = "black") +
  geom_vline(xintercept = cutoff_v, colour = "black") +
  geom_point(data = text_flag_leverage, colour = 'blue') +
  geom_text(data = text_flag_leverage, aes(label = case), 
             hjust = 1, vjust = 1,
              colour = 'blue') +
  geom_point(data = text_flag_influential, colour = 'red')+
  geom_text(data = text_flag_influential, aes(label = case),
            hjust = 1, vjust = 1, colour = 'red')+
  facet_wrap(~ tau_flag, scales = 'free_y') +
  xlab("Robust Distance") +
  ylab("|Residuals|")


## ---- engel-diagnosing-gcd
tau <- c(0.05, 0.1, 0.25,0.75, 0.9, 0.95)
y <- engel$foodexp
x <- cbind(1, engel$income)
case <- rep(1:length(y), length(tau))
GCD <- frame_mle(y, x, tau, error = 1e-06, iter = 10000,
                 method = 'cook.distance')
GCD_m <- cbind(case, GCD)
ggplot(GCD_m, aes(x = case, y = value )) +
  geom_point() +
  facet_wrap(~variable, scale = 'free_y') +
  geom_text(data = subset(GCD_m, value > mean(value) + sd(value)),
            aes(label = case), hjust = 0, vjust = 0, colour = 'red') +
  geom_point(data = subset(GCD_m, value > mean(value) + sd(value)),
            colour = 'red') +
  xlab("case number") +
  ylab("Generalized Cook Distance")

## ---- engel-diagnosing-qd
QD <- frame_mle(y, x, tau, error = 1e-06, iter = 10000,
                method = 'qfunction')
QD_m <- cbind(case, QD)
ggplot(QD_m, aes(x = case, y = value)) +
  geom_point() +
  facet_wrap(~variable, scale = 'free_y')+
  geom_text(data = subset(QD_m, value > mean(value) + sd(value)),
            aes(label = case), hjust = 0, vjust = 0, colour = 'red') +
  geom_point(data = subset(QD_m, value > mean(value) + sd(value)),
             colour = 'red') +
  xlab('case number') +
  ylab('Qfunction Distance')

## ---- engel-diagnosing-prob
prob <- frame_bayes(y, x, tau, M =  500, burn = 100, 
                    method = 'bayes.prob')

kl <- frame_bayes(y, x, tau, M = 500, burn = 100,
                  method = 'bayes.kl')
prob_m <- cbind(case, prob)
ggplot(prob_m, aes(x = case, y = value )) +
  geom_point(alpha = 0.3) +
  geom_text(data = subset(prob_m, value > mean(value) + 2*sd(value)),
            aes(label = case), hjust = 1, vjust = 1, colour = 'red') +
  geom_point(data = subset(prob_m, value > mean(value) + 2*sd(value)),
             colour = 'red') +
  facet_wrap(~variable, scale = 'free') +
  xlab("case number") +
  ylab("Mean probability of posterior distribution")
## ---- engel-diagnosing-kl
kl_m <- cbind(case, kl)
ggplot(kl_m, aes(x = case, y = value )) +
  geom_point(alpha = 0.3) +
  geom_text(data = subset(kl_m, value > mean(value) + 2*sd(value)),
            aes(label = case), hjust = 1, vjust = 0, colour = 'red') +
  geom_point(data = subset(kl_m, value > mean(value) + 2*sd(value)),
             colour = 'red') +
  facet_wrap(~variable, scale = 'free') +
  xlab("case number") +
  ylab("Kullback-Leibler")

## ---- baseball-plot
data(baseball)
baseball <- baseball %>% filter(Salary != 'NA')
baseball$case <- 1:nrow(baseball)
baseball_long <- baseball %>%
  select(c(case, Name, Homer, Years, Salary)) %>%
  gather(variable, value, -c(case, Name, Salary))
text_flag_homer <- baseball %>%
  select(c(case, Name, Homer, Salary)) %>%
  filter(Homer > 35 | Salary > 2000) 
colnames(text_flag_homer) <- c("case", "Name", "value", "Salary")
text_flag_homer$variable <- "Homer"
text_flag_years <- baseball %>%
  select(c(case, Name, Years, Salary)) %>%
  filter(Years > 20 | Salary > 2000)
colnames(text_flag_years) <- c("case", "Name", "value", "Salary")
text_flag_years$variable <- "Years"
text_flag <- rbind(text_flag_homer, text_flag_years)
ggplot(baseball_long, aes(x = value, y = Salary)) +
  geom_point() +
  geom_text(data = text_flag, aes(label = case), colour = 'red',
            hjust = 0.5, vjust = 1) +
  geom_point(data = text_flag, colour = 'red') +
  facet_wrap(~variable, scales = "free_x") +
  xlab("x") +
  ylab("Salary")
## ---- baseball-diagnosing-robust
tau <- c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95)
object <- rq(Salary ~ Homer + Years, data = baseball, 
             tau = tau, method = 'br')
plot_distance <- frame_distance(object, 
                                tau = c(0.05, 0.1, 0.25,
                                        0.75, 0.9, 0.95))
distance <- plot_distance[[1]]
cutoff_v <- plot_distance[[2]]; cutoff_v
cutoff_h <- plot_distance[[3]]; cutoff_h
distance$residuals <- abs(distance$residuals)
n <- nrow(object$model)
case <- rep(1:n, length(tau))
distance <- cbind(case, distance)
tau_f <- paste("tau", tau, sep="")
text_flag <- 1:length(cutoff_h) %>%
  map(function(i){
    distance %>%
      filter(tau_flag == tau_f[i]) %>%
      filter(residuals > cutoff_h[i] | rd > cutoff_v)
  })
text_flag_d <- rbind(text_flag[[1]], text_flag[[2]], text_flag[[3]],
                     text_flag[[4]], text_flag[[5]], text_flag[[6]])
ggplot(distance, aes(x = rd, y = residuals)) +
  geom_point() +
  geom_hline(data = data.frame(tau_flag = paste("tau", tau, sep=""), 
                               cutoff_h = cutoff_h),   
             aes(yintercept = cutoff_h), colour = "blue") +
  geom_vline(xintercept = cutoff_v, colour = "blue") +
  geom_text(data = text_flag_d, aes(label = case), hjust = 1, vjust = 1,
            colour = 'red') +
  geom_point(data = text_flag_d, colour = 'red') +
  facet_wrap(~ tau_flag, scales = 'free_y') +
  xlab("Robust Distance") +
  ylab("|Residuals|")
## ---- baseball-diagnosing-gcd
tau <- c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95)
y <- baseball$Salary
x <- cbind(1, baseball$Homer, baseball$Years)
case <- rep(1:length(y), length(tau))
GCD <- frame_mle(y, x, tau, error = 1e-06, iter = 10000,
                 method = 'cook.distance')
GCD_m <- cbind(case, GCD)
ggplot(GCD_m, aes(x = case, y = value )) +
  geom_point() +
  facet_wrap(~variable, scale = 'free_y') +
  geom_text(data = subset(GCD_m, value > mean(value) + sd(value)),
            aes(label = case), hjust = 0, vjust = 0) +
  xlab("case number") +
  ylab("Generalized Cook Distance")
## ---- baseball-diagnosing-qd
QD <- frame_mle(y, x, tau, error = 1e-06, iter = 10000,
                method = 'qfunction')
QD_m <- cbind(case, QD)
ggplot(QD_m, aes(x = case, y = value)) +
  geom_point() +
  facet_wrap(~variable, scale = 'free_y')+
  geom_text(data = subset(QD_m, value > mean(value) + 1.8*sd(value)),
            aes(label = case), hjust = 0, vjust = 0) +
  xlab('case number') +
  ylab('Qfunction Distance')
## ---- baseball-diagnosing-prob
prob <- frame_bayes(y, x, tau, M =  500, burn = 100, 
                    method = 'bayes.prob')
prob_m <- cbind(case, prob)
ggplot(prob_m, aes(x = case, y = value )) +
  geom_point() +
  facet_wrap(~variable, scale = 'free') +
  geom_text(data = subset(prob_m, value > mean(value) + 2*sd(value)),
            aes(label = case), hjust = 0, vjust = 0) +
  xlab("case number") +
  ylab("Mean probability of posterior distribution")
## ---- baseball-diagnosing-kl
kl <- frame_bayes(y, x, tau, M = 50, burn = 10,
                  method = 'bayes.kl')
kl_m <- cbind(case, kl)
ggplot(kl_m, aes(x = case, y = value )) +
  geom_point() +
  facet_wrap(~variable, scale = 'free') +
  geom_text(data = subset(kl_m, value > mean(value) + 2*sd(value)),
            aes(label = case), hjust = 0, vjust = 0) +
  xlab("case number") +
  ylab("Kullback-Leibler")
## ---- trout-plot
fish_habit1 <- nlrq(density ~ exp(beta0+beta1*wdratio), data=trout, 
                    start = list(beta0=0, beta1 = 0),
                    tau=0.05, trace=TRUE)
fish_habit2 <- nlrq(density ~ exp(beta0+beta1*wdratio), data=trout, 
                    start = list(beta0=0, beta1 = 0),
                    tau=0.1, trace=TRUE)
fish_habit3 <- nlrq(density ~ exp(beta0+beta1*wdratio), data=trout, 
                    start = list(beta0=0, beta1 = 0),
                    tau=0.25, trace=TRUE)
fish_habit4 <- nlrq(density ~ exp(beta0+beta1*wdratio), data=trout, 
                    start = list(beta0=0, beta1 = 0),
                    tau=0.5, trace=TRUE)
fish_habit5 <- nlrq(density ~ exp(beta0+beta1*wdratio), data=trout, 
                    start = list(beta0=0, beta1 = 0),
                    tau=0.75, trace=TRUE)
fish_habit6 <- nlrq(density ~ exp(beta0+beta1*wdratio), data=trout, 
                    start = list(beta0=0, beta1 = 0),
                    tau=0.9, trace=TRUE)
wdratio_pre <- seq(from = 8, to = 55, length = 71)
line1 <- data.frame(wdratio = wdratio_pre, 
                    density = predict(fish_habit1, 
                                      newdata=list(x=wdratio_pre)),
                    line_flag = 'tau=0.05')
line2 <- data.frame(wdratio = wdratio_pre, 
                    density = predict(fish_habit2, 
                                      newdata=list(x=wdratio_pre)),
                    line_flag = 'tau=0.1')
line3 <- data.frame(wdratio = wdratio_pre, 
                    density = predict(fish_habit3, 
                                      newdata=list(x=wdratio_pre)),
                    line_flag = 'tau=0.25')
line4 <- data.frame(wdratio = wdratio_pre, 
                    density = predict(fish_habit4, 
                                      newdata=list(x=wdratio_pre)),
                    line_flag = 'tau=0.5')
line5 <- data.frame(wdratio = wdratio_pre, 
                    density = predict(fish_habit6, 
                     newdata=list(x=wdratio_pre)),
                    line_flag = 'tau=0.75')
line6 <- data.frame(wdratio = wdratio_pre, 
                    density = predict(fish_habit7, 
                                      newdata=list(x=wdratio_pre)),
                    line_flag = 'tau=0.9')
Dat_line <- rbind(line1, line2, line3, line4, line5, line6)
Dat$case <- 1:nrow(Dat)
ggplot(Dat, aes(wdratio, density)) +
  geom_point() +
  geom_text(aes(label = case)) +
  geom_line(data = Dat_line, aes(wdratio, density,colour = line_flag),
            size = 1)
## --- trout-diagnosing-robust
x <- matrix(trout$wdratio, ncol = 1)
resid1 <- matrix(resid(fish_habit1), ncol = 1)
resid2 <- matrix(resid(fish_habit2), ncol = 1)
resid3 <- matrix(resid(fish_habit3), ncol = 1)
resid4 <- matrix(resid(fish_habit4), ncol = 1)
resid5 <- matrix(resid(fish_habit5), ncol = 1)
resid6 <- matrix(resid(fish_habit6), ncol = 1)
resid <- cbind(resid1, resid2, resid3, resid4, resid5, resid6)
tau = c(0.05, 0.1, 0.25,0.5, 0.75, 0.9)
plot_distance <- frame_distance_implement(x, resid, 
                                          tau = c(0.05, 0.1, 0.25,
                                                  0.5, 0.75, 0.9))
distance <- plot_distance[[1]]
cutoff_v <- plot_distance[[2]]; cutoff_v
cutoff_h <- plot_distance[[3]]; cutoff_h
distance$residuals <- abs(distance$residuals)
case <- 1:nrow(trout)
distance <- cbind(case, distance)
head(distance)
tau_f <- paste("tau", tau, sep="")
text_flag <- 1:length(cutoff_h) %>%
  map(function(i){
    distance %>%
      filter(tau_flag == tau_f[i]) %>%
      filter(residuals > cutoff_h[i] | rd > cutoff_v)
  })
text_flag_leverage <- rbind(text_flag[[1]], text_flag[[2]], text_flag[[3]],
                            text_flag[[4]], text_flag[[5]], text_flag[[6]])
text_flag2 <- 1:length(cutoff_h) %>%
  map(function(i){
    distance %>%
      filter(tau_flag == tau_f[i]) %>%
      filter(residuals > cutoff_h[i] & rd > cutoff_v)
  })

text_flag_influential <- rbind(text_flag2[[1]], 
                               text_flag2[[2]], 
                               text_flag2[[3]],
                               text_flag2[[4]], 
                               text_flag2[[5]], 
                               text_flag2[[6]])
ggplot(distance, aes(x = rd, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(data = data.frame(tau_flag = paste("tau", tau, sep=""), 
                               cutoff_h = cutoff_h),   
             aes(yintercept = cutoff_h), colour = "black") +
  geom_vline(xintercept = cutoff_v, colour = "black") +
  geom_point(data = text_flag_leverage, colour = 'blue') +
  geom_text(data = text_flag_leverage, aes(label = case), 
            hjust = 1, vjust = 1,
            colour = 'blue') +
  geom_point(data = text_flag_influential, colour = 'red')+
  geom_text(data = text_flag_influential, aes(label = case),
            hjust = 0, vjust = 1, colour = 'red')+
  facet_wrap(~ tau_flag, scales = 'free_y') +
  xlab("Robust Distance") +
  ylab("|Residuals|")

