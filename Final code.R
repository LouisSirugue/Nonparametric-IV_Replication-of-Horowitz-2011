####################################################################################################
####################################################################################################
##                                                                                                ##
##                              Applied Nonparametric IV Estimation                               ##
##                                                                                                ##
##                       Evseeva Uliana - Lajoumard Foucauld - Sirugue Louis                      ##
##                                                                                                ##
####################################################################################################
####################################################################################################

rm(list=ls())

#install.packages("ggplot2", "MASS", "pracma", "dplyr")
library(ggplot2)
library(MASS)
library(pracma)
library(dplyr)

########################################
#Ill-posed inverse problem illustration#
########################################

#Setup
######

#1. Define the psy function
psi <- function(j, z) {
  psi_z <- vector(mode = "numeric", length = 0)
  if (j == 1) {
    psi_z <- matrix(1, length(z), 1)
  } else {
    for (i in (1:length(z))) {
      newval <- sqrt(2) * cos((j - 1) * pi * z[i])
      psi_z <- c(psi_z, newval)
    }
  }
  return(psi_z)
}

#2. Define the function lambda(j)
lambda <- function(j) {
  if (j == 1) {
    lambdaj <- 1
  } else {
    lambdaj <- 0.2 * ((j - 1)^(-4))
  }
  return(lambdaj)
}

#3. Define the function for Fourier coefficients g(j)
g <- function(j) {
  if (j == 1) {
    gj <- 0.5
  } else {
    gj <- sqrt(2) * (((-1)^(j - 1)) - 1) * ((pi * (j - 1))^(-2))
  }
  return(gj)
}

#4. Define the function c(j) 
cvalue <- function(j) {
  cj <- g(j) * sqrt(lambda(j))
  return(cj)
}

#5. Generate V and W
set.seed(1)
V <- rnorm(n = 10000, mean = 0, sd = sqrt(0.01))
W <- runif(n = 10000, min = 0, max = 1)

#Computations
#############

#1. Compute all the psi(W)
for (j in (1:10)) {
  assign(paste("psi", j, sep = "_"), psi(j, W))
}

#2. Compute the 10 true values of Y
for (j in (1:10)) {
  assign(paste("Y", j, sep = "_"), V)
  for (i in (1:j)) {
    assign(paste("Y", j, sep ="_"), 
           get(paste("Y", j, sep ="_")) 
           + (cvalue(i) * get(paste("psi", i, sep = "_"))))
  }
}

#3. Compute the estimated standard errors
se_1 <- coef(summary(lm(Y_1 ~ psi_1)))[1,"Std. Error"] / sqrt(lambda(1))
se_2 <- coef(summary(lm(Y_2 ~ psi_1 + psi_2)))[2,"Std. Error"] / sqrt(lambda(2))
se_3 <- coef(summary(lm(Y_3 ~ psi_1 + psi_2 + psi_3)))[3,"Std. Error"] / sqrt(lambda(3))
se_4 <- coef(summary(lm(Y_4 ~ psi_1 + psi_2 + psi_3 + psi_4)))[4,"Std. Error"] / sqrt(lambda(4))
se_5 <- coef(summary(lm(Y_5 ~ psi_1 + psi_2 + psi_3 + psi_4 + psi_5)))[5,"Std. Error"] / sqrt(lambda(5))
se_6 <- coef(summary(lm(Y_6 ~ psi_1 + psi_2 + psi_3 + psi_4 + psi_5 + psi_6)))[6,"Std. Error"] / sqrt(lambda(6))
se_7 <- coef(summary(lm(Y_7 ~ psi_1 + psi_2 + psi_3 + psi_4 + psi_5 + psi_6 + psi_7)))[7,"Std. Error"] / sqrt(lambda(7))
se_8 <- coef(summary(lm(Y_8 ~ psi_1 + psi_2 + psi_3 + psi_4 + psi_5 + psi_6 + psi_7 + psi_8)))[8,"Std. Error"] / sqrt(lambda(8))
se_9 <- coef(summary(lm(Y_9 ~ psi_1 + psi_2 + psi_3 + psi_4 + psi_5 + psi_6 + psi_7 + psi_8 + psi_9)))[9,"Std. Error"] / sqrt(lambda(9))
se_10 <- coef(summary(lm(Y_10 ~ psi_1 + psi_2 + psi_3 + psi_4 + psi_5 + psi_6 + psi_7 + psi_8 + psi_9 + psi_10)))[10,"Std. Error"] / sqrt(lambda(10))
se <- c(se_1, se_2, se_3, se_4, se_5, se_6, se_7, se_8, se_9, se_10)

#4. Compute the 10 true coefficients
coeffs <- vector(mode = "numeric", length = 0)
for (i in (1:10)) {
  assign(paste("coeff", i, sep = "_"), abs(g(i)))
  coeffs <- c(coeffs, get(paste("coeff", i, sep = "_")))
}

#5. Add with variance of .05 rather than .01
V <- rnorm(n = 10000, mean = 0, sd = sqrt(0.05))
for (j in (1:10)) {
  assign(paste("psi", j, sep = "_"), psi(j, W))
}

for (j in (1:10)) {
  assign(paste("Y", j, sep = "_"), V)
  for (i in (1:j)) {
    assign(paste("Y", j, sep ="_"), 
           get(paste("Y", j, sep ="_")) 
           + (cvalue(i) * get(paste("psi", i, sep = "_"))))
  }
}

se_1 <- coef(summary(lm(Y_1 ~ psi_1)))[1,"Std. Error"] / sqrt(lambda(1))
se_2 <- coef(summary(lm(Y_2 ~ psi_1 + psi_2)))[2,"Std. Error"] / sqrt(lambda(2))
se_3 <- coef(summary(lm(Y_3 ~ psi_1 + psi_2 + psi_3)))[3,"Std. Error"] / sqrt(lambda(3))
se_4 <- coef(summary(lm(Y_4 ~ psi_1 + psi_2 + psi_3 + psi_4)))[4,"Std. Error"] / sqrt(lambda(4))
se_5 <- coef(summary(lm(Y_5 ~ psi_1 + psi_2 + psi_3 + psi_4 + psi_5)))[5,"Std. Error"] / sqrt(lambda(5))
se_6 <- coef(summary(lm(Y_6 ~ psi_1 + psi_2 + psi_3 + psi_4 + psi_5 + psi_6)))[6,"Std. Error"] / sqrt(lambda(6))
se_7 <- coef(summary(lm(Y_7 ~ psi_1 + psi_2 + psi_3 + psi_4 + psi_5 + psi_6 + psi_7)))[7,"Std. Error"] / sqrt(lambda(7))
se_8 <- coef(summary(lm(Y_8 ~ psi_1 + psi_2 + psi_3 + psi_4 + psi_5 + psi_6 + psi_7 + psi_8)))[8,"Std. Error"] / sqrt(lambda(8))
se_9 <- coef(summary(lm(Y_9 ~ psi_1 + psi_2 + psi_3 + psi_4 + psi_5 + psi_6 + psi_7 + psi_8 + psi_9)))[9,"Std. Error"] / sqrt(lambda(9))
se_10 <- coef(summary(lm(Y_10 ~ psi_1 + psi_2 + psi_3 + psi_4 + psi_5 + psi_6 + psi_7 + psi_8 + psi_9 + psi_10)))[10,"Std. Error"] / sqrt(lambda(10))
se2 <- c(se_1, se_2, se_3, se_4, se_5, se_6, se_7, se_8, se_9, se_10)

#6. Plot
df <- data.frame(estimates = coeffs,
                 se = se,
                 se2 = se2,
                 xabs = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

ggplot(df, aes(xabs)) + 
  geom_line(aes(y = estimates), size = 1) + 
  geom_point(aes(y = estimates), size = 1.5) +
  geom_line(aes(y = se), linetype = "dashed" , size = 1) +
  geom_point(aes(y = se), size = 1.5) +
  geom_line(aes(y = se2), linetype = "dotted" , size = 1) +
  geom_point(aes(y = se2), size = 1.5) +
  scale_x_discrete(name="j", limits=c(1, 10)) +
  scale_y_continuous(name="Fourier Coefficients and Standard Deviations", limits=c(0, 0.5)) +
  theme_bw() + theme(text = element_text(size=16))



####################################################################################################
####################################################################################################
##                                    Introduction                                                ##
####################################################################################################
####################################################################################################

##########
#1. Setup#
##########

#Parametrize the simulation.
n <- 10000
mean_x <- 10
var_x <- 1
mean_w <- 15
var_w <- 1
endogeneity <- 0.7
relevance <- 0.9
cov_xw <- relevance*sqrt(var_x)*sqrt(var_w)

#Generate the variables needed to model endogeneity due to omitted variable.
set.seed(1)
nu_epsilon <- mvrnorm(n, c(0, 0), matrix(c(1, endogeneity, endogeneity, 1), 2, 2))
nu <- nu_epsilon[,1]
epsilon <- nu_epsilon[,2]
x_w <- mvrnorm(n, c(mean_x, mean_w), matrix(c(var_x, cov_xw, cov_xw, var_w), 2, 2))
x <- x_w[,1]
w <- x_w[,2]
x_obs <- x + nu

#Check that endogeneity is modeled properly.
Coefficients <- c(round(summary(lm(x_obs~epsilon))$coefficient[2,1], 3), 
                  round(summary(lm(x_obs~w))$coefficient[2,1], 3), 
                  round(summary(lm(w~epsilon))$coefficient[2,1], 3), 
                  round(summary(lm(x~epsilon))$coefficient[2,1], 3), 
                  round(summary(lm(x~nu))$coefficient[2,1], 3))
P_values <- c(round(summary(lm(x_obs~epsilon))$coefficient[2,4], 3), 
              round(summary(lm(x_obs~w))$coefficient[2,4], 3), 
              round(summary(lm(w~epsilon))$coefficient[2,4], 3), 
              round(summary(lm(x~epsilon))$coefficient[2,4], 3), 
              round(summary(lm(x~nu))$coefficient[2,4], 3))
check_table <- data.frame(Coefficients, P_values, 
                          row.names = c("Observed X on epsilon", "Observed X on W", "W on epsilon", 
                                        "True X on epsilon", "True X on nu"))
print(check_table)

#################################
#2. Linear relationship and TSLS#
#################################

#Illustrate the bias of OLS and of how TSLS overcomes it.
y <- (3 * x) + epsilon
xhat <- fitted.values(lm(x_obs ~ w))
Coefficients <- c(round(summary(lm(y ~ x))$coefficients[2,1], 2), 
                  round(summary(lm(y ~ x_obs))$coefficients[2,1], 2), 
                  round(summary(lm(y ~ xhat))$coefficients[2,1], 2))
P_values <- c(round(summary(lm(y ~ x))$coefficients[2,4], 2), 
              round(summary(lm(y ~ x_obs))$coefficients[2,4], 2), 
              round(summary(lm(y ~ xhat))$coefficients[2,4], 2))
#ivmodel(y,x,w) for corrected standard errors (<2e-16 as well)
linear_results <- data.frame(Coefficients, P_values, 
                             row.names = c("OLS with true X", "OLS with observed X", "TSLS"))
print(linear_results)

#Plot the results (not included in the final paper).
#OLS <- (summary(lm(y ~ x_obs))$coefficients[2,1] * x) + summary(lm(y ~ x_obs))$coefficients[1,1]
#TSLS <- (summary(lm(y ~ xhat))$coefficients[2,1] * x) + summary(lm(y ~ xhat))$coefficients[1,1]
#linearplot_data <- data.frame(y, x, OLS, TSLS)
#ggplot(linearplot_data, aes(x)) + 
#  geom_point(aes(y = y), alpha = 0.1) +
#  geom_line(aes(y = OLS), color = "red") +
#  geom_line(aes(y = TSLS), color = "blue") +
#  scale_x_continuous(name = "True X") +
#  scale_y_continuous(name = "Y")

#Check that this was not pure luck (not included in the final paper).
#betaOLS <- vector(mode = "numeric", length = 0)
#beta2SLS <- vector(mode = "numeric", length = 0)
#for (i in 1:1000){
#  mc_nu_epsilon <- mvrnorm(n, c(0, 0), matrix(c(1, endogeneity, endogeneity, 1), 2, 2))
#  mc_x_w <- mvrnorm(n, c(mean_x, mean_w), matrix(c(var_x, relevance, relevance, var_w), 2, 2))
#  mc_nu <- mc_nu_epsilon[,1]
#  mc_epsilon <- mc_nu_epsilon[,2]
#  mc_x <- mc_x_w[,1]
#  mc_w <- mc_x_w[,2]
#  
#  mc_x_obs <- mc_x + mc_nu
#  mc_y <- (3 * mc_x) + mc_epsilon
#  mc_xhat <- fitted.values(lm(mc_x_obs ~ mc_w))
#  
#  beta2SLS_i <- summary(lm(mc_y ~ mc_xhat))$coefficients[2,1]
#  beta2SLS <- c(beta2SLS, beta2SLS_i)
#}
#beta_data <- data.frame(beta2SLS)
#ggplot(beta_data, aes(x=beta2SLS)) + 
#  geom_histogram(binwidth=0.006) +
#  scale_x_continuous(name = "Estimated coefficients") +
#  scale_y_continuous(name = "Density")



####################################################################################################
####################################################################################################
##                                   Non-parametric estimation                                    ##
####################################################################################################
####################################################################################################

####################################################
#1. Define step by step the Legendre basis function#
####################################################

rm(list=ls())

#Procedure to normalize vector to unit interval.
normalize <- function(x) x <- (x - min(x))/(max(x) - min(x))

#Procedure to get Legendre basis function.
bfunc <- function(x,nf) {
  z <- 2*x - 1
  rz <- length(z)
  psi <- matrix(0, nrow=rz, ncol=10)
  psi[,1] <- replicate(rz,1)
  psi[,2] <- (sqrt(3))*z
  psi[,3] <- 0.5*(sqrt(5))*(3*(z^2) - 1)
  psi[,4] <- 0.5*(sqrt(7))*(5*(z^3) - 3*z)
  psi[,5] <- 0.125*3*(35*(z^4) - 30*(z^2) + 3)
  psi[,6] <- 0.125*(sqrt(11))*(63*(z^5) - 70*(z^3) + 15*z)
  psi[,7] <- 0.0625*(sqrt(13))*(231*(z^6) - 315*(z^4) + 105*(z^2) - 5)
  psi[,8] <- 0.0625*(sqrt(15))*(429*(z^7) - 693*(z^5) + 315*(z^3) - 35*z)
  psi[,9] <- (1/128)*(sqrt(17))*(6435*(z^8) - 12012*(z^6) + 6930*(z^4) - 1260*(z^2) + 35)
  psi[,10] <- (1/128)*(sqrt(19))*(12155*(z^9) - 25740*(z^7) + 18018*(z^5) - 4620*(z^3) + 315*z)
  psi <- psi[,1:(nf+1)]
  return(psi)
}

#Procedure to find bootstrap uniform confidence intervals.
bootstr <- function(bhat, invaht, delbar, data=wxy, degree=3, clvl=0.95, nboot=1000, xx) {
  clvl <- clvl*nboot
  maxt <- vector(mode="numeric", length=nboot)
  for (iboot in 1:nboot) {
    #Draw bootstrap sample
    vboot <- data[sample(1:nrow(data), nrow(data), replace=T),]
    
    #Select data
    wboot <- vboot[,1]
    xboot <- vboot[,2]
    yboot <- vboot[,3]
    xboot <- normalize(xboot)
    wboot <- normalize(wboot)
    
    #Estimate Fourier coefficients
    xxb <- bfunc(xboot, degree)
    wwb <- bfunc(wboot, degree)
    bhbt <- solve(t(wwb) %*% xxb) %*% t(wwb) %*% yboot
    
    uboot <- yboot - xxb %*% bhat
    delbt <- t(wwb) %*% uboot
    swig <- xx %*% invaht %*% delbt - delbar
    tstat <- abs(swig)
    
    maxt[iboot] <- max(tstat)
  }
  maxt <- sort(maxt)
  return(maxt[clvl])
}

#Procedure to get curve and CI for nonparametric IV with Legendre basis.
npiv_legendre <- function(y, x, w, deg=4, n_graph=100, n_boot=1000, conf_level=0.95, bootstrap=T, xrange=NULL) {
  #y --> (vector) explained variable
  #x --> (vector) endogenous regressor
  #w --> (vector) instrument for x
  #deg --> (scalar) degree of approximation
  #n_graph --> (scalar) number of points for the graph
  
  #Transform X and W to Unit Interval
  n <- length(y)
  x_norm <- normalize(x)
  w_norm <- normalize(w)
  
  #Create grid for graphs
  if (is.null(xrange)) {
    del <- (max(x_norm) - min(x_norm))/(n_graph+1)
    x0 <- min(x_norm) + del
    grid <- seq(from=x0, by=del, length.out=n_graph)
    gr0 <- grid*(max(x) - min(x)) + min(x)
  } else {
    xrange_1 <- (xrange[1]-min(x))/(max(x) - min(x))
    xrange_2 <- (xrange[2]-min(x))/(max(x) - min(x))
    grid <- seq(from=xrange_1, to=xrange_2, length.out=n_graph)
    gr0 <- seq(from=xrange[1], to=xrange[2], length.out=n_graph)
  }
  #Estimate Fourier coefficients
  xx <- bfunc(x_norm, deg)
  ww <- bfunc(w_norm, deg)
  xgrid <- bfunc(grid, deg)
  
  invaht <- solve(t(ww) %*% xx)
  bhat <- invaht %*% t(ww) %*% y
  yhat <- xgrid %*% bhat
  
  #Use bootstrap to get confidence band
  uhat <- y - xx %*% bhat
  delbar <- xx %*% invaht %*% colMeans(apply(ww, 2, function(x) uhat*x))   # here is where Horowitz made a mistake, using xgrid instead of xx
  wxy <- data.frame(w,x,y)
  
  #Bootstrap s.e., if specified. Return final data
  if (bootstrap==T) {
    cval <- bootstr(bhat, invaht, delbar, data=wxy, degree=deg, clvl=conf_level, nboot=n_boot, xx=xx)
    return(data.frame(x=gr0, yhat, ylow=yhat-cval, yhigh=yhat+cval))
  } else {
    return(data.frame(x=gr0, yhat))
  }
}

#Example from Horowitz
yz <- read.csv2(file = file.choose(), header=F)
yz <- apply(yz, 2, as.numeric)
z <- yz[,6] ; yz <- yz[z<=3,]
y <- yz[,5] ; x <- yz[,3] ; w <- yz[,7]
data <- npiv_legendre(y, x, w, deg=3, n_boot=1000)
ggplot(data, aes(x=x)) + geom_line(aes(y=yhat)) +
  geom_line(aes(y=ylow), linetype=2) +
  geom_line(aes(y=yhigh), linetype=2) +
  theme_bw() +
  scale_x_continuous("Class size", limits=c(20,40)) +
  scale_y_continuous("Test score", breaks = c(72,76,80,84), limits=c(70,86))



####################################################################################################
####################################################################################################
##                                Non-parametrics and precision                                   ##
####################################################################################################
####################################################################################################

####################################################################
#1. Apply Legendre on the previously defined nonlinear relationship#
####################################################################

n <- 2000
mean_x <- 10
var_x <- 1
mean_w <- 15
var_w <- 1
endogeneity <- 0.7
relevance <- 0.9
cov_xw <- relevance*sqrt(var_x)*sqrt(var_w)

set.seed(13)
nu_epsilon <- mvrnorm(n, c(0, 0), matrix(c(1, endogeneity, endogeneity, 1), 2, 2))
nu <- nu_epsilon[,1] ; epsilon <- nu_epsilon[,2]
x_w <- mvrnorm(n, c(mean_x, mean_w), matrix(c(var_x, relevance, relevance, var_w), 2, 2))
x <- x_w[,1] ; w <- x_w[,2]
x_obs <- x + nu
y <- 250 + ((x-10)^3) + epsilon

for (i in 1:4) {
  data <- npiv_legendre(y, x_obs, w, deg=i, n_boot=1000)
  p <- ggplot(data, aes(x=x)) +
    geom_line(data=data.frame(x=data$x, y=250+((data$x-10)^3)), aes(x,y), color="orange", size=2) +
    geom_point(data=data.frame(x_obs, y), aes(x_obs,y), alpha=0.1, color="darkblue") +
    geom_line(aes(y=yhat), size=2) +
    geom_line(aes(y=ylow), linetype=3, size=2) +
    geom_line(aes(y=yhigh), linetype=3, size=2) +
    xlim(5.5,15) + ylim(80, 450) +
    theme_bw() + theme(text = element_text(size=30))
  if (i>1) {
    p <- p + theme(axis.title.y = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks = element_blank(),
                   axis.text.x = element_text(color="white"))
  }
  assign(paste("plot", i, sep = ""), p)
}
plot1
plot2
plot3
plot4

#########################
#2. Heavy vs. thin tails#
#########################

### GAUSSIAN SCENARIO ###
# Run the estimation on 1000 simulations
for (z in 1:1000) {
  nu_epsilon <- mvrnorm(n, c(0, 0), matrix(c(1, endogeneity, endogeneity, 1), 2, 2))
  nu <- nu_epsilon[,1] ; epsilon <- nu_epsilon[,2]
  x_w <- mvrnorm(n, c(mean_x, mean_w), matrix(c(var_x, relevance, relevance, var_w), 2, 2))
  x <- x_w[,1] ; w <- x_w[,2]
  x_obs <- x + nu
  y <- 250 + ((x-10)^3) + epsilon
  if (z==1) {
    data <- cbind(sample = z, npiv_legendre(y, x_obs, w, deg=3, bootstrap=F, xrange=c(4,15)))
  } else {
    data <- rbind(data, cbind(sample = z, npiv_legendre(y, x_obs, w, deg=3, bootstrap=F, xrange=c(4,15))))
  }
}

# For every x, keep the mean fitted-y and the interval [ylow, yhigh] which contains 95% of fitted values
plotdf <- data %>% group_by(x) %>% summarize(yhigh = quantile(yhat, 0.95), ylow = quantile(yhat, 0.05), yhat=mean(yhat)) %>% ungroup()
ggplot(plotdf) + theme_bw() +
  geom_ribbon(aes(x, ymin=ylow, ymax=yhigh), alpha=0.3, color="gray") +
  geom_line(data=data.frame(x=data$x, y=250+((data$x-10)^3)), aes(x,y), color="orange", size=1) +
  geom_line(aes(x,yhat), size=1) + scale_y_continuous(limits=c(-150, 465)) + scale_x_continuous("Normally distributed x")


### STUDENT SCENARIO ###
# Run the estimation on 1000 simulations
for (z in 1:1000) {
  nu_epsilon <- mvrnorm(n, c(0, 0), matrix(c(1, endogeneity, endogeneity, 1), 2, 2))
  nu <- nu_epsilon[,1] ; epsilon <- nu_epsilon[,2]
  df <- 4
  s_norm <- mvrnorm(n, c(0,0), matrix(c(var_x, relevance, relevance, var_w)*(df-2), 2))
  x_w <- t(t(s_norm * sqrt(1/rchisq(n,df))) + c(mean_x, mean_w))
  x <- x_w[,1] ; w <- x_w[,2]
  x_obs <- x + nu
  y <- 250 + ((x-10)^3) + epsilon
  if (z==1) {
    data <- cbind(sample = z, npiv_legendre(y, x_obs, w, deg=3, bootstrap=F, xrange=c(4,15)))
  } else {
    data <- rbind(data, cbind(sample = z, npiv_legendre(y, x_obs, w, deg=3, bootstrap=F, xrange=c(4,15))))
  }
}

# For every x, keep the mean fitted-y and the interval [ylow, yhigh] which contains 95% of fitted values
plotdf <- data %>% group_by(x) %>% summarize(yhigh = quantile(yhat, 0.95), ylow = quantile(yhat, 0.05), yhat=mean(yhat)) %>% ungroup()
ggplot(plotdf) + theme_bw() +
  geom_ribbon(aes(x, ymin=ylow, ymax=yhigh), alpha=0.3, color="gray") +
  geom_line(data=data.frame(x=data$x, y=250+((data$x-10)^3)), aes(x,y), color="orange", size=1) +
  geom_line(aes(x,yhat), size=1) + scale_y_continuous(limits=c(-150, 465)) + scale_x_continuous("Student distributed x (df=4)")



####################################################################################################
####################################################################################################
##                                     Relaxing assumptions                                       ##
####################################################################################################
####################################################################################################

#################################################################
#1 Slightly endogenous and increasingly weak instrument and TSLS#
#################################################################
n <- 10000
endogeneity <- 0.7
mean_x <- 10
var_x <- 1
mean_w <- 15
var_w <- 1
bias <- 0.2

meancoeffs <- vector(mode = "numeric", length = 0)
meancsts <- vector(mode = "numeric", length = 0)

set.seed(1)
for (i in (seq(0.1, 0.9, 0.1))){
  TSLScoeffs <- vector(mode = "numeric", length = 0)
  TSLScsts <- vector(mode = "numeric", length = 0)
  for (j in (1:100)){
    relevance <- i
    cov_xw <- relevance*sqrt(var_x)*sqrt(var_w)
    
    
    nu_epsilon_x_w <- mvrnorm(n, c(0, 0, mean_x, mean_w), matrix(c(1, endogeneity, 0, 0, endogeneity, 1, 0, bias, 0, 0, 1,  
                                                                   cov_xw, 0, bias, cov_xw, 1), 4, 4))
    nu <- nu_epsilon_x_w[,1]
    epsilon <- nu_epsilon_x_w[,2]
    x <- nu_epsilon_x_w[,3]
    w <- nu_epsilon_x_w[,4]
    x_obs <- x + nu
    
    y <- (3 * x) + epsilon
    xhat <- fitted.values(lm(x_obs ~ w))
    
    TSLScoeffs <- c(TSLScoeffs, summary(lm(y ~ xhat))$coefficients[2,1]) 
    TSLScsts <- c(TSLScsts, summary(lm(y ~ xhat))$coefficients[1,1])
  }
  meancoeff <- mean(TSLScoeffs)
  meancoeffs <- c(meancoeffs, meancoeff)
  meancst <- mean(TSLScsts)
  meancsts <- c(meancsts, meancst)
}

x <- seq(5, 15, 0.1)
true <- 3 * x
data_weak <- data.frame(x, true)
for (i in (1:length(meancoeffs))){
  TSLS <- (meancoeffs[i] * x) + meancsts[i]
  data_weak <- data.frame(data_weak, TSLS)
}

ggplot(data_weak, aes(x)) + 
  geom_line(aes(y = TSLS.8), color = "grey90", size = 1.5) +
  geom_line(aes(y = TSLS.7), color = "grey80", size = 1.5) +
  geom_line(aes(y = TSLS.6), color = "grey70", size = 1.5) +
  geom_line(aes(y = TSLS.5), color = "grey60", size = 1.5) +
  geom_line(aes(y = TSLS.4), color = "grey50", size = 1.5) +
  geom_line(aes(y = TSLS.3), color = "grey40", size = 1.5) +
  geom_line(aes(y = TSLS.2), color = "grey30", size = 1.5) +
  geom_line(aes(y = TSLS.1), color = "grey20", size = 1.5) +
  geom_line(aes(y = TSLS), color = "grey10", size = 1.5) +
  geom_line(aes(y = true), color = "red", size = 1.5) +
  labs(x = "Generated data", y = "TSLS estimation") +
  theme_bw() + theme(text = element_text(size=17))

###########################################################################
#2 Slightly endogenous and increasingly weak instrument and Nonparametrics#
###########################################################################
Np_data <- data.frame(seq(1,100,1))

set.seed(1)
for (i in (seq(0.1, 0.9, 0.1))){
  x_vect <- vector(mode = "numeric", length = 0)
  y_vect <- vector(mode = "numeric", length = 0)
  relevance <- i
  cov_xw <- relevance*sqrt(var_x)*sqrt(var_w)

  nu_epsilon_x_w <- mvrnorm(n, c(0, 0, mean_x, mean_w), matrix(c(1, endogeneity, 0, 0, endogeneity, 1, 0, bias, 0, 0, 1,  
                                                                 cov_xw, 0, bias, cov_xw, 1), 4, 4))
  nu <- nu_epsilon_x_w[,1]
  epsilon <- nu_epsilon_x_w[,2]
  x <- nu_epsilon_x_w[,3]
  w <- nu_epsilon_x_w[,4]
  x_obs <- x + nu
  
  y <- 250 + ((x-10)^3) + epsilon
  Np_data <- data.frame(Np_data, npiv_legendre(y, x_obs, w, deg=4, n_boot=100))
}

x <- seq(4.12, 16, 0.12)
true <- 250 + ((x-10)^3)
Np_data <- data.frame(Np_data, x, true)

x_sequence <- seq(2, (dim(Np_data)[2] - 1), 4)
all_xs <- vector(mode = "numeric", length = 0)
for (i in x_sequence){
  all_xs <- c(all_xs, Np_data[,i])
}

all_unique_x <- unique(all_xs)
graph_data <- data.frame(all_unique_x)

y_sequence <- seq(3, (dim(Np_data)[2]), 4)
for (i in y_sequence){
  oldy <- Np_data[,i]
  corresponding_x <- Np_data[,i-1]
  newy <- rep(NA, length(all_unique_x))
  indices <- match(corresponding_x, all_unique_x)
  newy[indices] <- oldy
  graph_data <- data.frame(graph_data, newy)
}

ylow_sequence <- seq(4, (dim(Np_data)[2] - 3), 4)
for (i in ylow_sequence){
  oldylow <- Np_data[,i]
  corresponding_x <- Np_data[,i-2]
  newylow <- rep(NA, length(all_unique_x))
  indices <- match(corresponding_x, all_unique_x)
  newylow[indices] <- oldylow
  graph_data <- data.frame(graph_data, newylow)
}

yhigh_sequence <- seq(5, (dim(Np_data)[2] - 2), 4)
for (i in yhigh_sequence){
  oldyhigh <- Np_data[,i]
  corresponding_x <- Np_data[,i-3]
  newyhigh <- rep(NA, length(all_unique_x))
  indices <- match(corresponding_x, all_unique_x)
  newyhigh[indices] <- oldyhigh
  graph_data <- data.frame(graph_data, newyhigh)
}

truey <- with(graph_data, interp1(all_unique_x, newy.9, all_unique_x, "linear"))
y1 <- with(graph_data, interp1(all_unique_x, newy, all_unique_x, "linear"))
y2 <- with(graph_data, interp1(all_unique_x, newy.1, all_unique_x, "linear"))
y3 <- with(graph_data, interp1(all_unique_x, newy.2, all_unique_x, "linear"))
y4 <- with(graph_data, interp1(all_unique_x, newy.3, all_unique_x, "linear"))
y5 <- with(graph_data, interp1(all_unique_x, newy.4, all_unique_x, "linear"))
y6 <- with(graph_data, interp1(all_unique_x, newy.5, all_unique_x, "linear"))
y7 <- with(graph_data, interp1(all_unique_x, newy.6, all_unique_x, "linear"))
y8 <- with(graph_data, interp1(all_unique_x, newy.7, all_unique_x, "linear"))
y9 <- with(graph_data, interp1(all_unique_x, newy.8, all_unique_x, "linear"))
yl1 <- with(graph_data, interp1(all_unique_x, newylow, all_unique_x, "linear"))
yl2 <- with(graph_data, interp1(all_unique_x, newylow.1, all_unique_x, "linear"))
yl3 <- with(graph_data, interp1(all_unique_x, newylow.2, all_unique_x, "linear"))
yl4 <- with(graph_data, interp1(all_unique_x, newylow.3, all_unique_x, "linear"))
yl5 <- with(graph_data, interp1(all_unique_x, newylow.4, all_unique_x, "linear"))
yl6 <- with(graph_data, interp1(all_unique_x, newylow.5, all_unique_x, "linear"))
yl7 <- with(graph_data, interp1(all_unique_x, newylow.6, all_unique_x, "linear"))
yl8 <- with(graph_data, interp1(all_unique_x, newylow.7, all_unique_x, "linear"))
yl9 <- with(graph_data, interp1(all_unique_x, newylow.8, all_unique_x, "linear"))
yh1 <- with(graph_data, interp1(all_unique_x, newyhigh, all_unique_x, "linear"))
yh2 <- with(graph_data, interp1(all_unique_x, newyhigh.1, all_unique_x, "linear"))
yh3 <- with(graph_data, interp1(all_unique_x, newyhigh.2, all_unique_x, "linear"))
yh4 <- with(graph_data, interp1(all_unique_x, newyhigh.3, all_unique_x, "linear"))
yh5 <- with(graph_data, interp1(all_unique_x, newyhigh.4, all_unique_x, "linear"))
yh6 <- with(graph_data, interp1(all_unique_x, newyhigh.5, all_unique_x, "linear"))
yh7 <- with(graph_data, interp1(all_unique_x, newyhigh.6, all_unique_x, "linear"))
yh8 <- with(graph_data, interp1(all_unique_x, newyhigh.7, all_unique_x, "linear"))
yh9 <- with(graph_data, interp1(all_unique_x, newyhigh.8, all_unique_x, "linear"))

new_data <- data.frame(all_unique_x, truey, y1, y2, y3, y4, y5, y6, y7, y8, y9, yl1, yl2, yl3, 
                       yl4, yl5, yl6, yl7, yl8, yl9, yh1, yh2, yh3, yh4, yh5, yh6, yh7, yh8, yh9)

#Plot the results
ggplot(new_data, aes(x = all_unique_x)) +
  #geom_line(aes(y = y1), alpha = 0.6) +
  #geom_line(aes(y = yl1), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = yh1), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = y2), alpha = 0.6) +
  #geom_line(aes(y = yl2), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = yh2), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = y3), alpha = 0.6) +
  #geom_line(aes(y = yl3), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = yh3), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = y4), alpha = 0.6) +
  #geom_line(aes(y = yl4), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = yh4), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = y5), alpha = 0.6) +
  #geom_line(aes(y = yl5), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = yh5), alpha = 0.6, linetype = 2) +
  geom_line(aes(y = y6), alpha = 0.6, color = "grey30", size = 1.5) +
  #geom_line(aes(y = yl6), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = yh6), alpha = 0.6, linetype = 2) +
  geom_line(aes(y = y7), alpha = 0.6, color = "grey50", size = 1.5) +
  #geom_line(aes(y = yl7), alpha = 0.6, linetype = 2, color = "blue") +
  #geom_line(aes(y = yh7), alpha = 0.6, linetype = 2, color = "blue") +
  geom_line(aes(y = y8), alpha = 0.6, color = "grey70", size = 1.5) +
  #geom_line(aes(y = yl8), alpha = 0.6, linetype = 2, color = "green") +
  #geom_line(aes(y = yh8), alpha = 0.6, linetype = 2, color = "green") +
  geom_line(aes(y = y9), alpha = 0.6, color = "grey90", size = 1.5) +
  #geom_line(aes(y = yl9), alpha = 0.6, linetype = 2, color = "red") +
  #geom_line(aes(y = yh9), alpha = 0.6, linetype = 2, color = "red") +
  geom_line(aes(y = truey), color = "black", size = 1.5) +
  labs(x = "Simulated data", y = "Nonparametric estimation") +
  theme_bw() + theme(text = element_text(size=17))
  

ggplot(new_data, aes(x = all_unique_x)) +
  #geom_line(aes(y = y1), alpha = 0.6) +
  #geom_line(aes(y = yl1), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = yh1), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = y2), alpha = 0.6) +
  #geom_line(aes(y = yl2), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = yh2), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = y3), alpha = 0.6) +
  #geom_line(aes(y = yl3), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = yh3), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = y4), alpha = 0.6) +
  #geom_line(aes(y = yl4), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = yh4), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = y5), alpha = 0.6) +
  #geom_line(aes(y = yl5), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = yh5), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = y6), alpha = 0.6, color = "grey30") +
  #geom_line(aes(y = yl6), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = yh6), alpha = 0.6, linetype = 2) +
  geom_line(aes(y = y7), alpha = 0.6, color = "blue", size = 1.5) +
  geom_line(aes(y = yl7), alpha = 0.6, linetype = 2, color = "blue", size = 1.5) +
  geom_line(aes(y = yh7), alpha = 0.6, linetype = 2, color = "blue", size = 1.5) +
  geom_line(aes(y = y8), alpha = 0.6, color = "green", size = 1.5) +
  geom_line(aes(y = yl8), alpha = 0.6, linetype = 2, color = "green", size = 1.5) +
  geom_line(aes(y = yh8), alpha = 0.6, linetype = 2, color = "green", size = 1.5) +
  geom_line(aes(y = y9), alpha = 0.6, color = "red", size = 1.5) +
  geom_line(aes(y = yl9), alpha = 0.6, linetype = 2, color = "red", size = 1.5) +
  geom_line(aes(y = yh9), alpha = 0.6, linetype = 2, color = "red", size = 1.5) +
  geom_line(aes(y = truey), color = "black", size = 1.5) +
  labs(x = "Simulated data", y = "Nonparametric estimation") +
  theme_bw() + theme(text = element_text(size=17))

######################################
#3 Applied on a utility-like function#
######################################
endogeneity <- 0.7
mean_x <- 6
var_x <- 1
mean_w <- 11
var_w <- 1
bias <- 0
relevance <- .9

set.seed(1)
nu_epsilon_x_w2 <- mvrnorm(n, c(0, 0, mean_x, mean_w), matrix(c(1, endogeneity, 0, 0, endogeneity, 1, 0, bias, 0, 0, 1,  
                                                               relevance, 0, bias, relevance, 1), 4, 4))
nu2 <- nu_epsilon_x_w2[,1]
epsilon2 <- nu_epsilon_x_w2[,2]
x2 <- nu_epsilon_x_w2[,3]
w2 <- nu_epsilon_x_w2[,4]
x_obs2 <- x2 + nu2

y2 <- (15 * log(x2 + 1)) + epsilon2
data_utility <- npiv_legendre(y2, x_obs2, w2, deg=3, n_boot=100)

#How it looks like when the assumptions hold
ggplot(data_utility, aes(x=x)) + 
  geom_point(data=data.frame(x_obs2, y2), aes(x_obs2,y2), alpha=0.1) +
  geom_line(data=data.frame(x=data_utility$x, y=(15 * log(data_utility$x + 1))), aes(x,y), color="green", size = 1.5) +
  geom_line(aes(y=yhat), color = "blue", size = 1.5) +
  geom_line(aes(y=ylow), linetype=2, size = 1.5) +
  geom_line(aes(y=yhigh), linetype=2, size = 1.5) +
  labs(x = "Generated data", y = "Nonparametric estimation") +
  theme_bw() + theme(text = element_text(size=17))

#Relaxing the assumptions
Np_data <- data.frame(seq(1,100,1))

set.seed(1)
for (i in (seq(0.1, 0.9, 0.1))){
  x_vect <- vector(mode = "numeric", length = 0)
  y_vect <- vector(mode = "numeric", length = 0)
  relevance <- i
  cov_xw <- relevance*sqrt(var_x)*sqrt(var_w)
  
  nu_epsilon_x_w <- mvrnorm(n, c(0, 0, mean_x, mean_w), matrix(c(1, endogeneity, 0, 0, endogeneity, 1, 0, bias, 0, 0, 1,  
                                                                 cov_xw, 0, bias, cov_xw, 1), 4, 4))
  nu <- nu_epsilon_x_w[,1]
  epsilon <- nu_epsilon_x_w[,2]
  x <- nu_epsilon_x_w[,3]
  w <- nu_epsilon_x_w[,4]
  x_obs <- x + nu
  
  y <- (15 * log(x + 1)) + epsilon
  Np_data <- data.frame(Np_data, npiv_legendre(y, x_obs, w, deg=3, n_boot=100))
}

x <- seq(0.12, 12, 0.12)
true <- (15 * log(x + 1))
Np_data <- data.frame(Np_data, x, true)

x_sequence <- seq(2, (dim(Np_data)[2] - 1), 4)
all_xs <- vector(mode = "numeric", length = 0)
for (i in x_sequence){
  all_xs <- c(all_xs, Np_data[,i])
}

all_unique_x <- unique(all_xs)
graph_data <- data.frame(all_unique_x)

y_sequence <- seq(3, (dim(Np_data)[2]), 4)
for (i in y_sequence){
  oldy <- Np_data[,i]
  corresponding_x <- Np_data[,i-1]
  newy <- rep(NA, length(all_unique_x))
  indices <- match(corresponding_x, all_unique_x)
  newy[indices] <- oldy
  graph_data <- data.frame(graph_data, newy)
}

ylow_sequence <- seq(4, (dim(Np_data)[2] - 3), 4)
for (i in ylow_sequence){
  oldylow <- Np_data[,i]
  corresponding_x <- Np_data[,i-2]
  newylow <- rep(NA, length(all_unique_x))
  indices <- match(corresponding_x, all_unique_x)
  newylow[indices] <- oldylow
  graph_data <- data.frame(graph_data, newylow)
}

yhigh_sequence <- seq(5, (dim(Np_data)[2] - 2), 4)
for (i in yhigh_sequence){
  oldyhigh <- Np_data[,i]
  corresponding_x <- Np_data[,i-3]
  newyhigh <- rep(NA, length(all_unique_x))
  indices <- match(corresponding_x, all_unique_x)
  newyhigh[indices] <- oldyhigh
  graph_data <- data.frame(graph_data, newyhigh)
}

truey <- with(graph_data, interp1(all_unique_x, newy.9, all_unique_x, "linear"))
y1 <- with(graph_data, interp1(all_unique_x, newy, all_unique_x, "linear"))
y2 <- with(graph_data, interp1(all_unique_x, newy.1, all_unique_x, "linear"))
y3 <- with(graph_data, interp1(all_unique_x, newy.2, all_unique_x, "linear"))
y4 <- with(graph_data, interp1(all_unique_x, newy.3, all_unique_x, "linear"))
y5 <- with(graph_data, interp1(all_unique_x, newy.4, all_unique_x, "linear"))
y6 <- with(graph_data, interp1(all_unique_x, newy.5, all_unique_x, "linear"))
y7 <- with(graph_data, interp1(all_unique_x, newy.6, all_unique_x, "linear"))
y8 <- with(graph_data, interp1(all_unique_x, newy.7, all_unique_x, "linear"))
y9 <- with(graph_data, interp1(all_unique_x, newy.8, all_unique_x, "linear"))
yl1 <- with(graph_data, interp1(all_unique_x, newylow, all_unique_x, "linear"))
yl2 <- with(graph_data, interp1(all_unique_x, newylow.1, all_unique_x, "linear"))
yl3 <- with(graph_data, interp1(all_unique_x, newylow.2, all_unique_x, "linear"))
yl4 <- with(graph_data, interp1(all_unique_x, newylow.3, all_unique_x, "linear"))
yl5 <- with(graph_data, interp1(all_unique_x, newylow.4, all_unique_x, "linear"))
yl6 <- with(graph_data, interp1(all_unique_x, newylow.5, all_unique_x, "linear"))
yl7 <- with(graph_data, interp1(all_unique_x, newylow.6, all_unique_x, "linear"))
yl8 <- with(graph_data, interp1(all_unique_x, newylow.7, all_unique_x, "linear"))
yl9 <- with(graph_data, interp1(all_unique_x, newylow.8, all_unique_x, "linear"))
yh1 <- with(graph_data, interp1(all_unique_x, newyhigh, all_unique_x, "linear"))
yh2 <- with(graph_data, interp1(all_unique_x, newyhigh.1, all_unique_x, "linear"))
yh3 <- with(graph_data, interp1(all_unique_x, newyhigh.2, all_unique_x, "linear"))
yh4 <- with(graph_data, interp1(all_unique_x, newyhigh.3, all_unique_x, "linear"))
yh5 <- with(graph_data, interp1(all_unique_x, newyhigh.4, all_unique_x, "linear"))
yh6 <- with(graph_data, interp1(all_unique_x, newyhigh.5, all_unique_x, "linear"))
yh7 <- with(graph_data, interp1(all_unique_x, newyhigh.6, all_unique_x, "linear"))
yh8 <- with(graph_data, interp1(all_unique_x, newyhigh.7, all_unique_x, "linear"))
yh9 <- with(graph_data, interp1(all_unique_x, newyhigh.8, all_unique_x, "linear"))

new_data <- data.frame(all_unique_x, truey, y1, y2, y3, y4, y5, y6, y7, y8, y9, yl1, yl2, yl3, 
                       yl4, yl5, yl6, yl7, yl8, yl9, yh1, yh2, yh3, yh4, yh5, yh6, yh7, yh8, yh9)

#Plot the results
ggplot(new_data, aes(x = all_unique_x)) +
  #geom_line(aes(y = y1), alpha = 0.4) +
  #geom_line(aes(y = yl1), alpha = 0.4, linetype = 2) +
  #geom_line(aes(y = yh1), alpha = 0.4, linetype = 2) +
  #geom_line(aes(y = y2), alpha = 0.4) +
  #geom_line(aes(y = yl2), alpha = 0.4, linetype = 2) +
  #geom_line(aes(y = yh2), alpha = 0.4, linetype = 2) +
  #geom_line(aes(y = y3), alpha = 0.4) +
  #geom_line(aes(y = yl3), alpha = 0.4, linetype = 2) +
  #geom_line(aes(y = yh3), alpha = 0.4, linetype = 2) +
  #geom_line(aes(y = y4), alpha = 0.4) +
  #geom_line(aes(y = yl4), alpha = 0.4, linetype = 2) +
  #geom_line(aes(y = yh4), alpha = 0.4, linetype = 2) +
  geom_line(aes(y = y5), alpha = 0.6, color = "grey10", size = 1.5) +
  #geom_line(aes(y = yl5), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = yh5), alpha = 0.6, linetype = 2) +
  geom_line(aes(y = y6), alpha = 0.6, color = "grey30", size = 1.5) +
  #geom_line(aes(y = yl6), alpha = 0.6, linetype = 2) +
  #geom_line(aes(y = yh6), alpha = 0.6, linetype = 2) +
  geom_line(aes(y = y7), alpha = 0.6, color = "grey50", size = 1.5) +
  #geom_line(aes(y = yl7), alpha = 0.6, linetype = 2, color = "blue") +
  #geom_line(aes(y = yh7), alpha = 0.6, linetype = 2, color = "blue") +
  geom_line(aes(y = y8), alpha = 0.6, color = "grey70", size = 1.5) +
  #geom_line(aes(y = yl8), alpha = 0.6, linetype = 2, color = "green") +
  #geom_line(aes(y = yh8), alpha = 0.6, linetype = 2, color = "green") +
  geom_line(aes(y = y9), alpha = 0.6, color = "grey90", size = 1.5) +
  #geom_line(aes(y = yl9), alpha = 0.6, linetype = 2, color = "red") +
  #geom_line(aes(y = yh9), alpha = 0.6, linetype = 2, color = "red") +
  geom_line(aes(y = truey), color = "black", size = 1.5) +
  labs(x = "Generated data", y = "Nonparametric estimation") +
  scale_y_continuous(limits=c(0, 60)) +
  scale_x_continuous(limits=c(0, 12)) +
  theme_bw() + theme(text = element_text(size=17))

ggplot(new_data, aes(x = all_unique_x)) +
  #geom_line(aes(y = y1), alpha = 0.4) +
  #geom_line(aes(y = yl1), alpha = 0.4, linetype = 2) +
  #geom_line(aes(y = yh1), alpha = 0.4, linetype = 2) +
  #geom_line(aes(y = y2), alpha = 0.4) +
  #geom_line(aes(y = yl2), alpha = 0.4, linetype = 2) +
  #geom_line(aes(y = yh2), alpha = 0.4, linetype = 2) +
  #geom_line(aes(y = y3), alpha = 0.4) +
  #geom_line(aes(y = yl3), alpha = 0.4, linetype = 2) +
  #geom_line(aes(y = yh3), alpha = 0.4, linetype = 2) +
  #geom_line(aes(y = y4), alpha = 0.4) +
  #geom_line(aes(y = yl4), alpha = 0.4, linetype = 2) +
  #geom_line(aes(y = yh4), alpha = 0.4, linetype = 2) +
  #geom_line(aes(y = y5), alpha = 0.4) +
  #geom_line(aes(y = yl5), alpha = 0.4, linetype = 2) +
  #geom_line(aes(y = yh5), alpha = 0.4, linetype = 2) +
  #geom_line(aes(y = y6), alpha = 0.4) +
  #geom_line(aes(y = yl6), alpha = 0.4, linetype = 2) +
  #geom_line(aes(y = yh6), alpha = 0.4, linetype = 2) +
  geom_line(aes(y = y7), alpha = 0.6, color = "blue", size = 1.5) +
  geom_line(aes(y = yl7), alpha = 0.6, linetype = 2, color = "blue", size = 1.5) +
  geom_line(aes(y = yh7), alpha = 0.6, linetype = 2, color = "blue", size = 1.5) +
  geom_line(aes(y = y8), alpha = 0.6, color = "green", size = 1.5) +
  geom_line(aes(y = yl8), alpha = 0.6, linetype = 2, color = "green", size = 1.5) +
  geom_line(aes(y = yh8), alpha = 0.6, linetype = 2, color = "green", size = 1.5) +
  geom_line(aes(y = y9), alpha = 0.6, color = "red", size = 1.5) +
  geom_line(aes(y = yl9), alpha = 0.6, linetype = 2, color = "red", size = 1.5) +
  geom_line(aes(y = yh9), alpha = 0.6, linetype = 2, color = "red", size = 1.5) +
  geom_line(aes(y = truey), color = "black", size = 1.5) +
  scale_x_continuous(limits = c(0,12)) +
  labs(x = "Generated data", y = "Nonparametric estimation") +
  theme_bw() + theme(text = element_text(size=17))

####################################################################################################
####################################################################################################
##                                           APPENDIX                                             ##
####################################################################################################
####################################################################################################

####################################
#1. Nonlinear relationship and TSLS# (Not included in the final paper)
####################################

##Parametrize the simulation.
#n <- 10000
#mean_x <- 10
#var_x <- 1
#mean_w <- 15
#var_w <- 1
#endogeneity <- 0.7
#relevance <- 0.9
#cov_xw <- relevance*sqrt(var_x)*sqrt(var_w)
#
##Generate the variables needed to model endogeneity due to omitted variable.
#set.seed(1)
#nu_epsilon <- mvrnorm(n, c(0, 0), matrix(c(1, endogeneity, endogeneity, 1), 2, 2))
#nu <- nu_epsilon[,1]
#epsilon <- nu_epsilon[,2]
#x_w <- mvrnorm(n, c(mean_x, mean_w), matrix(c(var_x, cov_xw, cov_xw, var_w), 2, 2))
#x <- x_w[,1]
#w <- x_w[,2]
#x_obs <- x + nu
#
##Plot
#xhat <- fitted.values(lm(x_obs ~ w))
#y <- 250 + ((x-10)^3) + epsilon
#TSLS <- (summary(lm(y ~ xhat))$coefficients[2,1] * x) + summary(lm(y ~ xhat))$coefficients[1,1]
#nonlinear_data <- data.frame(x, y, TSLS)
#ggplot(nonlinear_data, aes(x)) +
#  geom_point(aes(y=y), alpha = 0.1) +
#  geom_line(aes(y = TSLS), color = "red") +
#  labs(x = "True X", y = "TSLS estimation")

##################################
#2. B-splines as a basis function#
##################################

#Define step by step the B-spline basis function
################################################

#The following function takes as arguments a number of knots (m + 1) and the degree k of the splines, 
#and returns the set of (m + 1 + 2k) knots with the interior knots evenly spaced on the unit interval 
#[0,1].
knot_vector <- function(mp1, k) {
  #Generate the interior knots
  knot_set <- 0
  nk1 <- mp1 - 1
  for (i in (1:nk1)) {
    lastval <- knot_set[length(knot_set)]
    newval <- lastval + (1 / nk1)
    knot_set <- c(knot_set, newval)
  }
  #Extend the lower bound
  lbound_ext <- -1/nk1
  for (i in (1:(k - 1))) {
    lastval <- lbound_ext[length(lbound_ext)]
    newval <- lastval - (1 / nk1)
    lbound_ext <- c(lbound_ext, newval)
  }
  lbound_ext <- rev(lbound_ext)
  #Extend the upper bound
  ubound_ext <- 1 + (1/nk1)
  for (i in (1:(k - 1))) {
    lastval <- ubound_ext[length(ubound_ext)]
    newval <- lastval + (1 / nk1)
    ubound_ext <- c(ubound_ext, newval)
  }
  
  knot_set <- c(lbound_ext, knot_set, ubound_ext)
  return(knot_set)
}

#Example of knot-set for uniform B-splines of degree 3 and with 4 interior knots:
knot_vector(4, 3)

#The following function takes as arguments the order of the function k and its set of knots generated 
#with the previously defined function and provides what corresponds to the term inside the brackets 
#in the equation for each p. 
spline_mat <- function(k, knot_set, del){
  nb_knots <- length(knot_set)
  nb_p <- nb_knots - k - 1
  smat <- matrix(rep(0, len=(nb_knots*nb_p)), nrow = nb_knots) 
  for (row in (1:nb_knots)){
    for (col in (1:nb_p)){
      col_max <- col + k + 1
      if ((row >= col) & (row <= col_max)){
        smat[row,col] <- 1 
      }
      sign <- 1
      if (smat[row,col] > 0){
        for (column in (col:col_max)){
          if (column != row){
            sign <- sign * (del / (column - row))
          }
        }
      }
      smat[row,col] <- smat[row,col] * sign
    }
  }
  return(smat)
}

#Finally, the following function takes as an input the x variable, the index of the spline to be 
#computed, the degree of the function, the set of knots and what is inside the brackets in the 
#formula. As an output, it returns Bp(x), that is, the pth B-spline.
pth_spline <- function(variable, p, k, knot_set, spline_mat){
  #Compute x minus the knot for each knot.
  xminusxi <- matrix(rep(0, len=(length(knot_set)*length(variable))), 
                     nrow = length(knot_set))
  for (i in (1:length(knot_set))) {
    for (j in (1:length(variable))) {
      xminusxi[i, j] <- variable[j] - knot_set[i]
    }
  }
  #Raise it to the power k.
  xminusxik <- abs(xminusxi)^k
  term <- matrix(rep(0, len=(length(knot_set)*length(variable))), 
                 nrow = length(knot_set))
  #Compute to product of what was define in the previous function and (x - xi)^k.
  for (i in (1:dim(term)[1])) {
    for (j in (1:dim(term)[2])) {
      if (xminusxi[i, j] > 0) {
        col <- spline_mat[,p]
        term[i,j] <- xminusxik[i,j]*col[i]
      }
    }
  }
  #Apply the sum operator.
  pth_spline <- vector(mode = "numeric", length = 0)
  for (i in (1:(dim(term)[2]))){
    newval <- sum(term[,i])
    pth_spline <- c(pth_spline, newval)
  }
  return(pth_spline)
}

#The following function combines all what was set up above in one function that takes as an intput a 
#variable x, the number of knots (m + 1) and the degree k of the function. It returns as an output a 
#matrix of dimension (N * (m + k)) containing each Bp(xi); i = 1,...,N.
B_spline <- function(variable, mp1, k) {
  
  knot_set <- knot_vector(mp1, k)
  #smat <- spline_mat(k,knot_set, ((mp1 - 1) / 2)) #Uncomment this line and comment the following to get Horowitz's estimation instead
  smat <- spline_mat(k,knot_set, 1) 
  nb_p <- length(knot_set) - k - 1
  
  B_spline <- matrix(rep(0, length(variable)*nb_p), 
                     nrow = length(variable))
  for (p in (1:nb_p)) {
    B_spline[,p] <- pth_spline(variable, p, k, knot_set, smat)
  }
  return(B_spline)
}

#For instance, Horowitz uses B-splines as a basis function to estimate nonparametrically an Engel 
#curve as follows.
data1 <- read.table(file = file.choose(), header = FALSE, sep = "")
data <- data1
y <- data[,1]
x <- data[,2]
w <- data[,3]
#Convert variables to unit interval
delx <- max(x) - min(x)
minx <- min(x)
x <- (x - minx)/delx
w <- (w - min(w))/(max(w) - min(w))
w <- (w - min(w))/(max(w) - min(w))
#Apply the basis function
xx <- 5 * B_spline(x, 4, 3) 
ww <- 5 * B_spline(w, 4, 3) 

an <- 10^(-9) #Stabilization Parameter
denominator <- (crossprod(xx, ww)%*%crossprod(ww,xx))/(length(x)^2)
numerator <- (crossprod(xx,ww)%*%t(ww))/(length(x)^2)
nb_col_den = dim(denominator)[1]
denominator <- denominator + (an * diag(nb_col_den)) 
Ghat <- solve(denominator) %*% numerator %*% y
yhat <- xx %*% Ghat
xx <- (delx * x) + minx
results <- data.frame(yhat,xx)
results <- results[order(yhat),]

#Plot the curve
ggplot(results, aes(results$xx)) + 
  geom_line(aes(y = yhat), size = 1.5) + 
  scale_x_continuous(name="Logarithm of Total Expenditures", limits=c(3.5, 8)) +
  scale_y_continuous(name="Estimated Expenditure Share", limits=c(0, 0.3)) +
  theme_bw() + theme(text = element_text(size=17))
