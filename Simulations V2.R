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

#install.packages("ggplot2", "MASS")
library(ggplot2)
library(MASS)

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

#Generate the variables needed to model endogeneity due to omitted variable.
set.seed(1)
nu_epsilon <- mvrnorm(n, c(0, 0), matrix(c(1, endogeneity, endogeneity, 1), 2, 2))
nu <- nu_epsilon[,1]
epsilon <- nu_epsilon[,2]
x_w <- mvrnorm(n, c(mean_x, mean_w), matrix(c(var_x, relevance, relevance, var_w), 2, 2))
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
linear_results <- data.frame(Coefficients, P_values, 
                             row.names = c("OLS with true X", "OLS with observed X", "TSLS"))
print(linear_results)

#Plot the results.
OLS <- (summary(lm(y ~ x_obs))$coefficients[2,1] * x) + summary(lm(y ~ x_obs))$coefficients[1,1]
TSLS <- (summary(lm(y ~ xhat))$coefficients[2,1] * x) + summary(lm(y ~ xhat))$coefficients[1,1]
linearplot_data <- data.frame(y, x, OLS, TSLS)
ggplot(linearplot_data, aes(x)) + 
  geom_point(aes(y = y), alpha = 0.1) +
  geom_line(aes(y = OLS), color = "red") +
  geom_line(aes(y = TSLS), color = "blue") +
  scale_x_continuous(name = "True X") +
  scale_y_continuous(name = "Y")

#Check that this was not pure luck.
betaOLS <- vector(mode = "numeric", length = 0)
beta2SLS <- vector(mode = "numeric", length = 0)
for (i in 1:1000){
  mc_nu_epsilon <- mvrnorm(n, c(0, 0), matrix(c(1, endogeneity, endogeneity, 1), 2, 2))
  mc_x_w <- mvrnorm(n, c(mean_x, mean_w), matrix(c(var_x, relevance, relevance, var_w), 2, 2))
  mc_nu <- mc_nu_epsilon[,1]
  mc_epsilon <- mc_nu_epsilon[,2]
  mc_x <- mc_x_w[,1]
  mc_w <- mc_x_w[,2]
  
  mc_x_obs <- mc_x + mc_nu
  mc_y <- (3 * mc_x) + mc_epsilon
  mc_xhat <- fitted.values(lm(mc_x_obs ~ mc_w))
  
  beta2SLS_i <- summary(lm(mc_y ~ mc_xhat))$coefficients[2,1]
  beta2SLS <- c(beta2SLS, beta2SLS_i)
}
beta_data <- data.frame(beta2SLS)
ggplot(beta_data, aes(x=beta2SLS)) + 
  geom_histogram(binwidth=0.006) +
  scale_x_continuous(name = "Estimated coefficients") +
  scale_y_continuous(name = "Density")

####################################
#3. Nonlinear relationship and TSLS#
####################################

y <- 250 + ((x-10)^3) + epsilon
TSLS <- (summary(lm(y ~ xhat))$coefficients[2,1] * x) + summary(lm(y ~ xhat))$coefficients[1,1]
nonlinear_data <- data.frame(x, y, TSLS)
ggplot(nonlinear_data, aes(x)) +
  geom_point(aes(y=y), alpha = 0.1) +
  geom_line(aes(y = TSLS), color = "red") +
  scale_x_continuous(name = "True X") +
  scale_y_continuous(name = "TSLS estimation")

####################################################################################################
####################################################################################################
##                                2. Non-parametric estimation                                    ##
####################################################################################################
####################################################################################################

#########################
#2.1 B-spline estimation#
#########################


#2.1.1 Define step by step the B-spline basis function
######################################################

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
  smat <- spline_mat(k,knot_set, 1)#((mp1 - 1) / 2))
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
  geom_line(aes(y = yhat)) + 
  scale_x_continuous(name="Logarithm of Total Expenditures", limits=c(3.5, 8)) +
  scale_y_continuous(name="Estimated Expenditure Share", limits=c(0, 0.3))

#2.1.2 Apply B-splines on the previously defined nonlinear relationship
#######################################################################
set.seed(1)
nu_epsilon <- mvrnorm(n, c(0, 0), matrix(c(1, endogeneity, endogeneity, 1), 2, 2))
nu <- nu_epsilon[,1]
epsilon <- nu_epsilon[,2]
x_w <- mvrnorm(n, c(mean_x, mean_w), matrix(c(var_x, relevance, relevance, var_w), 2, 2))
x <- x_w[,1]
w <- x_w[,2]
x_obs <- x + nu
xhat <- fitted.values(lm(x_obs ~ w))

y <- 250 + ((x-10)^3) + epsilon
delx <- max(x_obs) - min(x_obs)
minx <- min(x_obs)
norm_x <- (x_obs - minx)/delx
norm_w <- (w - min(w))/(max(w) - min(w))
xx <- 5 * B_spline(norm_x, 3, 2)
ww <- 5 * B_spline(norm_w, 3, 2)
an <- 10^(-9)

denominator <- (crossprod(xx, ww)%*%crossprod(ww,xx))/(length(x)^2)
numerator <- (crossprod(xx,ww)%*%t(ww))/(length(x)^2)
nb_col_den = dim(denominator)[1]
denominator <- denominator + (an * diag(nb_col_den)) 
Ghat <- solve(denominator) %*% numerator %*% y
yhat <- xx %*% Ghat

#Estimate TSLS
xhat <- fitted.values(lm(x_obs ~ w))
TSLS <- (summary(lm(y ~ xhat))$coefficients[2,1] * x_obs) + summary(lm(y ~ xhat))$coefficients[1,1]
trueyhat <- 250 + (((x_obs-10)^3) / 2) + epsilon
spline_data <- data.frame(yhat, x_obs, trueyhat, TSLS)

#Plot the results
ggplot(spline_data, aes(x = x_obs)) +
  geom_point(aes(y = trueyhat), color = "black", alpha = 0.1) +
  geom_line(aes(y = TSLS), color = "red", size = 1.1) +
  geom_line(aes(y = yhat), color = "blue", size = 1.1) +
  scale_x_continuous(name = "True X") +
  scale_y_continuous(name = "B-spline estimation") 

#2.1.3 Vary the degree of the function to show that if we don't know the true relationship it could 
#be misleading

xx <- 5 * B_spline(norm_x, 3, 3)
ww <- 5 * B_spline(norm_w, 3, 3)
an <- 10^(-9)
denominator <- (crossprod(xx, ww)%*%crossprod(ww,xx))/(length(x)^2)
numerator <- (crossprod(xx,ww)%*%t(ww))/(length(x)^2)
nb_col_den = dim(denominator)[1]
denominator <- denominator + (an * diag(nb_col_den)) 
Ghat <- solve(denominator) %*% numerator %*% y
yhat_deg3 <- xx %*% Ghat

spline_data2 <- data.frame(spline_data, yhat_deg3, yhat_deg4)
ggplot(spline_data2, aes(x = x_obs)) +
  geom_point(aes(y = trueyhat), color = "black", alpha = 0.1) +
  geom_line(aes(y = yhat), color = "red") +
  geom_line(aes(y = yhat_deg3), color = "blue") +
  scale_x_continuous(name = "True X") +
  scale_y_continuous(name = "B-spline estimations")

##################################
#2.2 Nonparametrics and precision#
##################################

#2.2.1 Illustration of parametric's "false precision"
#####################################################

spline_data <- data.frame(spline_data)
#ADD CONFIDENCE INTERVALS AND PLOT IT)
ggplot(spline_data, aes(x = x_obs)) +
  geom_point(aes(y = trueyhat), color = "black", alpha = 0.1) +
  geom_line(aes(y = TSLS), color = "red", size = 1.1) +
  geom_line(aes(y = yhat), color = "blue", size = 1.1) +
  scale_x_continuous(name = "True X") +
  scale_y_continuous(name = "B-spline estimation") 

#2.2.2 See how B-splines (shape + CI) react with data that are increasingly noisy (reduce the 
#relative variance of x compared to the one of the error term as we did last time we met)

#2.1.3 Vary the steepness of the tails of the function with both normal and cauchy distributed x to 
#illustrate the implications of having few data points for some range of x

#2.2.4 maybe we could also plot the b-splines of different orders each with some confidence 
#intervals, because if the b-splines that are not accurate have relatively much larger CI than the 
#accurate one (degree 2 with 3 knots), we could put the conclusions of 2.1.3 into perspective.

####################################################################################################
####################################################################################################
##                                   3. Relaxing assumptions                                      ##
####################################################################################################
####################################################################################################

##################################################
#3.1 Weak and increasingly endogenous instruments#
##################################################

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
    
    nu_epsilon_x_w <- mvrnorm(n, c(0, 0, mean_x, mean_w), matrix(c(1, endogeneity, 0, 0, endogeneity, 1, 0, bias, 0, 0, 1,  
                                                                   relevance, 0, bias, relevance, 1), 4, 4))
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
  geom_line(aes(y = true), color = "red") +
  geom_line(aes(y = TSLS), color = "grey10") +
  geom_line(aes(y = TSLS.1), color = "grey20") +
  geom_line(aes(y = TSLS.2), color = "grey30") +
  geom_line(aes(y = TSLS.3), color = "grey40") +
  geom_line(aes(y = TSLS.4), color = "grey50") +
  geom_line(aes(y = TSLS.5), color = "grey60") +
  geom_line(aes(y = TSLS.6), color = "grey70") +
  geom_line(aes(y = TSLS.7), color = "grey80") +
  geom_line(aes(y = TSLS.8), color = "grey90")

#3.2 Do the same with B-splines and see how the shape and the confidence intervals are affected

#3.3 Introduce Legendre polynomials (a bit like in 2.1.1)  and do the same with Legendre polynomials 
#and maybe a different nonlinear relationship such as y = 0.5 * sqrt(x) for instance 
