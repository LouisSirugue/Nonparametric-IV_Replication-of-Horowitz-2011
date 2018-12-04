####################################################################
#                                                                  #
#                   The ill-posed inverse problem                  #
#                                                                  #
####################################################################

rm(list = ls())
set.seed(1)

####################################################################
############################# SETUP ################################
####################################################################

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
V <- rnorm(n = 10000, mean = 0, sd = sqrt(0.01))
W <- runif(n = 10000, min = 0, max = 1)

####################################################################
######################### COMPUTATIONS #############################
####################################################################

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

####################################################################
############################## PLOT ################################
####################################################################

library(ggplot2)
df <- data.frame(estimates = coeffs,
                 se = se,
                 se2 = se2,
                 xabs = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

ggplot(df, aes(xabs)) + 
  geom_line(aes(y = estimates), size = 1) + 
  geom_point(aes(y = estimates), size = 2.5) +
  geom_line(aes(y = se), linetype = "dashed" , size = 1) +
  geom_point(aes(y = se), size = 2.5) +
  geom_line(aes(y = se2), linetype = "dotted" , size = 1) +
  geom_point(aes(y = se2), size = 2.5) +
  scale_x_discrete(name="j", limits=c(1, 10)) +
  scale_y_continuous(name="Fourier Coefficients and Standard Deviations", limits=c(0, 0.5)) 
  # + ggtitle("The ill-posed inverse problem: an illustration")
  