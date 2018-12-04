##########################################################
#                                                        #
#                  Generate Figure 2                     #
#                                                        #
##########################################################

rm(list=ls())

##########################################################
#         Define the values of the parameters            #
#                and prepare the data                    #
##########################################################

#Parametrization
rparm <- 10^(-9) #Stabilization Parameter
nbase <- 4
ksp <- 3 #Order of the spline
nk0 <- 4 #No. of spline knots in [-1,1]

#Data import
data <- read.table(file = "Figure_2_data.txt", header = FALSE, sep = "")
y <- data[,1]
x <- data[,2]
w <- data[,3]
n <- length(x)

#Transform X and W to unit interval
delx <- max(x) - min(x)
minx <- min(x)
x <- (x - minx)/delx
w <- (w - min(w))/(max(w) - min(w))

##########################################################
#     Define functions used in the main program          #
##########################################################

#Generate n knot points [0,1] for B-splines of order k
#Create a vector of numbers from 0 to 1 with increment 1/nk1, 
#then k numbers below 0 and above 1 with increment 1/nk1
ksplines <- function(nk1, k) {
  
  kn0 <- 0
  for (i in (1:nk1)) {
    lastval <- kn0[length(kn0)]
    newval <- lastval + (1 / nk1)
    kn0 <- c(kn0, newval)
  }
  
  temp0 <- -1/nk1
  for (i in (1:(k - 1))) {
    lastval <- temp0[length(temp0)]
    newval <- lastval - (1 / nk1)
    temp0 <- c(temp0, newval)
  }
  temp0 <- rev(temp0)
  
  temp1 <- 1 + (1/nk1)
  for (i in (1:(k - 1))) {
    lastval <- temp1[length(temp1)]
    newval <- lastval + (1 / nk1)
    temp1 <- c(temp1, newval)
  }
  
  mylist <- c(temp0, kn0, temp1)
  return(mylist)
  
}

#Generate matrix identifying the non-zero terms of the spline summand
#and the product coefficient for each spline component. The p^{th} 
#column of *smat* corresponds to the p^{th} B-spline.
spmat <- function(k, xkt, del) {
  nb_rows <- length(xkt) #number of elements in the vector generated with kspline
  nb_cols <- nb_rows - k - 1
  smat <- matrix( rep( 0, len=(nb_rows*nb_cols)), nrow = nb_rows) #0-matrix (nk rows and np col)
  
  #Compute the product that gives the sign and factor of the term
  for (row_i in (1:nb_rows)){
    
    for (col_j in (1:nb_cols)){
      #replace the elements of the 0 matrix by 1 if:
      #its line >= column & line =< column + k + 1
      col_max <- col_j + k + 1
      if ((row_i >= col_j) & (row_i<= col_max)){
        smat[row_i,col_j] <- 1 
      }
      
      sign <- 1
      
      if (smat[row_i,col_j] > 0){
        for (column in (col_j:col_max)){
          if (column != row_i){
            sign <- sign * (del / (column - row_i))
          }
        }
      }
      
      smat[row_i,col_j] <- smat[row_i,col_j] * sign
      
    }
  }
  return(smat)
}

#Evaluate the p^{th} B-spline at the values of x
b_spline <- function(variable, p, k, xkt, smat){
  mymat <- matrix(rep( 0, len=(length(xkt)*length(variable))), nrow = length(xkt))
  for (i in (1:length(xkt))) {
    for (j in (1:length(variable))) {
      mymat[i, j] <- variable[j] - xkt[i]
    }
  }
  
  mymat2 <- abs(mymat)^k

  term <- matrix(rep( 0, len=(length(xkt)*length(variable))), nrow = length(xkt))
  for (i in (1:dim(mymat)[1])) {
    for (j in (1:dim(mymat)[2])) {
      if (mymat[i, j] > 0) {
        col <- smat[,p]
        term[i,j] <- mymat2[i,j]*col[i]
      }
    }
  }
  
  #sum of each column of term
  vect <- vector(mode = "numeric", length = 0)
  for (i in (1:(dim(term)[2]))){
    newval <- sum(term[,i])
    vect <- c(vect, newval)
  }
  return(vect)
}

#Evaluate the spline basis at each value of X or W
spbase <- function(variable, nk, k, n) {
  nk1 <- nk - 1
  del <- nk1 / 2
  xkt <- ksplines(nk1, k)
  smat <- spmat(k,xkt,del)
  nkk <- length(xkt)
  np <- nkk - k - 1
  
  #Evaluate spline basis at each x value
  z <- matrix(rep(0, n*np), nrow = n)
  for (ip in (1:np)) {
    z[,ip] <- b_spline(variable, ip, k, xkt, smat)
  }
  return(z)
}

##########################################################
#            Estimate Fourier coefficients               #
##########################################################
xx <- 5 * spbase(x, nk0, ksp, n)
ww <- 5 * spbase(w, nk0, ksp, n)

qmat <- (crossprod(xx, ww)%*%crossprod(ww,xx))/(n^2)
nmat <- (crossprod(xx,ww)%*%t(ww))/(n^2)
nb_col_qmat = dim(qmat)[1]
qmat <- qmat + (rparm * diag(nb_col_qmat)) 
bhat <- solve(qmat) %*% nmat %*% y

yhat <- xx %*% bhat
xx <- (delx * x) + minx

results <- data.frame(yhat,xx)
results <- results[order(yhat),]

library(ggplot2)
ggplot(results, aes(results$xx)) + 
  geom_line(aes(y = yhat)) + 
  scale_x_continuous(name="Logarithm of Total Expenditures", limits=c(3.5, 8)) +
  scale_y_continuous(name="Estimated Expenditure Share", limits=c(0, 0.3)) 
  # + ggtitle("Estimated Engel curve for food")
