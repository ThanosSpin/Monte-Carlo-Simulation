
####Clear Environment####
rm(list = ls())

####Call Libraries####

library(doParallel)##Parallel Processing
library(psych)#describe function
library(knitr)##Create .md, .html, and .pdf files
library(markdown)##Create .md, .html, and .pdf files


#### Parallel procesing for running performance ####

# Set number of cores for parallel processing
no_cores <- detectCores() - 2 

# Register number of cores for parallel processing
registerDoParallel(no_cores)

# Validate number of cores for parallel processing 
getDoParWorkers()

#### Monte Carlo in 2d - Estimation the area of circle ####

# MC function in 2 dimensions ([0, d] = interval of the Uniform Distribution, n = number of samples) function
mc.2d <-function(x, y, d, n)
  {  
  accept <- 0 # set counter to zero

  data <- c(x, y) # starting vector

  accepts <- matrix(data, 1, 2) # set matrix for the sampling data
  
    for (i in 2:n) {
      y <- 2*runif(1, 0, d) - 1
      x <- 2*runif(1, 0 ,d) - 1
      
      if((x^2 + y^2) <= 1) { 
    accepts <- rbind(accepts, c(x,y))
    accept <- accept + 1
  }
}
alpha <- accept/n # probability to be inside the circle

par(mfrow=c(1, 2)) # set 2 graphs per window

plot(density(accepts[ , 1]), main = "") # probability density

index <- seq(-1, 1, 0.1)
lines(index, 0.5*index/index, col = "green") # add green line to the mean of Uniform distribution, interval: [-1, 1]

plot(accepts, xlab = "x accepts", ylab = "y accepts") # accepted points (points insode the circle)

return(list(alpha, accepts)) # returns the probability to be inside the circle and the accepted points
}

#### Checking the system time ####
dev.off() # close the active windows

set.seed(15) # seed for regeneration
system.time(mc.2d(0, 0, 1, 1000)) # starting point: (x,y) = (0, 0), [0, 1]: Uniform Distribution interval, n = 1000 samples

set.seed(16)# seed for regeneration
system.time(mc.2d(0, 0, 1, 10000)) # starting point: (x,y) = (0, 0), [0, 1]: Uniform Distribution interval, n = 10000 samples

set.seed(19)# seed for regeneration
system.time(mc.2d(0, 0, 1, 50000)) # starting point: (x,y) = (0, 0), [0, 1]: Uniform Distribution interval, n = 50000 samples

set.seed(20)# seed for regeneration
system.time(mc.2d(0, 0, 1, 100000)) # starting point: (x,y) = (0, 0), [0, 1]: Uniform Distribution interval, n = 100000 samples

set.seed(21)# seed for regeneration
system.time(mc.2d(0, 0, 1, 150000)) # starting point: (x,y) = (0, 0), [0, 1]: Uniform Distribution interval, n = 150000 samples



#### calculate the mean of the generated samples (MC acceptance ratio - probability) to verify where the chain converges ####
mc_meanx_2d <- function(n, d, m) { # n = samples, [0 , d] : Uniform distribution interval, m - iterations
  
  y <- 2*runif(1, 0, d) - 1
  z <- 2*runif(1, 0, d) - 1
  x <- c()
  acc <- c()  
  for (i in 1:m) {
    x[i] <- mc.2k(y, z, d, n)[[1]] #Insert the data frame to a vector
    acc <-  mc.2k(y, z, d, n)[[2]][ , 1]*mc.2d(y, z, d, n)[[2]][ , 2] # joint Uniform distribution f(x)*f(y)
    }
  z <- mean(x) #finding the mean of the acceptance rate - probability to be inside the circle 
  return(list(x, acc)) # returns the acceptance rates of the process 
}


set.seed(14) 
start.time <- Sys.time() # starting time
system.time(z <- mc_meanx_2d(10000, 1, 100))
end.time <- Sys.time() # ending time
(time.taken <- end.time - start.time) # running time

# The mean value of the 100 independent samples from the Uniform distribution is the expected value
# of the distribution
# (f(x1) + ... + f(x100))/100 =.. E[f(x)] (Law of Large Numbers - LLN)

# We calculate the expected value of the 100 independent samples from the uniform distribution
# in 2 dimensions

dev.off() # close graphs per window
mean(z[[2]]) # The expected value of the Uniform distribution, true mean of U ~[ -1, 1] : μ = (-1 + 1)/2 = 0


index <- seq(-1, 1, 0.1)

plot(density(z[[2]]), main ="") # plot probability density of z - vector f(x,y)
abline(v = mean(z[[2]]), col = "red", lwd = 1, lty = 3) # add the line of mean value
text(mean(z[[2]]), 0.1 , round(mean(z[[2]]), 3)) # add the mean value as text
# lines(index, 0.5*index/index, col = "green") # straight line at the mean value of the Unifrom distribution


plot(density(z[[1]]), main ="") # plot probability density of z - vector (MC acceptance ratio)
abline(v = mean(z[[1]]), col = "red", lwd = 1, lty = 3) # add the line of mean value
text(mean(z[[1]]), 0.1 , round(mean(z[[1]]), 3)) # add the mean value as text


mean(z[[1]])# The expected value of the MC acceptance ratio - probability for a point to be inside the circle
dev.off()


hist(z[[1]], probability = TRUE, ylim = c(0, 120))# histogram of probabilities z (f(x1),..., f(x100))
lines(density(z[[1]]),col="red") # add line of probability density

sd(z[[2]]) # standard deviation of joint distribution of accepted values

# We know that the estimation of the integrate is |B| = prob_estimated*(2^k) where k = 2 are the dimensions.
# So the |B| estimation is mean(z)*(2^2)

#### Estimate the integrate of B (2 dimensions) ####
(B_estimated <- mean(z[[1]])*(2^2))


# We know from the calculation of the integrate with polar coordinates (2 dimensions) 
# that that the area of probability density function upon 2d sphere (B_true) is pi/2*(2)
(B_true <- pi)

# We estimated with Monte Carlo sampling from Uniform Distribution 
# that the integrate of B is mean_estimated*(2^2) = mean(z)*4 = 3.141796 close to the true volume of 2d sphere (pi*R^2, 
# pi for R = 1)
# Correspondingly, the estimated area of 2 dimension's sphere (B_estimated): is mean(z)*2^2 
B_estimated == 2^2*mean(z[[1]])

# The error of the integrate estimation is B_estimated - B_true 
(estimated_error <- abs(B_true - B_estimated))
round(estimated_error, 4) # round value in four decimals

# The relative error of estimation of area B  is 0.02%
(estimated_error_perc <- round(abs(1 - (B_estimated/B_true)), 4)*100)

# It is a precise algorithm (estimation error very low < 1%)

#### Monte Carlo in k - dimensions (Volume estimation of hyper-sphere) ####

# MC function in k - dimensions 
mc.d <-function(n, d, k) # k = dimensions, U ~ [0, d], n = samples 
{  
  accept <- 0 # set counter to zero
  
  x <- (2*runif(k, 0, d) - d) # U ~ [-1,1]
  
  data <- c(x) # assign data to a vector
  
  accepts <- matrix(data, 1, k) # matrix of 1 row and k - columns
  
  for (i in 1:n) {
    x <- t(matrix(2*runif(k, 0, d) - d)) # simulate 2 values from Uniform distribution
    if(sum(x^2) <= 1) { # if the candidate point is inside the circle
      accepts <- rbind(accepts, c(x)) # add the candidate point to the matrix
      accept <- accept + 1 # add one to counter
    }
  }
  
  alpha <- accept/n # probability to be inside the k - dimensions sphere
  x1 <- accepts # accepted points to vector
  
  par(mfrow=c(1, 2)) # 2 graphs per window
  plot(density(accepts[ , 1]), main = "") # density plot of 1 dimension (accepted points)
  
  abline(v = mean(accepts[, 1]), col = "red", lwd = 1, lty = 3) # add the line of mean value
  text(mean(accepts[ , 1]), 0.1 , round(mean(accepts[ , 1]), 3)) # add the mean value as text
  
  plot(accepts, xlab = "x accepts", ylab = "y accepts") # plot accepted points (2 dimensions)
    
  return(list(alpha, accepts)) # the function returns probability and accepted points
  }


set.seed(17) # seed for regeneration
(mc.d(10000, 1, 10)) # 10000 samples in 10 dimensions


#### Calculate the mean of these values to verify where the probability converges ####
mc_meanx_d <- function(n, d, k, m) { # n = number of samples, U ~ [0 ,d], k = dimensions, m = iterations
  
  x  <- c()
  for (i in 1:m) {
    x[i] <- mc.d(n, d, k)[[1]]
    
    # Insert the data frame to a vector
  }
  z <- mean(x)
  
  #Finding the average of the acceptance rate 
  return(x) # probability of accepting points inside the hyper-sphere
}

set.seed(18) # seed for regeneration
system.time(z <- mc_meanx_d(10000, 1, 10, 100)) # running process for 10000 simulations, U ~ [-1, 1], k = 10, m = 100

# The mean value of the 100 independent samples from the Uniform distribution is the expected value
# of the distribution
# (f(x1) + ... + f(x100))/100 = E[f(x)] (LLN)

# I calculate the expected value of the 100 independent samples from the uniform distribution
# in 10 dimensions

dev.off() # close active windows (graphs)
mean(z)# The expected mean value of the distribution
hist(z, probability = TRUE, ylim = c(0, 1000), xlim = c(0, 0.005 ), main = "MC acceptance ratio, k - 10 dimensions")
# histogram of probabilities z (f(x1),..., f(x100))
lines(density(z),col="red") # The histogram of probability density is similar to normal distribution,

# We know that the estimation of the integrate is |B| = prob_estimated*(2^d) where d = 10 are the dimensions.
# So the |B| estimation is mean(z)*(2^10)

#### Estimate the integrate of B (10 dimensions) ####
B_estimated <- mean(z)*(2^10)
B_estimated

# We know from the calculation of the integrate with polar coordinates (10 dimensions) 
# that that the true area of sphere is 2.55
(B_true <- (pi^(10/2))/gamma(10/2 + 1))

# We estimated with Monte Carlo sampling of Uniform distribution 
# that integrate of B is 2.516679 close to the true volume of 10d sphere

#### The error of the integrate estimation is B_estimated - B_true ####
round((estimated_error <- abs(B_true - B_estimated)), 3)

### The relative error of estimation of area B  is 0.99%
(estimated_error_perc <- round(abs(1 - (B_estimated/B_true)), 4)*100)


#### Calculation of the number of samples needed for estimation accuracy in three decimals ####

sd(z) # standard deviation from 100 iterations of 10.000 samples

sd(z)/sqrt(10000) #estimation error

cumprob <- pnorm(z, mean(z), sd(z)) #cumulative probability of normal distribution
# with mean = mean(z) and sd = sd(z)

sd(cumprob) # standard deviation of cumulative probability for expected mean value of
# MC acceptance ratio

sd(cumprob)/sqrt(10000) #estimaton error < 4*10^-3

(round((sd(cumprob)/sqrt(10000))/(sd(z)/sqrt(10000)), 0)) # for precision > 2*10^-3 (or error < 2*10^-3)
# we have to increase six hundred times the samples of the iterating process, that means, 
# 6 million samples in each iterating process



#### Calculate the running time depending on the number of samples (dimensions k = 10) ####

mc.10d_run_time <-function(number_accepts, d, k)
{  
  
  accept <- 0 # counter for number of accepts
  data <- c(2*runif(k, 0, d) - 1) # vector of simulated point
  accepts <- matrix(data, 1, k)  # matrix includes the first point of the process
  
    n <- 0 # counter for number of steps
      
    while(accept < number_accepts + 1) { # while loop (accept < number of accepts)
        n <- n + 1 # incremental counter for number of steps
        x <- t(matrix(2*runif(k, 0, d) - 1)) # new simulated point
       if(sum(x^2) <= 1) { # if x is inside the sphere
      accepts <- rbind(accepts, x) # x is inserted into the matrix of accepts
      accept <- accept + 1 # and accept counter is increasing by 1
      
      }
      }
    alpha <- accept/n # probability of accepting points
  return(n) # return number of steps
  }


set.seed(22) # seed for regeneration
number_steps <- (mc.10d_run_time(1000, 1, 10)) # run process for 1000 accepted points, U ~ [-1, 1], and k = 10 dimensions


start.time <- Sys.time() # starting time
set.seed(31) # seed for regeneration

for (i in seq(2000, 10000, by = 1000)) {
number_steps <- cbind(number_steps, mc.10d_run_time(i, 1, 10)) # iterating process for i - accepted points
}
end.time <- Sys.time() # ending time
(time.taken <- end.time - start.time) # running time


slope <- data.frame(t(number_steps)) # assign results to a data frame (df)
head(slope) # first rows of df

c(seq(1, nrow(slope), by = 1)) # sequence 1 to 10, e.g. 1, 2, 3, .., 10

slope$accepts <- c(seq(1000, nrow(slope)*1000, by = 1000)) # add column with number of accepted points 1000, 2000, etc.

names(slope)[which(names(slope) %in% c("t.number_steps."))] <- c("number_steps") # set column name to new column


plot(slope$number_steps, slope$accepts, xlab = "Number of Steps", ylab = "Accepts", main = "Slope") # plot the (x,y) columns
# of df


y <- lm(slope$accepts ~ slope$number_steps, data = slope) # linear model between accepts and number of steps
summary(y) # summary of linear model

y$coefficients[2]# the coefficient of linear regression is the expected mean value of probability density function

slope$Probability <- slope$accepts/slope$number_steps # expected mean value of probability density function 
slope$Exp_Time <- round(1/slope$Probability, 0) # expected mean time for having a point inside the sphere 
# (return to the same state X = 1)
slope$Constant <- 1/slope$accepts # constant value
mean(slope$Constant)

qqnorm(slope$accepts) # normal probability q-q plot
qqline(slope$accepts) # added line in normal probability q-q plot


predict(y, newdata = slope) # predicted number of accepted points based on the number of steps (linear model)
slope$accepts # number of accepted points


#### The expected Time T ####

# The expected Time T Monte Carlo has to run (number of samples) in order to have X accepted points inside the circle

c <- 1/X # calculated in slope data frame
d <- 10 # dimensions

# For instance for X = 20,000 (accepts - points inside the 10d sphere)
X <- 20000

# The slope of the regresion model between accepts and number of steps is the
# expected mean value of the probability density function

(expected_time <- round(c*X/y$coefficients[2], 0))# the expected mean time for each accepting point 

# As a result the expected running time in order to have 20,000 accepting points is: 
expected_time*20000


new_row <- c(expected_time*X, X, 1/expected_time, expected_time, c) # new row of slope data frame (values for each column)

#### Adding the predicted number of steps in slope data frame ####
slope <- rbind.data.frame(slope, new_row)


#### Finding the variance of expected value in k - 10d ####
asym_var <- function(l, n, d, k, m){
  
  x <- array(dim = c(m, l))
  for (i in 1:l) {
    x[, i] <- (mc_meanx_d(n, d, k, m))
  }
 y <- x
   return(y)
}


set.seed(123) # seed for regeneration

start.time <- Sys.time() # starting time
y <- asym_var(10, 10000, 1, 10, 100) # set values of the asym_var function to a vector (y)
end.time <- Sys.time() # ending time
(time.taken <- end.time - start.time) # running time

#### Histogram of expected value of y (10 iterations of 100 samples with 10000 observations, k - dimensions = 10) ####
hist(y, breaks = 5, freq = FALSE, ylim = c(0, 800))
lines(density(y),col="red")# The histogram of probability density is similar to normal distribution,
apply(y, 2, var) # variance for each dimension

summary(y) # summary statistics of y


### Validation of histogram (probability density value)
h <- hist(y, breaks = 5, freq = FALSE)
sum(h$density) # 2000 distributed values to equivalent intervals

unique(zapsmall(diff(h$breaks))) # range for each interval
unique(zapsmall(diff(h$breaks)))*sum(h$density) # validation of probability density (multiplication has to be 1)

dev.off() # close open windows

### Normality function ####
normality <- function (x, k) { 
  y <- (x - (B_true/2^k))/sd(x)
return(y)
  }

### asymptotic variance in each iteration ###
apply(normality(y[ , ], 10), 2, var) 


### Histogram of y ####
hist(normality(y[ , ], 10) ,freq = FALSE, main = "Normalization of Expected Mean 
     Values (k - 10 dimensions)", xlab = "")

### Add curved line (density) ###
lines(density(normality(y[ , ], 10)),col="red")

### set seed for generating the same sample ###
set.seed(321)

### add line with sampling from Normal disribution ###
lines(density(rnorm(1000)), col = " green", lty = 3)

