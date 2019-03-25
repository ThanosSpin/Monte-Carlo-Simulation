
####Clear Environment####
rm(list = ls())

####Call Libraries####

library(doParallel)##Parallel Processing
library(psych)#describe function
library(knitr)##Create .md, .html, and .pdf files
library(markdown)##Create .md, .html, and .pdf files


####Parallel procesing for running performance####
####Set number of cores for parallel processing
no_cores <- detectCores() - 2 

####Register number of cores for parallel processing
registerDoParallel(no_cores)

####Validate number of cores for parallel processing
getDoParWorkers()

####Monte Carlo in 2d####

###Mc in 2d function###
mc.2d <-function(x, y, d, n)
  {  
  accept <- 0

  data <- c(x, y)

  accepts <- matrix(data, 1, 2)
    for (i in 2:n) {
      y <- 2*runif(1, 0, d) - 1
      x <- 2*runif(1, 0 ,d) - 1
      
      if((x^2 + y^2) <= 1) { 
    accepts <- rbind(accepts, c(x,y))
    accept <- accept + 1
  }
}
alpha <- accept/n
par(mfrow=c(1, 2))

# hist(accepts[ , 1], freq=FALSE, las = 1, main="", xlab="x", ylab="Probability density", ylim = c(0, 1))
plot(density(accepts[ , 1]), main = "")
index <- seq(-1, 1, 0.1)
lines(index, 0.5*index/index, col = "green")

# x2 <- seq(min(c(0, accepts[ , 1])), max(accepts[ , 1]), length.out=100)
# x2 <- x
# curve(dnorm(x, mean(accepts[ , 1]), sd(accepts[ , 1])), col = 2,lty = 2, add = TRUE)
plot(accepts, xlab="x accepts", ylab="y accepts")

return(list(alpha, accepts))
}

####Checking the system time####
dev.off()
set.seed(15)
system.time(mc.2d(0, 0, 1, 1000))
set.seed(16)
system.time(mc.2d(0, 0, 1, 10000))
set.seed(19)
system.time(mc.2d(0, 0, 1, 50000))
set.seed(20)
system.time(mc.2d(0, 0, 1, 100000))
set.seed(21)
system.time(mc.2d(0, 0, 1, 150000))



####calculate the mean of the generated samples (MC acceptance ratio) to verify where the chain converges####
mc_meanx_2d <- function(n, d, m) {
  
  y <- 2*runif(1, 0, d) - 1
  z <- 2*runif(1, 0, d) - 1
  x <- c()
  acc <- c()  
  for (i in 1:m) {
    x[i] <- mc.2d(y, z, d, n)[[1]] #Insert the data frame to a vector
    acc <-  mc.2d(y, z, d, n)[[2]][ , 1]*mc.2d(y, z, d, n)[[2]][ , 2] # joint Uniform distribution f(x)*f(y)
    }
  z <- mean(x) #finding the mean of the acceptance rate 
  return(list(x, acc)) # return the acceptance rates of the process
}


set.seed(14)
start.time <- Sys.time()
system.time(z <- mc_meanx_2d(10000, 1, 100))
end.time <- Sys.time()
(time.taken <- end.time - start.time)

###The mean value of the 100 independent samples from the distribution is the expected value
###of the distribution
###(f(x1) + ... + f(x100))/100 =.. E[f(x)]

####We calculate the expected value of the 100 independent samples from the uniform distribution####
####in 2 dimensions

dev.off() # close graphs per window
mean(z[[2]]) ###The expected value of the distribution


index <- seq(-1, 1, 0.1)

plot(density(z[[2]]), main ="") #plot probability density of z - vector f(x,y)
abline(v = mean(z[[2]]), col = "red", lwd = 1, lty = 3) # add the line of mean value
text(mean(z[[2]]), 0.1 , round(mean(z[[2]]), 3)) # add the mean value as text
#lines(index, 0.5*index/index, col = "green") # straight line at the mean value of the Unifrom distribution


plot(density(z[[1]]), main ="") #plot probability density of z - vector (MC acceptance ratio)
abline(v = mean(z[[1]]), col = "red", lwd = 1, lty = 3) # add the line of mean value
text(mean(z[[1]]), 0.1 , round(mean(z[[1]]), 3)) # add the mean value as text



mean(z[[1]])###The expected value of the MC acceptance ratio
dev.off()


hist(z[[1]], probability = TRUE, ylim = c(0, 120))###histogram of probabilities z (f(x1),..., f(x100))
lines(density(z[[1]]),col="red")####

sd(z[[2]]) # standard deviation of joint distribution of accepted values

####We know that the estimation of the integrate is |B| = prob_estimated*(2^d) where d = 2 are the dimensions.
####So the |B| estimation is mean(z)*(2^2)

####Estimate the integrate of B (2 dimensions)####
(B_estimated <- mean(z[[1]])*(2^2))


##We know from the calculation of the integrate with polar coordinates (2 dimensions) 
##that that the area of probability density function upon 2d sphere (B_true) is pi/4*(4)
(B_true <- pi)

##We estimated with Monte Carlo sampling of uniform distribution 
##that integrate of B is mean_estimated*(2^2) = mean(z)*4 = 3.141796 close to the true volume of 2d sphere (pi*R^2, pi for R = 1)
##Correspondingly, the area of probability density function upon 2d sphere (B_estimated) is mean(z)*2^2 
B_estimated == 2^2*mean(z[[1]])

####The error of the integrate estimation is B_estimated/B_true####
(estimated_error <- abs(B_true - B_estimated))
round(estimated_error, 4)

###The error of estimation of area B  is 0.0006  or 0.02%
(estimated_error_perc <- round(abs(1 - (B_estimated/B_true)), 4)*100)

###It is a precise algorithm (estimation error very low < 1%)

####Monte Carlo in k - dimensions####

###MC in k - dimensions###
mc.d <-function(n, d, k)
{  
  accept <- 0
  
  x <- (2*runif(k, 0, d) - d)
  data <- c(x)
  accepts <- matrix(data, 1, k)
  for (i in 1:n) {
    x <- t(matrix(2*runif(k, 0, d) - d)) 
    if(sum(x^2) <= 1) { 
      accepts <- rbind(accepts, c(x))
      accept <- accept + 1
    }
  }
  alpha <- accept/n
  x1 <- accepts
  
  par(mfrow=c(1, 2))
  plot(density(accepts[ , 1]), main = "")
  
  abline(v = mean(accepts[, 1]), col = "red", lwd = 1, lty = 3) # add the line of mean value
  text(mean(accepts[ , 1]), 0.1 , round(mean(accepts[ , 1]), 3)) # add the mean value as text
  
  plot(accepts, xlab="x accepts", ylab="y accepts")
    
  return(list(alpha, accepts))
  }


set.seed(17)
(mc.d(10000, 1, 10))


####calculate the mean of these values to verify where the probability converges####
mc_meanx_d <- function(n, d, k, m) {
  
  x  <- c()
  for (i in 1:m) {
    x[i] <- mc.d(n, d, k)[[1]]
    
    #Insert the data frame to a vector
  }
  z <- mean(x)
  
  #Finding the mean of the acceptance rate 
  return(x)
}

set.seed(18)
system.time(z <- mc_meanx_d(10000, 1, 10, 100))

###The mean value of the 100 independent samples from the distribution is the expected value
###of the distribution
###(f(x1) + ... + f(x100))/100 =.. E[f(x)]

####We calculate the expected value of the 100 independent samples from the uniform distribution####
####in 10 dimensions

dev.off()
mean(z)###The expected value of the distribution
hist(z, probability = TRUE, ylim = c(0, 1000), xlim = c(0, 0.005 ), main = "MC acceptance ratio, k - 10 dimensions")###histogram of probabilities z (f(x1),..., f(x100))
lines(density(z),col="red")####The histogram of probability density is similar to normal distribution,
# plot(density(z))
###the f(x) probabilities are normally distributed

# x1 <- seq(min(c(0,z)), max(z), length.out=100)
# lines(x2, dnorm(x2, mean(x1), sd(x1)), lty = 2, col =2)


####We know that the estimation of the integrate is |B| = prob_estimated*(2^d) where d = 10 are the dimensions.
####So the |B| estimation is mean(z)*(2^10)

####Estimate the integrate of B (10 dimensions)####
B_estimated <- mean(z)*(2^10)
B_estimated

##We know from the calculation of the integrate with polar coordinates (10 dimensions) 
##that that the area of B sphere is 2.55
(B_true <- (pi^(10/2))/gamma(10/2 + 1))

##We estimated with Monte Carlo sampling of uniform distribution 
##that integrate of B is 2.516679 close to the true volume of 10d sphere

####The error of the integrate estimation is B_estimated/B_true####
round((estimated_error <- abs(B_true - B_estimated)), 3)

###The error of estimation of area B  is 0.025 or 0.99%
(estimated_error_perc <- round(abs(1 - (B_estimated/B_true)), 4)*100)


####Calculation of the number of samples####

sd(z) # standard deviation from 100 iterations of 10.000 samples

sd(z)/sqrt(10000) #estimation error

cumprob <- pnorm(z, mean(z), sd(z)) #cumulative probability of normal distribution
# with mean = mean(z) and sd = sd(z)

sd(cumprob) # standard deviation of cumulative probability for expected mean value of
# MC acceptance ratio

sd(cumprob)/sqrt(10000) #estimaton error < 4*10^-3

(round((sd(cumprob)/sqrt(10000))/(sd(z)/sqrt(10000)), 0)) # for precision > 2*10^-3 (or error < 2*10^-3)
# we have to increase six hundred times the samples of the iterating process, that means, 
#6 million samples in each iterating process



####Calculate the running time depending on the number of samples (dimensions k = 10)####

mc.10d_run_time <-function(number_accepts, d, k)
{  
  
  accept <- 0 #counter for number of accepts
  data <- c(2*runif(k, 0, d) - 1) #vector of simulated point
  accepts <- matrix(data, 1, k)  #matrix includes the first point of the process
  
    n <- 0 # counter for number of steps
      
    while(accept < number_accepts + 1) { #while loop (accept < number of accepts)
        n <- n + 1 #incremental counter for number of steps
        x <- t(matrix(2*runif(k, 0, d) - 1)) # new simulated point
       if(sum(x^2) <= 1) { #if x is inside the sphere
      accepts <- rbind(accepts, x) #x is inserted into the matrix of accepts
      accept <- accept + 1 #and accept counter is increasing by 1
      
      }
      }
    alpha <- accept/n # probability of accepting points
  return(n) #return number of steps
  }


set.seed(22)
number_steps <- (mc.10d_run_time(1000, 1, 10))


start.time <- Sys.time()
set.seed(31)
for (i in seq(2000, 10000, by = 1000)) {
number_steps <- cbind(number_steps, mc.10d_run_time(i, 1, 10))
}
end.time <- Sys.time()
(time.taken <- end.time - start.time)


slope <- data.frame(t(number_steps))
head(slope)

c(seq(1, nrow(slope), by = 1))

slope$accepts <- c(seq(1000, nrow(slope)*1000, by = 1000))

names(slope)[which(names(slope) %in% c("t.number_steps."))] <- c("number_steps")

#plot(slope$number_steps, slope$accepts)
plot(slope$number_steps, slope$accepts, xlab = "Number of Steps", ylab = "Accepts", main = "Slope")


y <- lm(slope$accepts ~ slope$number_steps, data = slope)
summary(y)

y$coefficients[2]###the coeff of linear regression is the expected mean value of probability density function

slope$Probability <- slope$accepts/slope$number_steps # expected mean value of probability density function 
slope$Exp_Time <- round(1/slope$Probability, 0) # expected mean time for having a point inside the sphere (return to the same state X = 1)
slope$Constant <- 1/slope$accepts
mean(slope$Constant)

qqnorm(slope$accepts)###normal probability q-q plot
qqline(slope$accepts)###added line in normal probability q-q plot


predict(y, newdata = slope)
slope$accepts


####The expected Time T####

###The expected Time T Monte Carlo has to run (number of samples) in order to have X accepts

c <- 1/X ###calculated in slope data frame
d <- 10 ##dimensions

###For instance for X = 20,000 (accepts - points inside the 10d sphere)###
X <- 20000

###The slope of the regresion line between accepts and number of steps is the
###Expected value of the probability density function
(expected_time <- round(c*X/y$coefficients[2], 0))###the expected mean time for each accepting point 

###As a result the expected running time in order to have 20,000 accepting points is: 
expected_time*20000


new_row <- c(expected_time*X, X, 1/expected_time, expected_time, c) # new row of slope data frame (values for each column)

#### Adding the predicted number of steps in slope data frame ####
slope <- rbind.data.frame(slope, new_row)


####Finding the variance of expected value in k - 10d####
asym_var <- function(l, n, d, k, m){
  
  x <- array(dim = c(m, l))
  for (i in 1:l) {
    x[ ,i] <- (mc_meanx_d(n, d, k, m))
  }
 y <- x
   return(y)
}


set.seed(123)
start.time <- Sys.time()
y <- asym_var(10, 10000, 1, 10, 100)
end.time <- Sys.time()
(time.taken <- end.time - start.time)

####Histogram of expected value of y (10 iterations of 100 samples with 10000 observations, k - dimensions = 10)####
hist(y, breaks = 5, freq = FALSE, ylim = c(0, 800))
lines(density(y),col="red")####The histogram of probability density is similar to normal distribution,
apply(y, 2, var)

summary(y)###summary statistics of y


###Validation of histogram (probability density value)
h <- hist(y, breaks = 5, freq = FALSE)
sum(h$density)

unique(zapsmall(diff(h$breaks)))
unique(zapsmall(diff(h$breaks)))*sum(h$density)
# rm(h)
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


####Creating the script as a report in order to be able to read it easily####
# Create .md, .html, and .pdf files
knit("MCsampling_v2.Rmd")
markdownToHTML('MCsampling_v2.md', 'MCsampling_v2.html', options=c("use_xhml"))
