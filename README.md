
# Monte Carlo (MC) Simulations for Volume Estimation

The purpose of this repository is the use of MC Simulations for volume estimation. We focus on the shape of Sphere.

To begin with, we use MC in order to calculate the integral of a circle in two dimensions with radius R = 1. 
Then, we repeat the process for more dimensions. 

First, we calculate the integral in closed-form (on paper) and then, we simulate values from Uniform distribution, 
U ~ [0 , 1], based on the circle model x1^2 + x2^2 + ... + xn^2 <= 1, for R = 1. 

In addition, we implement a corresponding process for MCMC (MCMC repository).
