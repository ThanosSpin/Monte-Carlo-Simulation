
# Monte Carlo (MC) Simulations for Volume Estimation

The goal of this repository is the use of MC Simulations for volume estimation. 
I focus on the shape of Sphere.

To begin with, I use MC in order to estimate the integral of a circle in two dimensions with radius R = 1. 
Then, I repeat the process for more dimensions (k = 1,.., k = 200). 

First, I calculate the integral in closed-form (true volume) and then, I simulate values from Uniform distribution, 
U ~ [-1,1], based on the circle model x1^2 + x2^2 + ... + xk^2 <= 1, k - dimensions, for R = 1. 

In addition, I implement a corresponding process for MCMC (Markov Chain Monte Carlo Simulation). 

MCMC process will be uploaded in a separate repository.
