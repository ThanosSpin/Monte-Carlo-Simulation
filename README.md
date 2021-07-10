
# Monte Carlo (MC) Simulations for Volume Estimation of a d-dimensional Sphere

The goal of this repository is the use of MC Simulations for volume estimation of a d-dimensional Sphere.

To begin with, I use MC in order to estimate the integral of a circle in two dimensions with radius R = 1. 
Then, I repeat the process for more dimensions (d = 1,.., d = 200). 
This volume can be estimated by using independent Bernoulli trials for approximating the probability of a simulated vector of U[-1,1]^d (d-dimensions) to be inside the sphere (Sd). The test is if the simulated values x1^2 + x2^2 + ... + xd^2 are <= 1 or not, d - dimensions, for R = 1.
Knowing that P = |Sd|/|Cd| the volume of |Sd| can be estimated: |Sd| = P*|Cd| (|Sd|: Volume of Sphere, |Cd|: Volume of |Cd|) 

In addition, I implement a corresponding process for MCMC (Markov Chain Monte Carlo Simulation). 

MCMC process will be uploaded in a separate repository.
