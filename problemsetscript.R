################## SCRIPT: Problem set of non-parametric statistics
################## Students: Arturo Prieto Tirado and Sofia Sorbet Santiago



############# Loading required libraries
pacman::p_load(rotasym, ks, matlib, rootSolve)


###########################################################################
##########################   EXERCISE 1 ###################################
###########################################################################


############## a. 

# Load the dataset
data(sunspots_births, package = "rotasym") 

# Covariate of interest: phi
sunspots_phi = sunspots_births$phi

# The default h is the one obtained using DPI selector which can also be obtained by ks::hpi
kde_phi = ks::kde(x = sunspots_phi)

# Plotting the kde of phi
plot(kde_phi, lwd = 3, xlab = "Phi")



############## b. 

# The kdde determines the density derivative of order deriv.order (0 is also the density)
# using the DPI bandwidth adequate for that derivative order. 
# (remember that higher order derivatives need bigger bandwidths)
kde_derivative_phi = ks::kdde(x = sunspots_phi, deriv.order = 1)
est = kde_derivative_phi$estimate 
evpoints = kde_derivative_phi$eval.points 

# Plot the kdde
plot(kde_derivative_phi, lwd = 3, xlab = "Phi")

# Function that calculates the kdde for phi at a given point
# @params x vector of evaluation points
# @return The value of the kdde at x
solvekde = function(x){
  kde_derivative_phi = ks::kdde(x = sunspots_phi, deriv.order = 1, eval.points = x)
  return(kde_derivative_phi$estimate)
}

# Find the relative extremes using uniroot.all in interval (-0.5, 0.5)
relativeextremes = uniroot.all(solvekde, lower = -0.5, upper = 0.5)



############## c. 

# Covariate of interest: total_area 
area = sunspots_births$total_area

# Normal Scaled bandwidth selector
bw = bw.nrd(x = area)

# Kde for area using previous bw
kde_area = ks::kde(x = area, h = bw, adj.positive = 1)

# Plotting the kde
plot(kde_area, col = "red")



############## d.

# Simulating 10.000 samples from the kde of phi covariate 
samp_kde = ks::rkde(n = 1e4, fhat = kde_phi)

# Histogram of those samples 
hist(samp_kde, main="Histogram of phi kde sampling", xlab = "Phi", xlim = c(-1,1))





###########################################################################
##########################   EXERCISE 2 ###################################
###########################################################################


############## a.

# Function that determines the cubic estimator for a sample data and returns its
# value at given points
# @params x vector of evaluation points
# @params data sample data with response and predictors as columns
# @params h bandwidth
# @return The cubic estimator from data with bandwidth h evaluated at x
cubic_estimator = function(x, data, h){
  
  # Initialization of some parameters using the arguments of the function:
  Y = data[,2] # response variable
  n=length(Y)  # number of observations 
  p=3          # order of the polynomial estimator
  
  # Procedure: finding the estimated betas using the least squares
  # using the solution (3.12) and normal kernel for each evaluation point
  
  # Initialization of the estimation of regression function  
  m = rep(0, times = length(x)) # as many points as evaluation points
  
  # Estimation of regression function at each point 
  for(i in 1:length(x)){
    
    # Define the X matrix for the evaluation point x[i]
    X = matrix(nrow=n, ncol = p+1)  # initialization of the X matrix, size (nx(p+1))
    X[,1]=1                         # first column filled with 1s
    for(j in (2:(p+1))){
      X[,j]=(data[,1]-x[i])^(j-1)   # computing the matrix X
    }
    
    # Initialize vector of beta estimates
    beta = rep(0, times = (p+1))
    
    # Use of both scaled kernel and standard normal
    W = 1/h * diag(dnorm((data[,1]-x[i])/h))
    
    # Least squares solution
    beta = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Y
    
    # Recall that beta=(m,m',m'',m''') so:
    m[i] = beta[1]
  }
  # Return the estimation of the cubic estimator for each x
  return(m)
}



############## b.

# Testing the implementation: 
# estimating the regression function for given X and Y

# Required initializations 
set.seed(1) # seed for simulation                        
n=500       # sample size

# Simulations 
X_sim = rnorm(n, mean = 1, sd = 1) # X follows a N(1,1)
eps = rnorm(n, mean = 0, sd = 0.5) # epsilon follows a N(0, 0.5)
Y_sim = (X_sim - 1)^2 + eps        # Y = m(x) + epsilon 

# Parameters for the cubic_estimator function
data = cbind(X_sim, Y_sim) # matrix of size nx2 with simulated (X, Y)
evaluation_points = seq(min(X_sim), max(X_sim), by = 0.001) # evaluation points

# Use of cubic_estimator function  
Yhat = cubic_estimator(x = evaluation_points, data = data, h = 0.5)

# Plotting last results 
plot(X_sim, Y_sim, main = "Cubic Estimator", xlab = "x", ylab = "y")
lines(x = evaluation_points, y = Yhat, col = "red")




###########################################################################
##########################   EXERCISE 3 ###################################
###########################################################################



############## a.

# Loading the required dataset
temps7 = read.csv(url("https://raw.githubusercontent.com/egarpor/handy/master/datasets/temps-7.txt"))

# Plotting the histogram of working temperatures of problematic smartphones
hist(temps7$x, main = "Histogram of working temperatures of problematic smartphones", xlab = "Temperature")

# Function that computes the scaled gradient at a given point
# @params x evaluation point
# @params h bandwidth
# @params data vector of data points
# @returns the value of the scaled gradient at x
eta_hat = function(x, h, data){
  
  # Idea: estimate kde and obtain its gradient using numDeriv::grad
  # as numDeriv::grad takes as arguments another function, declare this function
  
  # Function that estimates the kde of the data at point x
  # @params x evaluation point
  # @returns the value of the kde at x
  fhatsolve=function(x){
    return(ks::kde(x = data, h = h, eval.points = x)$estimate)
  }
  
  # Use fhatsolve to obtain the value of the gradient in each x point  
  gradient = numDeriv::grad(func = fhatsolve, x = x)
  
  # kde of data using as bw the h argument
  fhat = fhatsolve(x)
  
  # Returning the value of the scaled gradient
  return(h^2*gradient/fhat)
}



############## b.

# Function that estimates the scaled gradient at given point using the interpolated spline 
# @params x evaluation points
# @params interpolated_spline fitted spline function
# @returns the value of the scaled gradient at x using interpolated spline
eta_hat_spline=function(x, interpolated_spline){
  # Use of the function interpolated_spline
  # interpolated_spline must be declare previously  using splinefun 
  return(interpolated_spline(x = x))
}


### Comparison between eta_hat and eta_hat_spline functions

## eta_hat_spline estimation

# Initialization of the evaluation points for the fit of the spline
xvalues = seq(0,60, by = 1)

# Evaluation of those points using eta_hat function
yvalues = eta_hat(xvalues, h = 0.5, data = temps7$x)

# Fit the spline function using the pairs (xvalues, yvalues)
interpolated_spline = splinefun(x = xvalues, y = yvalues)

# Initialization of the points in which evaluate both the interpolation
# and the original eta_hat function
est_grid = seq(0,60, by = 0.1)

# eta_hat_spline performance: evaluation of the scaled gradient using interpolation
# this function evaluates the gradient for the set of points est_grid
# internally uses interpolated_spline 
spline_est = eta_hat_spline(x = est_grid, interpolated_spline)

## eta_hat estimation
# eta_hat performance: evaluation of the scaled gradient 
# the same grid of evaluation points as for eta_hat_spline will be used for comparison purposes
eta_est = eta_hat(x = est_grid, h = 0.5, data = temps7$x)

## visual comparison
# Comparing the estimated scaled gradients superimposing the 
par(mfrow = c(1,1))
plot(x = est_grid, y = eta_est)
lines(x = est_grid, y = spline_est, col = "red")
legend("bottomright", legend = c("eta_hat", "eta_hat_spline"), 
       col = c("black", "red"), lty = c(1,1), cex = 0.7)



############## c.

# Function that applies the euler method to calculate the trajectory of a given initial point
# @params x initial point
# @params h bandwidth
# @params data vector of data points
# @params N number of iterations
# @params epsilon absolute tolerance
# @params xgrid sequence of points for evaluating function eta_hat
# @return the final value to which the inital point x converges (close to a mode)
euler = function(x, data, h, N, epsilon, xgrid){
  
  # gradient estimated values using eta_hat function
  grad_est = eta_hat(xgrid, h = h, data = data)
  
  # function for interpolation using as grid (xgrid, grad_est)
  interpolated_spline = splinefun(x = xgrid, y = grad_est)
  
  # initialization of parameters
  phi = numeric(length = (N + 1))
  phi[1] = x
  
  # Euler method
  for (t in 1:N) {
    # Estimated gradient interpolating the point
    Df = eta_hat_spline(x = phi[t], interpolated_spline)
    # t iteration of Euler algorithm
    phi[t + 1] = phi[t] + h * Df 
    # Condition for breaking the loop
    if( abs(phi[t + 1] - phi[t]) < epsilon) return(phi[t + 1])
  }
  return(NA)
}



############## d.

# Suitable bandwidth
# @params p number of dimensions
# @params r derivative order
# @params sigma covariance matrix
# @params n number of data points
# @return the normal scaled bandwidth matrix (or scalar in univariate case)

suitable_bw = function(p, r, sigma, n) {
  return((4/(p + 2 * r + 2))^(2 / (p + 2 * r + 4)) * n^(- 2 / (p + 2 * r + 4)) * sigma)
}

# x sequence for which obtaining their mode:
x = temps7$x

# Initialization of the modes for each point
mode = numeric(length = length(x))

# Normal scaled bandwidth: Computing a suitable bw for derivative estimation
p = 1 # number of dimensions 
r = 1 # derivative order
n = length(temps7$x) # sample size
sigma = var(temps7$x) # variance of the data
H = suitable_bw(p = p, r = r, n = n, sigma = sigma)  # bw for gradient

# Loop for computing the mode for each point from x
for(i in 1:length(x)){
  mode[i] = euler(x[i], data = temps7$x, h = H, N = 500, epsilon = 1e-4,
                  xgrid = seq(6.8,63, by = 1))
}

# Plotting the final results using histograms and the sample 
hist(mode, breaks = 20)
rug(mode)



############## e.

# use k means with k=2 clusters on the convergence points of euler
k = 2
kms = kmeans(x = na.omit(mode), centers = k)

# see the center of the cluster (average)
kms$centers

# associate defect cluster with high temperature
defect_cluster = which.max(kms$centers)

# calculate the proportion of defective phones
prop = sum(kms$cluster == defect_cluster)/length(kms$cluster)
prop