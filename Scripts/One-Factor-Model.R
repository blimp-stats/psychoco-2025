#-------------------------------------------------------------------#
# Automatically Derived Full Conditionals for Factored Regression 
#   Models with Missing Data and Latent Variables
#
# Psychoco 2025
#
# One-Factor-Model.R
#
# Copyright Brian T. Keller 2025 all rights reserved
# License: GPL-3
#-------------------------------------------------------------------#

### Generate Multivariate Normal
mrnorm <- function(n, mu, cov.mat) {
    p <- length(mu)
    if (p != nrow(cov.mat)) stop("covariance matrix not compatible with mean vector")
    mu.mat <- matrix(1, nrow = n) %*% mu
    X <- mu.mat + matrix(rnorm(p*n), ncol = p) %*% chol(cov.mat)
    return(X)
}
### Generates Covariance Matrix
covariance_mat <- function(sds = c(1, 1), rho = 0.3) {
    cov.mat <- matrix(rho, length(sds), length(sds)) # Create num_vars by num_vars matrix
    diag(cov.mat) <- 1 # Set diagonal to 1
    sd.mat <- diag(length(sds))
    diag(sd.mat) <- sds
    return( sd.mat %*% cov.mat %*% sd.mat) # Compute Covariance Matrix
}

# Set up
set.seed(1597231)
N <- 1000

# Generate one latent variable
latent <- rnorm(N, mean = 0, sd = 1)

# Set up y indicators
y <- latent + mrnorm(
    N,
    mu = 5:7,
    cov.mat = covariance_mat(rep(1,3), rho = 0.0)
)

# Create data
mydata <- data.frame(y = y)
names(mydata) <- paste0('y', 1:3)
# write.csv(mydata, 'mydata.csv', row.names = F)


# Estimate in R
source('Distribution-Solver.R')
source('Bayes-OneFactor.R')

# Create Starts
starts  <- list(
    y1 = create_linear('y1', c('fy'), b = c(0,1)),
    y2 = create_linear('y2', c('fy')),
    y3 = create_linear('y3', c('fy')),
    fy = create_latent('fy', intercept = F, fix_var = T)
)

# Set up latent starts
mydata$fy <- rnorm(N)

# Run Gibbs
results <- gibbs(starts, mydata, reps = 10000, seed = 197231)

# Look at results and plots
print(results, digits = 3)
# plot(results)

# Check with lavaan
lav_mdl <- lavaan::sem(
    'f_y =~ NA*y1 + y2 + y3; f_y ~~ 1*f_y;',
    mydata,
    meanstructure = TRUE
)
lavaan::summary(lav_mdl)
