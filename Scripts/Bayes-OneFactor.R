#-------------------------------------------------------------------#
# Automatically Derived Full Conditionals for Factored Regression 
#   Models with Missing Data and Latent Variables
#
# Psychoco 2025
#
# Bayes-OneFactor.R
#
# Copyright Brian T. Keller 2025 all rights reserved
# License: GPL-3
#-------------------------------------------------------------------#


### Generate Multivariate Data
mrnorm <- function(n, mu, cov.mat) {
    p <- length(mu)
    if (p != nrow(cov.mat)) 
        stop("covariance matrix not compatible with mean vector")
    mu.mat <- matrix(1, nrow = n) %*% mu
    X <- mu.mat + matrix(rnorm(p*n), ncol = p) %*% chol(cov.mat)
    return(X)
}

# Build  n_effective and psrf 
source('diagnostics.R')

# Posterior mode
mode_v <- function(v){
    dens.y <- density(v)
    mode <- dens.y$x[order(dens.y$y,decreasing=T)][1]
    return(mode)
}

#  HPD: Monte Carlo method
emp_hpd <- function (theta) {
    alpha <- 0.95
    alpha <- min(alpha, 1 - alpha) 
    n <- length(theta) 
    L.start <- round(n * alpha) 
    theta <- sort(theta) 
    e <- theta[(n - L.start + 1):n] - theta[1:L.start] 
    ind <- which(e == min(e))[1] 
    return(c(`hpd 2.5%` = theta[ind], `hpd 97.5%` = theta[n - L.start + ind])) 
}

#### Gibbs Sampler ####

## bayes Classes
create_linear <- function(Y, X = c(), b = rep(0, length(X) + 1), s_e = 1) {
    structure(list(b = b, s_e = s_e),
              class = 'blinear', X = X, Y = Y)
}

create_latent <- function(
        Y, X = c(),
        intercept = FALSE,
        fix_var   = FALSE,
        b = rep(0, length(X) + ifelse(intercept, 1, 0)),
        s_e = 1
) {
    structure(
        list(b = b, s_e = s_e),
        class = 'blatent',
        X = X, Y = Y,
        intercept = intercept,
        fix_var = fix_var
    )
}

## Results class
as.results <- function(x,...) {
    structure(x, class = 'results')
}

summary.results <- function(object, ...) {
    object <- t(apply(object, 1, function(x){
        c(mean = mean(x),
          median = median(x),
          mode = mode_v(x),
          sd = sd(x),
          quantile(x, probs = c(0.025, 0.975)),
          emp_hpd(x),
          n_eff = n_eff(as.matrix(x)),
          psrf  = psrf(as.matrix(x))
        )
    }))
    return(object)
}

print.results <- function(x, ...) {
    print(summary(x), ...)
}

# Function for specific data set up
plot.results <- function(results, param = NULL, len = seq_len(NCOL(results)), ...) {
    if (is.null(param)) {
        for(p in rownames(results)) {
            plot(results, param = p, len, ...)
            invisible(readline(prompt="  Press [enter] to continue"))
        }
    } else {
        plot(results[param,len], type = 'l',
             ylab = 'Parameter', main = param,
             xlab = 'Iterations', ...)
    }
}

## predict method
predict.blinear <- function(object, data) {
    X <- cbind(1, as.matrix(data[,attr(object,'X'), drop = F]))
    return(X %*% object$b)
}

predict.blatent <- function(object, data) {
    X <- as.matrix(data[,attr(object,'X'), drop = F])
    if (attr(object,'intercept')) { X <- cbind(1, X) }
    if (NCOL(X) == 0) {
        return(rep(0, NROW(data)))
    }
    return(X %*% object$b)
}

## Set up generic
f <- function(model, ...) {
    UseMethod("f")
}

# Compute prediction
f.blinear <- function(model,...) {
    predict(model, data.frame(...))
}

## Set up generic
reg <- function(object, data) {
    UseMethod("reg")
}

### Set up one iteration
reg.blinear <- function(p, data) {
    
    # Get matrices
    X <- cbind(1, as.matrix(data[,attr(p,'X'), drop = F]))
    Y <- as.matrix(data[,attr(p,'Y'), drop = F])
    
    ## Draw Betas
    inv_XtX <- solve(crossprod(X))
    bhat <- inv_XtX %*% t(X) %*% Y
    S_b  <- p$s_e^2 * inv_XtX
    p$b <- c(mrnorm(1, c(bhat), S_b)) # Draw
    
    ## Draw Residual Variance
    eTe <- crossprod(Y - X %*% p$b)
    s2_e <- 1 / rgamma(1, nrow(Y) / 2, eTe / 2)
    p$s_e <- sqrt(s2_e)
    return(p)
}

reg.blatent <- function(p, data) {
    
    # Get matrices
    X <- as.matrix(data[,attr(p,'X'), drop = F])
    
    if (attr(p,'intercept')) { X <- cbind(1, X) }
    Y <- as.matrix(data[,attr(p,'Y'), drop = F])
    
    ## Draw Betas
    if (NCOL(X) != 0) {
        inv_XtX <- solve(crossprod(X))
        bhat <- inv_XtX %*% t(X) %*% Y
        S_b  <- p$s_e^2 * inv_XtX
        p$b <- c(mrnorm(1, c(bhat), S_b)) # Draw
    }
    ## Draw Residual Variance
    if (!attr(p,'fix_var')) {
        eTe <- crossprod(Y - X %*% p$b)
        s2_e <- 1 / rgamma(1, nrow(Y) / 2, eTe / 2)
        p$s_e <- sqrt(s2_e)
    }
    return(p)
}

# Solve distributions
dist <- compute_dist(
    y = 'mydata$y1',
    a = 'y1$b[1]',
    b = 'y1$b[2]',
    var_y  = 'y1$s_e^2',
    mu_x = 'predict(fy, mydata)',
    var_x  = 'fy$s_e^2'
)
dist <- compute_dist(
    y = 'mydata$y2',
    a = 'y2$b[1]',
    b = 'y2$b[2]',
    var_y  = 'y2$s_e^2',
    mu_x = dist$mean,
    var_x  = dist$variance
)
dist <- compute_dist(
    y = 'mydata$y3',
    a = 'y3$b[1]',
    b = 'y3$b[2]',
    var_y  = 'y3$s_e^2',
    mu_x = dist$mean,
    var_x  = dist$variance
)


### Function to Run Gibbs Sampler
gibbs <- function(starts, mydata, reps = 1000, seed ) {
    # Set Seed
    if(!missing(seed)) { set.seed(seed) }
    
    # Set up starts
    p <- lapply(seq_len(reps + 1), function(x) {
        list(
            y1 = NULL,
            y2 = NULL,
            y3 = NULL,
            y4 = NULL,
            fy = NULL
        )
    })
    p[[1]] <- starts
    
    # Run Sampler
    for (i in seq_len(reps) + 1) {
        
        # Set up environment
        env <- list2env(p[[i-1]])
        env$mydata <- mydata
        
        # Draw latent
        mydata$fy <- rnorm(
            NROW(mydata),
            dist$mean |> str2lang() |> eval(env),
            dist$variance |> str2lang() |> eval(env) |> sqrt()
        )

        # Draw Parameters
        p[[i]]$y1 <- reg(p[[i-1]]$y1, mydata)
        p[[i]]$y2 <- reg(p[[i-1]]$y2, mydata)
        p[[i]]$y3 <- reg(p[[i-1]]$y3, mydata)
        p[[i]]$fy <- reg(p[[i-1]]$fy, mydata)
    }
    
    # Obtain Values as Matrix (exclude starts [-1])
    results <- sapply(p[-seq_len(round(reps/2) + 1)], function(x) {
        with(x, {
            c(
                y1.b   = y1$b,
                y1.s_e = y1$s_e,
                y2.b   = y2$b,
                y2.s_e = y2$s_e,
                y3.b   = y3$b,
                y3.s_e = y3$s_e,
                y4.b   = y4$b,
                y4.s_e = y4$s_e,
                fy.b   = fy$b,
                fy.s_e = fy$s_e
            )
        })
    })
    
    # Return results
    return(as.results(results))
}


