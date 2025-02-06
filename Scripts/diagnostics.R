## Compute MCMC diagnostic information
# Copyright Brian Keller 2025, all rights reserved
# License: GPL-3

## Calculate effective sample size 
#' @param x a matrix with rows equal to the iterations and columns equal to the
#'          chains
n_eff <- function(x) {
    
    # Declare variables
    n <- NROW(x)
    m <- NCOL(x)
    
    # Compute Chain Means
    t_j <- colMeans(x)
    
    # Compute Autocovariance
    #' @param x a one dimensional vector
    chain_cov <- apply(x, 2,  function(x) {
        # Get size
        n <- NROW(x)
        
        # Add Padding
        x <- c(x - mean(x), rep(0, n))
        
        # Use FFT algorithm to compute autocors
        fx <- fft(x)
        
        # Need to normalize inverse because R doesn't do it
        tmp <- fft(Conj(fx) * fx, inverse = TRUE)  / length(fx)
        
        # Subset and Convert to Real
        x <- Re(tmp[seq_len(n)]) * 2 * n
        
        # Use Geyer 1992 estimate
        # see https://projecteuclid.org/euclid.ss/1177011137
        # Return covariance (NOTE: For cor would need to do x / x(0))
        x / (n * n * 2.0)
    })
    
    # Compute var+
    var_plus <- mean(chain_cov[1, ])
    
    # Compute average variance
    var_bar <- var_plus * (n / (n - 1))
    
    # Add between chain if needed
    if (m > 1) var_plus <- var_plus + var(t_j)
    
    # Set up
    r_t <- rep(0, n)
    r_odd <- r_even <- 1.0
    i <- 0
    
    # Compute rho hat
    r_t[1L] <- r_odd
    r_even   <- 1 - (var_bar - mean(chain_cov[2L, ])) / var_plus
    r_t[2L] <- r_even
    
    # Loop over
    for (i in seq(3, n - 4, by = 2)) {
        # Compute even and odd
        r_odd  <- 1 - (var_bar - mean(chain_cov[i    , ])) / var_plus
        r_even <- 1 - (var_bar - mean(chain_cov[i + 1, ])) / var_plus

        # Handle NA values        
        if (is.na(r_odd + r_even)) return(NaN)
        # Exit if negative
        if ((r_odd + r_even) < 0)  break
        
        # Update
        r_t[i]     <- r_odd
        r_t[i + 1] <- r_even
    }
    i_max <- i + 2
    if (r_odd > 0) r_t[i_max] <- r_odd
    
    if (i_max - 3 > 3) {
        for (i in seq(3, i_max - 3, by = 2)) {
            if ((r_t[i] + r_t[i + 1]) > (r_t[i - 1] + r_t[i - 2])) {
                r_t[i]     <- (r_t[i - 1] + r_t[i - 2]) / 2
                r_t[i + 1] <- r_t[i]
            }
        }
    }

    # Compute tau
    tau_h <- -1 + 2 * sum(r_t[seq_len(i_max - 1)]) + r_t[i_max]
    tau <- max(tau_h,  1 / log10(n * m)) 
    
    # Return effective sample size
    n * m / tau
}


# Compute psrf
#' @param x a matrix with rows equal to the iterations and columns equal to the
#'          chains
#' @param split_chain a logical (default TRUE) indicating if we should take
#'                    the split chain PSRF (i.e., taking the first half and
#'                    comparing to the second half).
psrf <- function(x, split_chain = TRUE) {
    # Coerce x to matrix
    x <- as.matrix(x)

    # Set up params
    n <- ifelse(split_chain, floor(NROW(x) / 2), NROW(x))
    m <- ifelse(split_chain, NCOL(x) * 2, NCOL(x))
    # If not enough draws error
    if (n < 2) stop("Not enough draws")

    # Create matrix of entire parameters
    chaindat <- do.call("cbind", apply(x, 2, function(x) {
        if (split_chain) {
            return(cbind(x[1:n], x[(n + 1):(2 * n)]))
        } else {
            return(x[1:n])
        }
    }, simplify = F))

    # Compute PSR
    t_j <- colMeans(chaindat)
    t_bar <- mean(t_j)

    # Compute between
    B <- (1 / (m - 1)) * sum((t_j - t_bar)^2)

    # Compute within
    W <- sum((sweep(chaindat, 2, t_j)^2) / (n - 1))
    W <- ((n - 1) / n) * (W / m)

    # Compute psr
    return(sqrt((W + B) / W))
}

