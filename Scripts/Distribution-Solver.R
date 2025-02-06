#-------------------------------------------------------------------#
# Automatically Derived Full Conditionals for Factored Regression 
#   Models with Missing Data and Latent Variables
#
# Psychoco 2025
#
# Distribution-Solver.R
#
# Copyright Brian T. Keller 2025 all rights reserved
# License: GPL-3
#-------------------------------------------------------------------#


## Compute Symbolic Distribution for x | y if
#   y ~ normal(a + b * x, var_y);
#   x ~ normal(mu_x, var_x);
#
compute_dist <- function(y, a, b, var_y, mu_x, var_x) {
    # Obtain substitutions
    y <- substitute(y)
    a <- substitute(a)
    b <- substitute(b)
    mu_x <- substitute(mu_x)
    var_y <- substitute(var_y)
    var_x <- substitute(var_x)
    
    # Check if they exist and obtain values, otherwise use name
    y <- tryCatch(eval(y), error = \(.) paste('(',deparse(y),')'))
    a <- tryCatch(eval(a), error = \(.) paste('(',deparse(a),')'))
    b <- tryCatch(eval(b), error = \(.) paste('(',deparse(b),')'))
    mu_x <- tryCatch(eval(mu_x), error = \(.) paste('(',deparse(mu_x),')'))
    var_y <- tryCatch(eval(var_y), error = \(.) paste('(',deparse(var_y),')'))
    var_x <- tryCatch(eval(var_x), error = \(.) paste('(',deparse(var_x),')'))
    
    # Solve
    list(
        mean = paste(
            '(', var_x, '*', b, '*', '(', y, '-', a, ')', '+',
            var_y, '*', mu_x, ')', '/',
            '(', var_x, '*', b, '*', b, '+', var_y, ')'
        ),
        variance = paste(
            '(', var_y, '*', var_x, ')', '/', '(', 
            var_x, '*', b, '*', b, '+', var_y, ')'
        )
    )
}
