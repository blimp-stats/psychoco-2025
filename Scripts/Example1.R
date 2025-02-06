#-------------------------------------------------------------------#
# Automatically Derived Full Conditionals for Factored Regression 
#   Models with Missing Data and Latent Variables
#
# Psychoco 2025
#
# Example1.R
#
# Copyright Brian T. Keller 2025 all rights reserved
# License: GPL-3
#-------------------------------------------------------------------#

# Load Solver Function
source('Distribution-Solver.R')

# Solve factorization
fx_z <- compute_dist(
    y = z_i,
    a = a_0,
    b = a_1,
    var_y = s_z^2,
    mu_x = g_0,
    var_x = s_x^2
)
fx_z |> lapply(str2lang)

fx_yz <- compute_dist(
    y = y_i,
    a = b_0 + b_2 * z_i,
    b = b_1 + b_3 * z_i,
    var_y = s_x^2,
    mu_x = fx_z$mean,
    var_x = fx_z$variance
)
fx_yz |> lapply(str2lang)
