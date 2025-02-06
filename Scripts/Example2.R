#-------------------------------------------------------------------#
# Automatically Derived Full Conditionals for Factored Regression 
#   Models with Missing Data and Latent Variables
#
# Psychoco 2025
#
# Example2.R
#
# Copyright Brian T. Keller 2025 all rights reserved
# License: GPL-3
#-------------------------------------------------------------------#

# Load Solver Function
source('Distribution-Solver.R')

## Factorization:
#
# f(Y_1 | eta) * f(Y_2 | eta) * f(Y_3 | eta) * f(eta)
#
## Models:
#
# f(Y_1 | eta) -> Y_1 ~ normal(nu_1 + l_1 * eta, s2_e1)
# f(Y_2 | eta) -> Y_2 ~ normal(nu_2 + l_2 * eta, s2_e2)
# f(Y_3 | eta) -> Y_3 ~ normal(nu_3 + l_3 * eta, s2_e3)
# f(eta)       -> eta ~ normal(alpha, psi2)
#
## Calculate distribution for f(eta | Y_1, Y_2, Y_3)
#
# Step 1: f(eta | Y_1)           \propto f(Y_1 | eta) * f(eta)
# Step 2: f(eta | Y_1, Y_2)      \propto f(Y_2 | eta) * f(eta | Y_1)
# Step 3: f(eta | Y_1, Y_2, Y_3) \propto f(Y_3 | eta) * f(eta | Y_1, Y_2)
#

# Step 1: f(eta | Y_1) propto f(Y_1 | eta) * f(eta)
dist <- compute_dist(
    y = Y_1,
    a = nu_1,
    b = l_1,
    var_y = s_e1^2,
    mu_x = alpha,
    var_x = psi2
)

# Step 2:  f(Y_2 | eta) * dist
dist <- compute_dist(
    y = Y_2,
    a = nu_2,
    b = l_2,
    var_y = s_e2^2,
    mu_x = dist$mean,
    var_x = dist$variance
)

# Step 3: f(Y_3 | eta) * dist
dist <- compute_dist(
    y = Y_3,
    a = nu_3,
    b = l_3,
    var_y = s_e3^2,
    mu_x = dist$mean,
    var_x = dist$variance
)

# Print as language
dist |> lapply(str2lang)

# $mean
# (((s_e2^2) * ((s_e1^2) * (psi2))/((psi2) * (l_1) * (l_1) + (s_e1^2)))/(((s_e1^2) * 
#     (psi2))/((psi2) * (l_1) * (l_1) + (s_e1^2)) * (l_2) * (l_2) + 
#     (s_e2^2)) * (l_3) * ((Y_3) - (nu_3)) + (s_e3^2) * (((s_e1^2) * 
#     (psi2))/((psi2) * (l_1) * (l_1) + (s_e1^2)) * (l_2) * ((Y_2) - 
#     (nu_2)) + (s_e2^2) * ((psi2) * (l_1) * ((Y_1) - (nu_1)) + 
#     (s_e1^2) * (alpha))/((psi2) * (l_1) * (l_1) + (s_e1^2)))/(((s_e1^2) * 
#     (psi2))/((psi2) * (l_1) * (l_1) + (s_e1^2)) * (l_2) * (l_2) + 
#     (s_e2^2)))/(((s_e2^2) * ((s_e1^2) * (psi2))/((psi2) * (l_1) * 
#     (l_1) + (s_e1^2)))/(((s_e1^2) * (psi2))/((psi2) * (l_1) * 
#     (l_1) + (s_e1^2)) * (l_2) * (l_2) + (s_e2^2)) * (l_3) * (l_3) + 
#     (s_e3^2))
# 
# $variance
# ((s_e3^2) * ((s_e2^2) * ((s_e1^2) * (psi2))/((psi2) * (l_1) * 
#     (l_1) + (s_e1^2)))/(((s_e1^2) * (psi2))/((psi2) * (l_1) * 
#     (l_1) + (s_e1^2)) * (l_2) * (l_2) + (s_e2^2)))/(((s_e2^2) * 
#     ((s_e1^2) * (psi2))/((psi2) * (l_1) * (l_1) + (s_e1^2)))/(((s_e1^2) * 
#     (psi2))/((psi2) * (l_1) * (l_1) + (s_e1^2)) * (l_2) * (l_2) + 
#     (s_e2^2)) * (l_3) * (l_3) + (s_e3^2))
