# Psychoco 2025 Talk

### Automatically Derived Full Conditionals for Factored Regression Models with Missing Data and Latent Variables

#### Author
Brian T. Keller

#### Abstract
This talk explores a computational approach for deriving conditional distributions of incomplete observations in a factored regression (i.e., sequential modeling) framework. This method allows full conditional distributions to be solved and implemented in a Gibbs sampler for a wide range of models, including products of latent and incomplete variables, latent centering in single and multilevel models, and dynamical system models with incomplete variables. I have implemented this algorithm in Blimp, a standalone software that estimates complex incomplete data models; however, for this talk, I demonstrate an R implementation of the algorithm and illustrate using functional programming techniques to have R solve the conditional distribution for incomplete predictors inside a Gibbs sampler.

