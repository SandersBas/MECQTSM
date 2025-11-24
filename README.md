# MECQTSM
This repository contains the accompanying toolkit for "_Measurement Error and Counterfactuals in Quantitative Trade and Spatial Models_" by Bas Sanders.

The folder **General algorithm toolkit** implements the high-level approach from Algorithm 1 in the paper.

- The function `GA_toolkit.m` takes as inputs the number of posterior draws to produce, the structural parameter and a vector of noisy data, and outputs the point estimate for gamma and a vector of posterior draws for gamma. It requires auxiliary functions `pi_post.m` and `g.m`.
- The function `pi_post.m` takes as input a vector of noisy data and outputs a posterior draw of the true data.
- The function `g.m` takes as inputs a vector of true data and the structural parameter and outputs a scalar counterfactual prediction.

The folder **Default approach toolkit** implements the default approach from Algorithm 3 in the paper.

- The function `DA_toolkit.m` takes as inputs the number of posterior draws to produce, the structural parameter, a vector of noisy data, a vector of probabilities of true zeros, a vector of probabilities of spurious zeros, a vector of measurement error variances, a vector of prior means, and a vector of prior variances. It outputs the point estimate for gamma and a vector of posterior draws for gamma. It also outputs a report that compares the normalized residuals histogram with standard normal pdf. It requires the auxiliary function "g.m".
- The function `g.m` takes as inputs a vector of true data and the structural parameter and outputs a scalar counterfactual prediction.

The folder **Mirror trade toolkit** uses the mirror trade dataset of Linsi, Burgoon, and Mügge (2023) and allows the researcher to choose countries and years for which they want to estimate the parameters of the prior and measurement error model.

- The function `MT_toolkit.m` takes as inputs a set of countries, a set of years to use for calibration, and a set of years to produce bootstrap draws for. It outputs a vector of noisy data, a vector of probabilities of true zeros, a vector of probabilities of spurious zeros, a vector of measurement error variances, a vector of prior means, and a vector of prior variances. It also outputs a report that computes the adjusted R-squared of the gravity model for the last year, and plots log flows against log distance after partitioning out fixed effects for the last year.
- The script `Example.m` illustrates the use of the function `MT_toolkit.m` for a small set of countries.
