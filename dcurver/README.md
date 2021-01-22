
<!-- README.md is generated from README.Rmd. Please edit that file -->
Davidian curves in R
====================

A Davidian curve defines a seminonparametric density, whose shape and flexibility can be tuned by easy to estimate parameters. Since a special case of a Davidian curve is the standard normal density, Davidian curves can be used for relaxing normality assumption in statistical applications \[1\].

This package provides the density function, the gradient of the loglikelihood and a random generator for Davidian curves:

1.  `ddc(x, phi)`, Davidian curve density function
2.  `rdc(n, phi)`, a random sampler for Davidian curves
3.  `dc_grad(x, phi)`, the gradient function of Davidian curves

The second argument `phi` of these functions is a vector, which contains the parameter(s) of a Davidian curve. The higher the number of parameters, the more flexible the density. The values of the parameres define the shape of the curve. Following \[1, 2\], this package provides support for Davidian curves up to 10 parameters.

References
==========

1.  Zhang, D., & Davidian, M. (2001). Linear mixed models with flexible distributions of random effects for longitudinal data. *Biometrics, 57*(3), 795-802.
2.  Woods, C. M., & Lin, N. (2009). Item response theory with estimation of the latent density using Davidian curves. *Applied Psychological Measurement, 33*(2), 102-117.
