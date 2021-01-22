
## elmNNRcpp ( Extreme Learning Machine )
<br>

The *elmNNRcpp* package is a reimplementation of *elmNN* using *RcppArmadillo* after the [*elmNN* package was archived](https://CRAN.R-project.org/package=elmNN). Based on the documentation of the *elmNN* it consists of,
*"Training and predict functions for SLFN ( Single Hidden-layer Feedforward Neural Networks ) using the ELM algorithm. The ELM algorithm differs from the traditional gradient-based algorithms for very short training times ( it doesn't need any iterative tuning, this makes learning time very fast ) and there is no need to set any other parameters like learning rate, momentum, epochs, etc."*. More details can be found in the package Documentation, Vignette and [blog-post](http://mlampros.github.io/2018/07/05/the_extreme_learning_machine_package/).
<br><br>

To install the package from CRAN use, 

```R

install.packages("elmNNRcpp")


```
<br>

and to download the latest version from Github use the *install_github* function of the devtools package,
<br><br>

```R

remotes::install_github('mlampros/elmNNRcpp')


```
<br>

Use the following link to report bugs/issues,
<br><br>

[https://github.com/mlampros/elmNNRcpp/issues](https://github.com/mlampros/elmNNRcpp/issues)

