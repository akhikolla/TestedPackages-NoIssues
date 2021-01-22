## EMbC 2.0.3

### Bug fixes

	- Solved error in function *stbc()* when using a solar-covariate due to reverse dependencies with maptools.
	
## EMbC 2.0.2

### Bug fixes

	- The parameter *showClst* is now working properly in function *pkml()* (thanks to Huan Xia for noting it).

## EMbC 2.0.1

### Bug fixes

	- Corrected a bug in the post-smoothing function (thanks to Simona Imperio for noting it). Also, we have added a method in the post-smoothing function for binClstStck instances.

## EMbC 2.0.0

### Major improvements

  - Some core functionalities have been implemented using Rcpp. This results in a significant increase in computational speed when using large datasets.

  - With the aim of enhancing the *EMbC* R-package as a real general purpose algorithm, the dependency with respect to the *move* R-package has been dropped. Thus, potential users coming from domains other than ecology are not forced to install packages that are not desired. However, *Move* objects from the *move* R-package can still be used directly as input data.

### Bug fixes

  - The *bg* parameter of the *sctr* function is working properly now.


## EMbC 1.9.4

### Bug fixes

  - When working with *move* objects, version 1.9.3 was picking GMT timestamps instead of using study local timestamps. This has been changed in version 1.9.4. Of note: the example in the vignette illustrating the use of *move* objects was affected by this error.

### Minor improvements

  - The command *view* was depicting the trajectory without taking into account the proportionality of the axes. This has been corrected in version 1.9.4.

  - In the multivarite case, commands *view()* and *sctr()* depict plots with a light-grey background colour to enhance the visibility of the data points. We added a *bg* parameter to allow changing this default behaviour.
