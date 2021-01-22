# YPInterimTesting 1.0.3
* Fix error occurred when .Random.seed does not exits on .GlobalEnv
* Change the name of the data from 'virtual_data' to 'virtual' 

# YPInterimTesting 1.0.2
* Fix error in setting random seed in the 'ypinterm.default' function

# YPInterimTesting 1.0.1
* Fix error when setting a random.seed

# YPInterimTesting 1.0.0

* Add error messages to `ypinterim.default` to help users detect common errors easily.

* Fixing errors when `critvalue` is supplied in the `ypinterim` function.

* the argument name of `spenfun` in the `ypinterim` function is changed to `spendfun`.

* When the alpha spent is less than 1E-4, the boundaries of that look and the corresponding nominal p-values are now obtained using a (conditional) normal distribution.

* When the value of the observed statistic is greater than `qnorm(1 - 1E-4/2)`, the nominal p-values are now obtained using a normal distribution.

* The default value for `repnum` is set to 1E4.

* The the results from `print` and `summary` of the `ypinterim` class now show the nominal p-values for the boundaries.

* Update the references on the help documents and the DESCRIPTION file as they are now available online.

* Add new examples to the `ypinterim` help documents in order to show how to utilize the argument `critvalue` in the function clearly. 

* Error fixed in the example data: the `virtual_data` has also been changed to have only one `group` vector.

# YPInterimTesting 0.1.0

* The initial version


