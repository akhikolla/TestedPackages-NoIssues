## rENA 0.2.1.2

#### Bugs
  * Bug fix in window accumulation

## rENA 0.2.1.1

#### Bugs
  * Removing test checking for output
  * Including 'webshot' as Suggests to fix for issues on CRAN

## rENA 0.2.1.0

#### Features
  * Projections use the centering vector of a provided rotation set. See `ena.make.set()`
  * Faster correlation function: `ena_correlation()`

## rENA 0.2.0.1
  * Fix for updated r-devel to 4.0.0

## rENA 0.2.0.0

#### Features

  * New `ena()` function for easier model generation
  * Updated default model object returned by all methods.  See help for `ena()`
  * Custom S3 methods for removing meta data from data.frames on the ena model.
      - e.g. as.matrix(set$line.weights)

#### Bugs

  * Fixed bug in accumulation code for forward windows
