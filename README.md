# mlintervals
R package for interpolating interval average data. Only method provided is a mean conservng quadratic spline that happens to be maximum likelihood estimator for a Winer data generating process.

## Installing
Package can be installad directly from source using devtools:

```
library(devtools)
install_github("c-foschi/mlintervals")
```

Package limSolve is required for tridiagonal system solving and will be automatically installed.

A small running example of how to use it:

```
M= c(0, 1, -1, 2, 4)

library(mlinterval)
mod= interv_fit(M)
mod
# Object of class interv_fit.
# Number of intervals evaluated: 5

plot(mod)
```
