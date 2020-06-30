# mapinterval
R package for interpolating interval average data. Only method provided is a mean conserving quadratic spline that happens to be the maximum a-posteriori estimator for a Wiener data generating process.

## Installing
Package can be installad directly from source using devtools:

```
library(devtools)
install_github("c-foschi/mapinterval")
```

Package limSolve is required for tridiagonal system solving and will be automatically installed.

A small running example of how to use it:

```
M= c(0, 1, -1, 2, 4)

library(mlintervals)
mod= interv_fit(M)
mod
# Object of class interv_fit.
# Number of intervals evaluated: 5

plot(mod)
```
![Rplot](https://user-images.githubusercontent.com/39349171/84865661-e56e5a00-b078-11ea-8187-7a4c495c8538.png)
