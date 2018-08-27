# aspline
## The Adaptive Splines Regression

This package implements [A-Spline](https://arxiv.org/abs/1808.01770) regression, an adaptive procedure for fitting splines with automatic selection of the knots.
One-dimentional B-Splines of any non-negative degree  can be fitted.
This method uses a penalization approach to compensate overfitting.
The penalty is proportional to the number of knots used to support the regression spline.
Consequently, the fitted splines will have knots only where necessary.

## Illustration

```r

```
