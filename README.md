# aspline
## What are Adaptive Splines?

This package implements [A-Spline](https://arxiv.org/abs/1808.01770) regression, an adaptive procedure for fitting splines with automatic selection of the knots.
One-dimentional B-Splines of any non-negative degree  can be fitted.
This method uses a penalization approach to compensate overfitting.
The penalty is proportional to the number of knots used to support the regression spline.
Consequently, the fitted splines will have knots only where necessary and the fitted model is sparser than comparable penalized spline regressions (for instance, P-splines).

## Illustration

Below is an illustration of the A-spline procedure using the data [helmet](https://github.com/goepp/aspline/blob/master/data/helmet.rda).
The thick line represents the fitted spline. 
The thin line represents the B-spline basis decomposition of the fitted curve.
```r
library(aspline)
library(tidyverse)
library(splines2)
data(helmet)
x <- helmet$x
y <- helmet$y
k <- 40
knots <- seq(min(x), max(x), length = k + 2)[-c(1, k + 2)]
pen <- 10 ^ seq(-4, 4, 0.25)
x_seq <- seq(min(x), max(x), length = 1000)
aridge <- aridge_solver(x, y, knots, pen)
a_fit <- lm(y ~ bSpline(x, knots = aridge$knots_sel[[which.min(aridge$ebic)]]))
X_seq <- bSpline(x_seq, knots = aridge$knots_sel[[which.min(aridge$ebic)]], intercept = TRUE)
a_basis <- (X_seq %*% diag(coef(a_fit))) %>%
  as.data.frame() %>%
  mutate(x = x_seq) %>%
  reshape2::melt(id.vars = "x", variable.name = "spline_n", value.name = "y") %>%
  as_tibble() %>%
  filter(y != 0)
a_predict <- data_frame(x = x_seq, pred = predict(a_fit, data.frame(x = x_seq)))
ggplot() +
  geom_point(data = helmet, aes(x, y), shape = 1) +
  geom_line(data = a_predict, aes(x, pred), size = 0.5) +
  geom_line(data = a_basis, aes(x, y, group = spline_n), linetype = 1, size = 0.1) +
  theme(legend.position = "none") +
  ylab("") +
  xlab("")
```
![alt text][helmet_p_spline]

For the sake of comparision, we display here the estimated P-spline with the same data.
The thin line also represent the B-spline basis decomposition.
```r 
p_fit <- mgcv::gam(y ~ s(x, bs = "ps", k = length(knots) + 3 + 1, m = c(3, 2)))
X <- bSpline(x_seq, knots = knots, intercept = TRUE)
p_basis <- (X %*% diag(coef(p_fit))) %>%
  as.data.frame() %>%
  mutate(x = x_seq) %>%
  reshape2::melt(id.vars = "x", variable.name = "spline_n", value.name = "y") %>%
  as_tibble() %>%
  filter(y != 0)
p_predict <- data_frame(x = x_seq, pred = predict(p_fit, data.frame(x = x_seq)))
ggplot() +
  geom_point(data = helmet, aes(x, y), shape = 1) +
  geom_line(data = p_predict, aes(x, pred), size = 0.5) +
  geom_line(data = p_basis, aes(x, y, group = spline_n), linetype = 1, size = 0.1) +
  theme(legend.position = "none") +
  ylab("") + xlab("")
```
![alt text][helmet_a_spline]

The stricking difference between the two methods is that A-spline fits a far sparser model than P-Spline.

[helmet_a_spline]: https://github.com/goepp/aspline/vignettes/helmet_a_spline.pdf
[helmet_p_spline]: https://github.com/goepp/aspline/vignettes/helmet_p_spline.pdf
