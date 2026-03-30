## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----fig.height = 4, fig.width = 6--------------------------------------------
set.seed(1)
n <- 50
num_groups <- 10
group <- rep(1:num_groups, each = 5) |> factor()
group_effects <- rnorm(length(group), 0, 0.3)
x <- runif(n, 0, 20)
coef_x <- 0.6
se <- runif(n, 0.1, 0.3)
fx <- 1.5 * sin(2 * pi * x / 12)
y <- 10 + coef_x * x + fx + 
  group_effects[group] +
  fx + rnorm(n, 0, 1) + rnorm(n, 0, se)
data <- data.frame(y, x, se, group)
plot(y ~ x, data = data, col = group, pch = 19, main = "Simulated data (color = group)")

## -----------------------------------------------------------------------------
library(splinemixmeta)
smm <- splinemixmeta( mgcv::s(x, bs = "cr"), y ~ x, 
                      se = se, manual_fixed = TRUE, 
                      data = data, random = ~ 1 | group)

## -----------------------------------------------------------------------------
summary(smm)

## -----------------------------------------------------------------------------
pred_spline_only <- predict(smm, include_smooths = TRUE, include_REs = FALSE, include_residuals = FALSE, type = "outcome")
pred_spline_and_groups <- predict(smm, include_smooths = TRUE, include_REs = TRUE, include_residuals = FALSE, type = "outcome")
# look at both together
head(cbind(pred_spline_only, pred_spline_and_groups))
# look at only fixed effect, which could be done with `predict`
# but here we illustrate a direct call to `blup`
blup(smm, level=0, vcov = TRUE, se = TRUE) |> head()

## ----fig.height = 4, fig.width = 6--------------------------------------------
plot(smm, ylab = "y", title = "Spline meta-regression fit")

