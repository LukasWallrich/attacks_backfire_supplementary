Simulation to test analysis method
================
Lukas Wallrich
15/05/2022

We propose to use an analysis inspired by segmented regression - but
different from its usual formulation - to assess whether violence
against refugees influences public attitudes.

Below we show a simulation that might help to clarify the model and
helps to establish its robustness. The parameters are: - a: the effect
of attitudes on violence - b: the effect of violence on attitudes -
noise: the variations in attitude change over time, not explained by
co-variates.

<details>
<summary>
Simulation code
</summary>

``` r
if (!require(pacman)) install.packages("pacman")
```

    ## Loading required package: pacman

``` r
pacman::p_load(faux, purrr, ggplot2, dplyr, furrr)
```

``` r
if (!require(pacman)) install.packages("pacman")
```

    ## Loading required package: pacman

``` r
pacman::p_load(faux, purrr, ggplot2, dplyr, furrr)
set.seed(12345)

# Set possible parameter values
params <- expand.grid(list(
  a = c(0, .05, .1, .15, .2, .25, .3),
  b = c(-.3, -.2, -.1, 0, .1, .2, .3),
  noise = c(.3, .5, .7, .9)
))


# Simulate datasets and estimate b, following approach in article

simulate <- function(a, b, noise, runs = 100, N = 500) {
  diffs <- vector("numeric", runs)
  ps <- vector("numeric", runs)
  att_cors <- vector("numeric", runs)

  for (i in 1:runs) {
    data <- rnorm_multi(N, 2, 0, 1, a, varnames = c("AT1", "VT1"))

    data$AT2 <- c(scale(data$AT1 + b * data$VT1 + noise * rnorm(N)))

    data$VT2 <- rnorm_pre(data$AT2, r = a)

    mod <- lm(AT2 ~ VT1 + VT2, data)

    att_cors[i] <- cor(data$AT1, data$AT2)
    diffs[i] <- mod %>%
      coef() %>%
      {
        .[2] - .[3]
      }

    ps[i] <- car::linearHypothesis(mod, "VT1 = VT2")$`Pr(>F)`[2]
  }

  data.frame(a = a, b = b, noise = noise, diff = diffs, p = ps, att_cor = att_cors)
}

permutations <- nrow(params)


library(magrittr)
```

    ## Warning: package 'magrittr' was built under R version 4.1.3

    ## 
    ## Attaching package: 'magrittr'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     set_names

``` r
plan(multisession, workers = 6)

i <- 0
sim_results_null <- params %>%
  filter(b == 0) %T>%
  {
    permutations <<- nrow(.)
  } %>%
  future_pmap_dfr(.options = future_options(seed = TRUE), function(...) {
    x <- data.frame(...)
    if (interactive()) {
      i <<- i + 1
      print(i / permutations)
    }
    simulate(x$a, x$b, x$noise, runs = 5000, N = 1000)
  })
```

    ## Warning: `future_options()` was deprecated in furrr 0.2.0.
    ## Please use `furrr_options()` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

``` r
i <- 0
sim_results_not_null <- params %>%
  filter(b != 0) %T>%
  {
    permutations <<- nrow(.)
  } %>%
  future_pmap_dfr(.options = future_options(seed = TRUE), function(...) {
    x <- data.frame(...)
    if (interactive()) {
      i <<- i + 1
      print(i / permutations)
    }
    simulate(x$a, x$b, x$noise, runs = 5000, N = 1000)
  })
```

</details>

# Type I error rates

When b is set to 0, i.e. when no effect of attacks on attitudes is
assumed, the simulation can be used to estimate potential Type I error
rates. This shows that the test is generally too liberal and that the
Type I error rate can be substantially inflated if attitudes are
unstable over time or if the correlation between municipality-level
attitudes and subsequent attacks is large.

Halving the alpha level to .025 leads to a Type I error rate below the
desired 5% for the most likely parameter combinations.

<details>
<summary>
Plot code
</summary>

``` r
p <- sim_results_null %>%
  mutate(att_cor = round(att_cor / 5, 2) * 5) %>%
  group_by(a, b, att_cor) %>%
  filter(n() > 100) %>%
  summarise(sig_share = mean(p < .025), below_alpha = sig_share < .05) %>%
  filter(b == 0) %>%
  ggplot(aes(x = a, y = att_cor, fill = sig_share * 100, col = below_alpha)) +
  geom_tile(
    lwd = 1.5,
    linetype = 1
  ) +
  scale_fill_gradient(low = "green", high = "red") +
  labs(
    title = "Simulation of Type I error rates", subtitle = "with alpha = .025", fill = "Type I error rate (%)",
    col = " ... below 5%?", y = "Correlation between attitudes over time"
  )
```

    ## `summarise()` has grouped output by 'a', 'b'. You can override using the
    ## `.groups` argument.

</details>

# Type II error rates

Note that we cannot plausibly estimate the power of the multi-level
regression, since we have no basis to estimate interclass correlations.
It will be somewhere between a linear model with 2400 participants (the
number of participants) and one with 140 participants (the number of
clusters). To get a sense of power, we simulated a sample size of 500.

<details>
<summary>
Plot code
</summary>

``` r
p <- sim_results_not_null %>%
  mutate(att_cor = round(att_cor / 5, 2) * 5) %>%
  group_by(a, b, att_cor) %>%
  filter(n() > 100) %>%
  summarise(type_2 = mean(p > .025), above_20 = type_2 < .20) %>%
  filter(!b == 0) %>%
  ggplot(aes(x = a, y = att_cor, fill = type_2 * 100, col = above_20)) +
  geom_tile(
    lwd = 1.5,
    linetype = 1
  ) +
  scale_fill_gradient(low = "green", high = "red") +
  labs(
    title = "Simulation of Type II error rates", subtitle = "with alpha = .025", fill = "Type II error rate (%)",
    col = " ... below 20%?", y = "Correlation between attitudes over time"
  ) +
  facet_wrap(vars(b))
```

    ## `summarise()` has grouped output by 'a', 'b'. You can override using the
    ## `.groups` argument.

</details>

## Bias in estimates

Overall, the estimates for b tend to be biased towards 0
(i.e. underestimated), yet the bias is small unless a and b are similar
in value or attitudes are unstable over time. The figures below provide
detail.
<details>
<summary>
Plot code
</summary>

``` r
estimates_summarised <- sim_results_not_null %>%
  mutate(att_cor = round(att_cor / 5, 2) * 5) %>%
  group_by(a = factor(a), b, att_cor) %>%
  filter(n() > 100) %>%
  summarise(est = mean(diff), .groups = "drop")

p <- estimates_summarised %>%
  ggplot(aes(x = att_cor, y = est, fill = a)) +
  geom_col(position = position_dodge()) +
  facet_wrap(vars(b)) +
  geom_hline(aes(yintercept = b)) +
  labs(title = "Estimates for b under different conditions", 
       subtitle = "Actual value of b is shown in panel title and black line", 
       x = "Correlation between attitudes over time", y = "Estimate for b")
```

</details>
<details>
<summary>
Plot code
</summary>

``` r
p <- estimates_summarised %>%
  group_by(a, b) %>%
  summarise(est = mean(est), .groups = "drop") %>%
  mutate(bias = est / b - 1) %>%
  ggplot(aes(x = a, y = b, fill = bias * 100)) +
  geom_tile(lwd = 1.5, linetype = 1) +
  scale_fill_gradient2(low = "pink", mid = "green", high = "red") +
  labs(title = "Bias in estimates for b depending on a", subtitle = "Averaged over all correlations between attitudes over time", 
       x = "a", y = "b", fill = "Bias in %", 
       caption = "Positive bias (orange) indicates that b is overestimated, \n while negative bias (pink) indicates that it is underestimated")
```

</details>
<details>
<summary>
Plot code
</summary>

``` r
p <- estimates_summarised %>%
  group_by(b, att_cor) %>%
  summarise(est = mean(est), .groups = "drop") %>%
  mutate(bias = est / b - 1) %>%
  ggplot(aes(x = att_cor, y = b, fill = bias * 100)) +
  geom_tile(lwd = 1.5, linetype = 1) +
  scale_fill_gradient2(low = "pink", mid = "green", high = "red") +
  labs(title = "Bias in estimates for b depending on a", subtitle = "Averaged over all values for a", 
       x = "Correlation between attitudes over time", y = "b", fill = "Bias in %", 
       caption = "Positive bias (orange) indicates that b is overestimated, \n while negative bias (pink) indicates that it is underestimated")
```

</details>
