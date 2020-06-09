
library(rstanarm)
library(bayesplot)
library(broom)
library(tidybayes)
library(modelr)
library(sjPlot)
dat <- data


stan_1 <- stan_glmer(Sensitivity/9 ~ Year_i + (1|Site/Plot/Tree) + (1|Year_i),
                      data = dat, family = betar,
                      chains = 4, cores = 4, adapt_delta = 0.99)

stan_1.5 <- stan_glmer(rich ~ I(Year-1996) + (1|Year_i),
                     data = data, family = poisson,
                     chains = 4, cores = 4, adapt_delta = 0.99)

stan_2 <- stan_glmer(rich ~ I(Year-1996) + (1|Site/Plot) + (1|Year_i),
                     data = data, family = poisson,
                     chains = 4, cores = 4, adapt_delta = 0.99)


stan_3 <- stan_glmer(rich ~ I(Year-1996) + (1|Year/Site/Plot),
                     data = data, family = poisson,
                     chains = 4, cores = 4, adapt_delta = 0.99)


fit <- stan_1
summary(fit)

#plot(fit, "trace")
pp_check(fit,plotfun = "stat", stat = "mean")
pp_check(fit, nreps = 200)

#plot model predictions and actuals data

(model_fit <- data %>% 
    data_grid(Year = seq_range(Year, n = 101))%>%
    add_predicted_draws(fit) %>%
    ggplot(aes(x = Year, y = rich)) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),
                    alpha = 1/2, colour = "black")
  +
    geom_point(data = data, colour = "darkseagreen4", size = 3) +
    scale_fill_brewer(palette = "Greys"))


data %>% group_by(Site) %>% 
  data_grid(Year_i) %>%
  add_fitted_draws(fit) %>%
  ggplot(aes(x = .value, y = Year, colour = Site)) +
  stat_pointintervalh(.width = c(.66, .95))

launch_shinystan(fit)

y_pred = posterior_predict(fit)

# show group-level (varying) parameters and group by parameter
posterior_vs_prior(pvp, pars = "varying",
                   group_by_parameter = TRUE, color_by = "vs")

# group by parameter and allow axis scales to vary across facets
posterior_vs_prior(pvp, regex_pars = "period",
                   group_by_parameter = TRUE, color_by = "none",
                   facet_args = list(scales = "free"))


newdata = dat1 %>% cbind(t(y_pred)) %>% gather(key = "Rep", value = "Value",
                                               -n_nh4, -rich.y)

ggplot(dat1, aes(Value, x = N.y)) + geom_violin(color = "blue",
                                                fill = "blue", alpha = 0.5) + geom_violin(data = dat1, aes(y = n_nh4,
                                                                                                           x = N.y), fill = "red", color = "red", alpha = 0.5)


######plotting regression lines ##
# Extract the (post-warmup) posterior draws
posterior1 <- as.matrix(fit)
colnames(posterior1)
means1 <- colMeans(posterior1)

# Take random 100 posterior draws of intercept and slope
# 100 isn't special, but enough to show uncertainty without
# making the graph unreadable
betas <- posterior1[sample(nrow(posterior1), 200), 1:2]


# Plot regression lines implied by the betas
blues <- color_scheme_get("brightblue")
mod1p1 <- ggplot(dat1, aes(x = n_nh4, y = rich.y)) +
  geom_point(color = "gray30") +
  geom_abline(
    intercept = betas[, 1], 
    slope = betas[, 2],
    color = blues[[2]], 
    size = 0.15, 
    alpha = 0.5
  ) +
  geom_abline(
    intercept = means1[1], 
    slope = means1[2],
    size = 1.25, 
    color = blues[[6]]
  ) +
  ylim(0, 10)

plot(mod1p1)