# ------------------------------------------------------------------------------
# Title: Statistical Rethinking Chaper 6: Haunted DAG
# ------------------------------------------------------------------------------

library(dagitty)
library(brms)
library(tidyverse)
library(tidybayes)

theme_set(theme_minimal())

# simulating selection bias ----------------------------------------------------

set.seed(1914)
n <- 200  # num grant proposals
p <- 0.1  # proportion to select

d <-
  # uncorrelated newsworthiness and trustworthiness
  tibble(newsworthiness  = rnorm(n, mean = 0, sd = 1),
         trustworthiness = rnorm(n, mean = 0, sd = 1)) %>% 
  # total_score
  mutate(total_score = newsworthiness + trustworthiness) %>% 
  # select top 10% of combined scores
  mutate(selected = ifelse(total_score >= quantile(total_score, 1 - p), TRUE, FALSE))

head(d)

d %>% 
  filter(selected == TRUE) %>% 
  select(newsworthiness, trustworthiness) %>% 
  cor()

# we'll need this for the annotation
text <-
  tibble(newsworthiness  = c(2, 1), 
         trustworthiness = c(2.25, -2.5),
         selected        = c(TRUE, FALSE),
         label           = c("selected", "rejected"))

d %>% 
  ggplot(aes(x = newsworthiness, y = trustworthiness, color = selected)) +
  geom_point(aes(shape = selected), alpha = 3/4) +
  geom_text(data = text,
            aes(label = label)) +
  geom_smooth(data = . %>% filter(selected == TRUE),
              method = "lm", fullrange = T,
              color = "lightblue", se = F, size = 1/2) +
  scale_color_manual(values = c("black", "lightblue")) +
  scale_shape_manual(values = c(1, 19)) +
  scale_x_continuous(limits = c(-3, 3.9), expand = c(0, 0)) +
  coord_cartesian(ylim = range(d$trustworthiness)) +
  theme(legend.position = "none")

# multicollinearity ------------------------------------------------------------
# predict height using both legs
n <- 100
set.seed(909)

d <- tibble(
  height = rnorm(n, mean = 10, sd = 2),
  leg_prop = runif(n, min = 0.4, max = 0.5)
  ) %>% 
  mutate(
    leg_left  = leg_prop * height + rnorm(n, mean = 0, sd = 0.02),
    leg_right = leg_prop * height + rnorm(n, mean = 0, sd = 0.02)
    )

# correlation between leg height
d %>%
  select(leg_left:leg_right) %>%
  cor() %>%
  round(digits = 4)

# mcmc model
b6.1 <- brm(
  data = d,
  family = gaussian(),
  formula = height ~ 1 + leg_left + leg_right, # linear model
  prior = c(
    prior(normal( 10 , 100 ), class = Intercept), # vague/bad prior
    prior(normal( 2 , 10), class = b),
    prior( exponential( 1 ), class = sigma)
  ),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 6,
  file = './r_code/fits/b06.01'
)
  
print(b6.1)

mcmc_plot(b6.1, 
          type = "intervals", 
          prob = .5, 
          prob_outer = .95,
          point_est = "mean") +
  labs(title = "The coefficient plot for the two-leg model",
       subtitle = "Holy smokes; look at the widths of those betas!") +
  theme(axis.text.y = element_text(hjust = 0),
        panel.grid.minor = element_blank(),
        strip.text = element_text(hjust = 0)) 

# note in a linear model of both legs, the posterior is a good sum of legs, which
# is reasonable for prediciton problems, but won't help us understand the 
# relationship; fine for prediction; hard for interpretation

pairs(b6.1, pars = parnames(b6.1)[2:3])


post %>% 
  ggplot(aes(x = b_leg_left + b_leg_right, y = 0)) +
  geom_halfeyeh(point_interval = median_qi, 
                fill = "steelblue", .width = .95) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(title = "Sum the multicollinear coefficients",
       subtitle = "Marked by the median and 95% PIs")

# extract posterior
post <- posterior_samples(b6.1)

post %>% 
  ggplot(aes(x = b_leg_left, y = b_leg_right)) +
  geom_point(color = "forestgreen", alpha = 1/10, size = 1/2)


post %>% 
  ggplot(aes(x = b_leg_left + b_leg_right, y = 0)) +
  geom_halfeyeh(point_interval = median_qi, 
                fill = "steelblue", .width = .95) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(title = "Sum the multicollinear coefficients",
       subtitle = "Marked by the median and 95% PIs")

# refitting the model with just one leg
b6.2 <- brm(
  data = d,
  family = gaussian(),
  formula = height ~ 1 + leg_left,
  prior = c(
    prior(normal( 10, 100 ), class = Intercept),
    prior(normal( 2, 10 ), class = b),
    prior(exponential( 1 ), class = sigma)
  ),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 6,
  file = './r_code/b06.02'
)

summary(b6.2)
# for evern unit increase in left leg height, hiehg increases by 2
mcmc_plot(b6.2)

# multicollinarity milk example ------------------------------------------------
# load milk data from rethinking 
data(milk)
d <- milk

# standarize k, f, l
d <-d %>% 
  mutate(kcal_s = standardize(kcal.per.g),
         fat_s = standardize(perc.fat),
         lactose_s = standardize(perc.lactose))
# view correlation
d %>% 
  select(kcal_s, fat_s, lactose_s) %>% 
  pairs(col = 'blue')

# univariate model k | f; 1 sd increase in fat results in 0.83 increase in kcal_s
b6.3 <- brm(
  data = d,
  family = gaussian(),
  formula = kcal_s ~ 1 + fat_s,
  prior = c(
    prior( normal( 0 , 0.2 ), class = Intercept),
    prior( normal( 0 , 0.5 ), class = b),
    prior( exponential( 1 ), class = sigma)
  ),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 6,
  file = './r_code/fits/b06.03'
)

summary(b6.3)

mcmc_plot(b6.3)

posterior_summary(b6.3) %>% round(digits = 2)
# univariate model k | l; 1 sd increase in lactose resutls in 0.9 decrease in kcal_s

b6.4 <- brm(
  data = d,
  family = gaussian(),
  formula = kcal_s ~ 1 + lactose_s,
  prior = c(
    prior( normal( 0 , 0.2 ), class = Intercept),
    prior( normal( 0 , 0.5 ), class = b),
    prior( exponential( 1 ), class = sigma)
  ),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 6,
  file = './r_code/fits/b06.04'
)

summary(b6.4)

mcmc_plot(b6.4)

posterior_summary(b6.4) %>% round(digits = 2)
# add both to same model; can update brm model 6.4
b6.5 <- update(
  object = b6.4,
  newdata = d,
  formula = kcal_s ~ 1 + fat_s + lactose_s,
  seed = 6,
  file = './r_code/fits/b06.05'
)

summary(b6.5)
posterior_summary(b6.5) %>% round(digits = 2)
# page 168 summary
# get posterior samples
post <- posterior_samples(b6.5) %>% 
  mutate(fat_lactose_sum_post = b_fat_s + b_lactose_s) %>% 
  select(fat_lactose_sum_post)


ggplot(post, aes(x=fat_lactose_sum_post)) +
  geom_density(fill = 'lightblue', alpha = 0.5)

# post-treatment bias ----------------------------------------------------------
# mediation

# simulated 
# how many plants would you like?
n <- 100

set.seed(71)
d <- tibble(
  h0 = rnorm(n, mean = 10, sd = 2), 
  treatment = rep(0:1, each = n / 2),
  fungus = rbinom(n, size = 1, prob = .5 - treatment * 0.4),
  h1 = h0 + rnorm(n, mean = 5 - 3 * fungus, sd = 1)
  )

precis(d)

# simulate prior proportion
sim_p <- rlnorm(1e4, 0, 0.25)
precis( data.frame( sim_p ) )

# model m6.6 based on proportion of height at baseline
b6.6 <- brm(
  data = d,
  family = gaussian(),
  formula = h1 ~ 0 + h0, # based on proportion
  prior = c(
    prior(lognormal( 0 , 0.25 ), class = b),
    prior(exponential( 1 ), class = sigma)
  ),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 6,
  file = './r_code/fits/b06.06'
)

summary(b6.6)

# fit mediated path model
b6.7 <- brm(
  data = d,
  family = gaussian(),
  bf(
    h1 ~ h0 * (a + t * treatment + f * fungus),
    a + t + f ~ 1,
    nl = TRUE
  ),
  prior = c(
    prior(lognormal( 0, 0.2 ), nlpar = a),
    prior(normal( 0, 0.5 ), nlpar = t),
    prior(normal( 0, 0.5 ), nlpar = f),
    prior(exponential( 1 ), class = sigma ) 
  ),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 7,
  file = './r_code/fits/b06.07'
  )

summary(b6.7)

# proper model without mediator
b6.8 <- brm(
  data = d,
  family = gaussian(),
  bf(
    h1 ~ h0 * (a + t*treatment),
    a + t ~ 1,
    nl = TRUE
    ),
  prior = c(
    prior(lognormal( 0 , 0.2), nlpar = a),
    prior(normal( 0 , 0.5 ), nlpar = t),
    prior(exponential( 1 ), class = sigma)
  ),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 6,
  file = './r_code/fits/b06.08'
  )

summary(b6.8)

# ended here: might add homework, but i got the gist of it. 

