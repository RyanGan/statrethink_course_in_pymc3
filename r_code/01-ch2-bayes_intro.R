# ------------------------------------------------------------------------------
# Title: Statistical Rethinking Chaper 2
# Week 1
# ------------------------------------------------------------------------------

library(tidyverse)
library(rethinking)

# Globe Example ----------------------------------------------------------------

# probability of 6 waters in nin toses under equal liklihood of water or land
dbinom(6, size = 9, prob = 0.5) # 0.16

# Grid Approximation Algorithm -------------------------------------------------
# 1. Define the grid. Decide how manuy points to use in estimating the posterior 
#    and make a list parameter values.
# 2. Compute the value of the prior at each parameter value on grid.
# 3. Compute the likelihood of each parameter value.
# 4. Compute the unstandardized posterior at each parameter value by multiplying
#    prior by likelihood
# 5. Standardize the posterior by dividing each value by the sum of all values

# 1. define the grid; using 101 points to match python script
reps <- 1000

p_grid <- seq(from = 0, to = 1, length.out = reps)

# 2. define prior; all 1s
prior <- rep(1, reps)

# 3. compute likelihood
likelihood <- dbinom(6, size = 9, prob = p_grid)

# 4. compute product of likelihood and prior
unstd_posterior <- likelihood * prior

# 5. stadardize the posterior so it sums to 1
posterior <- unstd_posterior / sum(unstd_posterior)

# check that it sums to 1
sum(posterior)
assertthat::are_equal(sum(posterior), 1)

# plot posterior probability of water
ggplot(tibble(p_grid, posterior), aes(x = p_grid, y = posterior)) +
  geom_point() +
  geom_line() +
  xlab('Probability of water grid') +
  ylab('Posterior probability') +
  theme_minimal()

# Replication of Figure 2.5
# replicate different priors in fig 2.5 
prior <- ifelse(p_grid < 0.5, 0, 1)

prior <- exp( -5*abs( p_grid - 0.5 ))

likelihood <- dbinom(6, size = 9, prob = p_grid)

unstd_posterior <- likelihood * prior

posterior <- unstd_posterior / sum(unstd_posterior)

ggplot(tibble(p_grid, posterior), aes(x = p_grid, y = posterior)) +
  geom_line() +
  xlab('Probability of water grid') +
  ylab('Posterior probability') +
  theme_minimal()

# Quadratic approximation ------------------------------------------------------
# uses function in rethinking package; i'll dig in to it sometime to understand
# also consider porting to python for fun at some point
# assumes guassian/normal distribution

# 1. Find the posterior mode
# 2. estimate the curvature of the peak

globe.qa <- quap(
  alist(
    W ~ dbinom( W+L, p), # binomial likelihood
    p ~ dunif(0,1) # uniform prior
  ) ,
   data = list(W=6, L=3)
)

# display summary of quadratic approximation
precis( globe.qa )

# Comparison with beta distribution
# beta distribution with shape w+1 and l+1
w <- 6
l <- 3
# beta given shape
curve( dbeta (x, w+1, l+1), from=0, to=1)
# quadratic approximation
curve( dnorm( x, 0.67, 0.16), lty=2, add = TRUE)

# MCMC Example =----------------------------------------------------------------
# number of samples
n <- 10000
# empty vector
p <- rep(NA, n)
# starting probability
p[1] <- 0.5
# water
w <- 6
# land
l <- 3

# mcmc sampler
for (i in 2:n){
  i <- 2
  p_new <- rnorm( 1, p[i-1] , 0.1)

  if ( p_new < 0 ) p_new <- abs( p_new )
  if ( p_new > 1) p_new <- 2 - p_new
  q0 <- dbinom( w, w+l , p[i-1] )
  q1 <- dbinom( w, w+l, p_new)
  p[i] <- ifelse( runif(1) < q1/q0, p_new , p[i-1])
}

runif(1)
# density of mcmc sampler
dens( p , xlim = c(0,1))
curve( dbeta( x, w+1, l+1 ), lty = 2, add =TRUE)

# Additional Problem Sets ------------------------------------------------------
# 1. Suppose the globe tossing data had turned out to be 8 water in 15 tosses. 
# Construct the posterior distribution, using grid approximation. 
# Use the same flat prior as before.



# 1. define the grid
reps <- 100

p_grid <- seq(from = 0, to = 1, length.out = reps)

# 2. define prior; all 1s
prior <- rep(1, reps)

# 3. compute likelihood
likelihood <- dbinom(8, size = 15, prob = p_grid)

# 4. compute product of likelihood and prior
unstd_posterior <- likelihood * prior

# 5. stadardize the posterior so it sums to 1
posterior <- unstd_posterior / sum(unstd_posterior)

# check that it sums to 1
sum(posterior)
assertthat::are_equal(sum(posterior), 1)

# plot posterior probability of water
ggplot(tibble(p_grid, posterior), aes(x = p_grid, y = posterior)) +
  geom_line() +
  xlab('Probability of water grid') +
  ylab('Posterior probability') +
  theme_minimal()


# 2. Start over in 1, but now use a prior that is zero below p = 0.5 
# and a constant above $p = 0.5$. This corresponds to prior information that 
# a majority of the Earth's surface is water. What difference does the better 
# prior make? If it helps, compare posterior distributions (using both priors) 
# to the true value p = 0.7.

reps <- 100

p_grid <- seq(from = 0, to = 1, length.out = reps)


# 2. define prior; all 1s
prior <- c(rep(0, reps/2), rep(1, reps/2))

# 3. compute likelihood
likelihood <- dbinom(8, size = 15, prob = p_grid)

# 4. compute product of likelihood and prior
unstd_posterior <- likelihood * prior

# 5. stadardize the posterior so it sums to 1
posterior2 <- unstd_posterior / sum(unstd_posterior)

# check that it sums to 1
sum(posterior2)
assertthat::are_equal(sum(posterior2), 1)

# plot uniform posterior and p = 0.5
ggplot() +
  geom_line(
    data = tibble(p_grid, posterior), 
    aes(x = p_grid, y = posterior), 
    color = 'black'
  ) +
  geom_line(
    data = tibble(p_grid, posterior2), 
    aes(x = p_grid, y = posterior2), 
    color = 'blue'
    ) +
  geom_vline(aes(xintercept = 0.7), color = 'red') +
  xlab('Probability of water grid') +
  ylab('Posterior probability') +
  theme_minimal()


# 3. This problem is more open-ended than the others. Feel free to collaborate
# on the solution. Suppose you want to estimate the Earth's proportion of water 
# very precisely. Specifically, you want the 99% percentile interval of the 
# posterior distribution of p to be only 0.05 wide. This means the distance 
# between the upper and lower bound of the interval should be 0.05. How many 
# times will you have to toss the globe to do this? I won't require a precise 
# answer. I'm honestly more interested in your approach.

# set function to take different sample sizes and probability 
grid_function <- function(p, n){
  # random binomial draw based on p and n
  k <- sum(rbinom(n, size = 1, prob = p ))

  p_grid <- seq(from = 0, to = 1, length.out = n)
  
  # 2. define prior; all 1s
  prior <- rep(1, n)
  
  # 3. compute likelihood
  likelihood <- dbinom(k, size = n, prob = p_grid)
  
  # 4. compute product of likelihood and prior
  unstd_posterior <- likelihood * prior
  
  # 5. stadardize the posterior so it sums to 1
  posterior <- unstd_posterior / sum(unstd_posterior)
  
  # sample from probability grid based on posterior
  samp <- sample(p_grid, size = 5000, replace = TRUE, prob = posterior)
  perc99 <- quantile(samp, prob = c(0.05, 0.995))

  df <- tibble( 
    n = n, 
    mean = mean(samp), 
    lower_99 = perc99[1], 
    upper_99 = perc99[2],
    diff = (perc99[2] - perc99[1]),
    less_than_5 = ifelse((perc99[2] - perc99[1]) < 0.05, TRUE, FALSE)
    )
  
  return(df)
}

# iterate through list of sample sizes
n_list <- list(20, 50, 100, 500, 1000, 2000, 3000, 4000, 5000)

posterior_df <- do.call(
  'rbind', 
  lapply(
    n_list, 
    grid_function, 
    p = 0.7 # supply p, n will be provided in the list
    )
)

# print posterior_df
print(posterior_df)

# plot
ggplot(posterior_df, aes(x = as.factor(n), y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower_99, ymax = upper_99), width = 0.2) +
  geom_hline(aes(yintercept = 0.7), color = 'red') +
  geom_text(aes(label = round(diff,2)), hjust = -.5) +
  ylab('Posterior mean and 99% credibility bounds') +
  xlab('n sample size') +
  theme_minimal()



