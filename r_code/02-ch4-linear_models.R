# ------------------------------------------------------------------------------
# Title: Statistical Rethinking Chaper 4; Linear models
# Week 2
# ------------------------------------------------------------------------------

library(tidyverse)
library(rethinking)

# 4.1  Normal by addition ------------------------------------------------------
# Point: Process that adds together random values form same distribution
# convergest to normal

pos <- replicate( 10000, sum (runif(16, -1, 1) ) )

plot(density(pos))

# 4.2 Normal by multiplication -------------------------------------------------
# Point: Multiplying small numbers is approximately the same as normal by addition
growth <- replicate(10000, prod ( 1 + runif(12, 0, 0.1) ) )
dens( growth, norm.comp = TRUE, )

# 4.3 Normal by log-multiplication ---------------------------------------------
# Point: Multiplying even big numbers that are log transformed is approximately
# the same as normal by addition

big <- replicate(10000, prod ( 1 + runif(12, 0, 0.5) ) )

# big plot; not normal
dens(big, norm.comp = TRUE)

# log normal; normal
log.big <- replicate(10000, log( prod ( 1 + runif(12, 0, 1) ) ) )

dens( log.big, norm.comp = TRUE)

# Model definition to Bayes Theorem --------------------------------------------
# water/land problem; see page 78 for bayes formula
w <- 6
n <- 9
# grid approximation
p_grid <- seq(from=0, to=1, length.out=100)
posterior <- dbinom(w, n, p_grid)*dunif(p_grid, 0, 1)
posterior <- posterior/sum(posterior)

plot(posterior, type = 'l')

# Guassian model for height ----------------------------------------------------
data("Howell1")
d <- Howell1

head(d)

# precis summary function
precis(d)

# filter to only adults
d2 <- d[d$age >= 18, ]

# density of height of adults
dens(d2$height)

# check prior
# modeling prior height of 178 cm with sd of 20
curve( dnorm( x, mean = 178, sd = 20), from = 100, to = 250)
# variance prior
curve( dunif( x, 0, 50 ), from = -10, to = 60)


# sample mu and sigma from 10k, mu 178 and sigma 20
sample_mu <- rnorm( 1e4 , mean = 178, sd = 20)
sample_sigma <- runif( 1e4 , min = 0 , max = 50)
# prior
prior_h <- rnorm( 1e4, sample_mu , sample_sigma)

dens(prior_h)

# sample mu 178 and sigma of 100; flatter prior
sample_mu <- rnorm( 1e4 , mean = 178, sd = 100)
# prior
prior_h <- rnorm( 1e4, sample_mu , sample_sigma)

dens(prior_h)


# 4.16 code
# vector of means from 150 to 160 by .101
mu.list <- seq( from = 150, to = 160, length.out = 100)
sigma.list <- seq( from = 7, to = 9, length.out = 100)
# matrix of mu and sigma combinations
post <- expand.grid( mu = mu.list, sigma = sigma.list)

# log likelihood of density of normal distribution function
post$LL <- sapply( 1:nrow(post), function(i) {
  sum( dnorm(d2$height, mean = post$mu[i], sd= post$sigma[i], log = TRUE ) )
  } )

# product? of log likelihood; top part of function?
post$prod <- (
  post$LL + dnorm( post$mu, mean = 178, sd = 20, log = TRUE) +
    dunif( post$sigma, min = 0, max = 50, log = TRUE)
)
# posterior distribution?
post$prob <- exp( post$prod - max(post$prod) )

# countour
contour_xyz( post$mu, post$sigma, post$prob )

# heatmap
image_xyz( post$mu, post$sigma, post$prob )

# resample based on posterior probability
sample.rows <- sample( 1:nrow(post) , size = 1e4, replace = T, prob = post$prob)
sample.mu <- post$mu[ sample.rows ]
sample.sigma <- post$sigma[ sample.rows ]

plot( sample.mu, sample.sigma, cex = 0.5, pch = 16, col = col.alpha(rangi2, 0.1))

dens( sample.mu )
dens( sample.sigma )

PI( sample.mu )
PI( sample.sigma )

# example on with smaller sample sizes, normality can't be assumed
d3 <- sample( d2$height, size = 20)

mu.list <- seq( from = 150, to = 170, length.out = 200)
sigma.list <- seq( from = 4, to = 20, length.out = 200)
# matrix of mu and sigma combinations
post2 <- expand.grid( mu = mu.list, sigma = sigma.list)

# log likelihood of density of normal distribution function
post2$LL <- sapply( 1:nrow(post2), function(i) {
  sum( dnorm(d3 , mean = post2$mu[i], sd= post2$sigma[i], log = TRUE ) )
} )

# product? of log likelihood; top part of function?
post2$prod <- (
  post2$LL + dnorm( post2$mu, mean = 178, sd = 20, log = TRUE) +
    dunif( post2$sigma, min = 0, max = 50, log = TRUE)
)
# posterior distribution?
post2$prob <- exp( post2$prod - max(post2$prod) )
# resample based on posterior probability
sample2.rows <- sample( 1:nrow(post2) , size = 1e4, replace = TRUE, prob = post2$prob)
sample2.mu <- post2$mu[ sample2.rows ]
sample2.sigma <- post2$sigma[ sample2.rows ]

plot( sample2.mu, sample2.sigma, cex = 0.5, pch = 16, col = col.alpha(rangi2, 0.1))
dens( sample2.sigma, norm.comp = TRUE)

# Quadratic approximation ------------------------------------------------------
# function list
flist <- alist(
  height ~ dnorm( mu, sigma ),
  mu ~ dnorm( 178, 20 ),
  sigma ~ dunif( 0 , 50 )
)


m4.1 <- quap( flist, data = d2 )

precis( m4.1 )

# linear model for checking; basically the same. neat
summary(lm(height ~ 1, data = d2))

# starting values for hill climb
start <- list(
  mu = mean(d2$height),
  sigma = sd(d2$height)
)

m4.1 <- quap( flist , data = d2, start = start )

precis( m4.1 )


curve( dnorm( x, 155, 7.7), from = 130, to = 175 )

# narrow prior
m4.2 <- quap(
  alist(
    height ~ dnorm( mu, sigma) ,
    mu ~ dnorm( 178, 0.1) ,
    sigma ~ dunif( 0, 50)
  ) ,
  data = d2 )

# not very good with a narrow prior; sigma (sd) is higher
precis( m4.2 )

curve( dnorm( x , 178, 25 ), 100, 250 )

# correlation matrix
vcov( m4.1 )
# variance (sigma^2)
diag( vcov( m4.1 ) )
# sd in model m4.1
sqrt( diag( vcov( m4.1 ) ) )
precis( m4.1 )

# sampling from multi-dimentional posterior
post <- extract.samples( m4.1 , n = 1e4)

head(post)

precis( post )

# relationship between height and weight

# produces a non-sensicle plot; point is to make beta log normal
set.seed(2971)
n <- 100
a <- rnorm( n, 178, 20)
b <- rnorm( n, 0 , 10)

plot(NULL, xlim = range(d2$weight), ylim = c(-100, 400),
     xlab = 'weight', ylab = 'height')
abline( h = 0, lty = 2)
abline( h = 272, lty = 1, lwd = 0.5)
mtext( "b ~ dnorm(0,10)")
# mean weight
xbar <- mean(d2$weight)
for ( i in 1:n) curve( a[i] + b[i]*(x - xbar) ,
  from=min(d2$weight), to = max(d2$weight), add = TRUE,
  col = col.alpha('black', 0.2) )

# b lognormal
set.seed(2971)
n <- 100
a <- rnorm( n, 178, 20)
# beta log normal
b <- rlnorm( n, 0 , 1)

plot(NULL, xlim = range(d2$weight), ylim = c(-100, 400),
     xlab = 'weight', ylab = 'height')
abline( h = 0, lty = 2)
abline( h = 272, lty = 1, lwd = 0.5)
mtext( "b ~ log(dnorm(0,10))")
# mean weight
xbar <- mean(d2$weight)
for ( i in 1:n) curve( a[i] + b[i]*(x - xbar) ,
                       from=min(d2$weight), to = max(d2$weight), add = TRUE,
                       col = col.alpha('black', 0.2) )


# Fit quap model to find posterior distribution
# mean of weight
xbar <- mean(d2$weight)

# simple regression model; every parameter needs a prior
m4.3 <- quap(
  alist(
    height ~ dnorm( mu , sigma ), # likelihood
    mu <- a + b*(weight - xbar), # linear relationship
    a ~ dnorm( 178, 20 ) , # normal prior of alpha/intercept
    b ~ dlnorm( 0 , 1) ,  # log normal prior of beta (change in height/weight)
    sigma <- dunif( 0, 50) # prior of sigma 
    )
  , data = d2)

precis( m4.3 )

m4.3b <- quap(
  alist(
    height ~ dnorm( mu , sigma ), # likelihood
    mu <- a + exp(log_b)*(weight - xbar), #  exp of b, linear relationship
    a ~ dnorm( 178, 20 ) , # normal prior of alpha/intercept
    log_b ~ dnorm( 0 , 1) ,  #normal prior of beta (change in height/weight)
    sigma <- dunif( 0, 50) # prior of sigma 
  )
  , data = d2)

precis( m4.3b )


# plotting the disribution
plot(height ~ weight, data = d2, col = rangi2)
# add posterior
post <- extract.samples( m4.3 )
a_map <- mean(post$a)
b_map <- mean(post$b)
curve( a_map + b_map*(x - xbar ), add = TRUE)

# posterior based on first 10 samples
N <- 10

dN <- d2[ 1:N , ]
mN <- quap(
  alist(
    height ~ dnorm( mu , sigma ), # likelihood
    mu <- a + b*(weight - xbar), # linear relationship
    a ~ dnorm( 178, 20 ) , # normal prior of alpha/intercept
    b ~ dlnorm( 0 , 1) ,  # log normal prior of beta (change in height/weight)
    sigma <- dunif( 0, 50) # prior of sigma 
  )
  , data = dN
)


# plot 20 samples of the quadratic approximation of the posterior
post <- extract.samples( mN, n = 20)

# display raw data and sample size
plot( dN$weight , dN$height ,
      xlim = range(d2$weight), ylim = range(d2$height),
      col = rangi2, xlab = 'weight', ylab = 'height')
mtext(concat('N = ', N))

for ( i in 1: 20)
  curve( post$a[i] + post$b[i]*(x-mean(dN$weight)) , 
         col = col.alpha('black', 0.3), add = TRUE)

# confidence of mean of line tightens with increase in data

# how does it look from a linear model? pretty similar
ggplot(dN, aes(x=weight, y=height)) +
  geom_point() +
  geom_smooth(method = 'lm', se = TRUE) 

# mean at 50
post <- extract.samples( m4.3 )
mu_at_50 <- post$a + post$b * ( 50 - xbar )

dens( mu_at_50, col = rangi2, lwd = 2, xlab = 'mu|weight=50')

PI( mu_at_50, prob = .89)

# prediction intervals; do this some other time --------------------------------


# Non-linearities --------------------------------------------------------------
# height and weight including children
ggplot(data = d, aes(weight, height )) +
  geom_point() +
  theme_classic()


# define standardize weight
d$weight_s <- (d$weight - mean(d$weight) ) / sd(d$weight)
# squared std weight
d$weight_s2 <- d$weight_s^2

m4.5 <- quap(
  alist(
    height ~ dnorm( mu, sigma ) , # likelihood of height is centered at mu with spread around sigma
    mu <- a + b1*weight_s + b2*weight_s^2 , # polynomial formula
    a ~ dnorm( 178, 20), # prior for a
    b1 ~ dlnorm( 0, 1 ), # log norm prior for std weight
    b2 ~ dnorm( 0, 1 ), # normal prior for sq weight
    sigma ~ dunif( 0, 50)
  ), 
  data = d
)

# estimates
precis( m4.5 )

# check standardized summary of weight
summary(d$weight_s)
# plot of polynomial
weight.seq <- seq( from = -2.2, to = 2, length.out = 30 )
pred_dat <- list( weight_s = weight.seq, weight_s2 = weight.seq^2 )
# mean
mu <- link( m4.5 , data = pred_dat )
# mean for each prediction value
mu.mean <- apply( mu, 2, mean )
mu.PI <- apply( mu , 2, PI, prob = 0.89 )
sim.height <- sim( m4.5, data = pred_dat )
height.PI <- apply( sim.height , 2, PI, prob = 0.89 )


plot(height ~ weight_s, d, col=col.alpha(rangi2, 0.5))
lines(weight.seq, mu.mean)
shade( mu.PI, weight.seq ) # credible interval
shade( height.PI , weight.seq ) # prediction interval

# plot ; standardized weight
ggplot(data = d, aes(standardize(weight), height )) +
  geom_point() +
  geom_line(
    data = tibble(weight.seq, mu.mean ), 
    aes(x = weight.seq, y=mu.mean), 
    color = 'darkblue'
    ) + 
  theme_classic()


# Splines ----------------------------------------------------------------------
# load cherry blossom data
data("cherry_blossoms")
d <- cherry_blossoms
precis( d )

plot( doy ~ year, d )

# complete cases
d2 <- d[ complete.cases(d$doy), ]
num_knots <- 15
# where to place knots
knot_list <- quantile( d2$year, probs = seq(0,1, length.out = num_knots) )
# view knot list
knot_list

# library splines
library(splines)
# b spline based on year
B <- bs(d2$year, knots = knot_list[-c(1, num_knots)], degree = 4, intercept = TRUE)

dim(B)

# plot basis function
plot( NULL, xlim = range(d2$year), ylim = c(0,1) , xlab = 'year', ylab = 'basis')
for ( i in 1:ncol(B) ) lines(d2$year , B[, i])


# spline model of average day of year cherry blossom as a function of year -----
m4.7 <- quap(
  alist(
    D ~ dnorm( mu , sigma ), # probability of day of year based on data
    mu <- a + B %*% w , # model of expected day of year intercept + sum of spline basis 
    a ~ dnorm(100, 15) , # intercept prior
    w ~ dnorm(0, 10) , # weight priors
    sigma ~ dexp(1) # variance prior
  ), data = list( D = d2$doy, B=B ) ,
  start = list( w = rep( 0 , ncol(B) ) )
)

# mean basis weights, but hard to interpret; easier to plot
precis( m4.7 , depth = 2)

# plot basis weights for expectation of doy over years
post <- extract.samples( m4.7 ) # extract samples to posterior
w <- apply( post$w , 2, mean ) # average weight at given year
plot( NULL , xlim = range(d2$year), ylim = c(-6, 6) , 
      xlab = 'year', ylab = 'basis * weight' )
for ( i in 1:ncol(B) ) lines( d2$year , w[i] * B[, i])

# link back to estimates of day of year for given year over time 
mu <- link( m4.7 )
mu_year <- apply(mu, 2, mean)
mu_PI <- apply( mu, 2, PI, 0.97 )
# plot points
plot( d2$year , d2$doy , col = col.alpha(rangi2, 0.3), pch = 16 ,
      xlab = 'Year', ylab = 'Day in Year')
lines(d2$year , mu_year , col = 'red') # average expected value
# plot estimates
shade( mu_PI , d2$year,  col = col.alpha('red', 0.2) ) # 97% credible interval

# Homework ---------------------------------------------------------------------
# predicted heights and 89% compatibility intervals/credible intervals

