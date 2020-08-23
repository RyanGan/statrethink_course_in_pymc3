# ------------------------------------------------------------------------------
# Title: Statistical Rethinking Chaper 5; Confounding
# Week 3 associated python notebook
# ------------------------------------------------------------------------------

library(tidyverse)
library(rethinking)

# divorce and age --------------------------------------------------------------
data("WaffleDivorce")

d <- WaffleDivorce

head(d)
# standardize values
d$D <- standardize( d$Divorce )
d$M <- standardize( d$Marriage )
d$A <- standardize( d$MedianAgeMarriage )

# sd of age at marriage
sd( d$MedianAgeMarriage )

m5.1 <- quap(
  alist(
    D ~ dnorm( mu, sigma ) , # likelihood
    mu <- a + bA * A, # relationship
    a ~ dnorm( 0 , 0.2 ), # alpha/intercept prior
    bA ~ dnorm( 0, 0.5 ), # beta median age marriage prior
    sigma ~ dexp( 1 ) # variance prior
  ), data = d )


precis( m5.1 )

# simulate priors
set.seed(10)
prior <- extract.prior( m5.1 )
mu <- link( m5.1, post = prior, data = list( A = c(-2, 2) ) )
plot( NULL , xlim = c(-2,2), ylim = c(-2,2) )
for ( i in 1:50 ) lines( c(-2, 2), mu[i, ], col = col.alpha('black', 0.4) )

# posterior predictions
A_seq <- seq( from = -3, to = 3.2 , length.out = 30 ) # sequence of age at marriage
mu <- link( m5.1, data = list(A = A_seq ) )
mu.mean <- apply( mu, 2, mean )
mu.PI <- apply( mu, 2, PI )

# plot priors
plot( D ~ A, data = d, col = rangi2 )
lines( A_seq, mu.mean, lwd = 2)
shade( mu.PI , A_seq )

# 1 sd increase in age at marriage decreases divorce rate
precis( m5.1 )

# linear model pretty comparable in slope
summary(lm(D ~ A, data = d ))

# marriage rate and divorce rate -----------------------------------------------

# modeling divorce rate dependent on marriage rate
m5.2 <- quap(
  alist(
    D ~ dnorm( mu, sigma ) , # likelihood of divorce
    mu <- a + bM * M , # function of mu conditional on marriage rate
    a ~ dnorm( 0 , 0.2 ),  # weakly informed prior on standarized intercept
    bM ~ dnorm( 0 , 0.5), # weakly informed prior on beta
    sigma ~ dexp( 1 )
  ), data = d )

precis( m5.2 )

summary(d$M)
sd( d$M )
# posterior predictions
M_seq <- seq( from = -2, to = 3 , length.out = 30 ) # sequence of age at marriage
mu <- link( m5.2, data = list(M = M_seq ) )
mu.mean <- apply( mu, 2, mean )
mu.PI <- apply( mu, 2, PI )

# plot priors
plot( D ~ M, data = d, col = rangi2 )
lines( M_seq, mu.mean, lwd = 2)
shade( mu.PI , M_seq )

# multiple regression with quadratic approximation -----------------------------

m5.3 <- quap(
  alist(
    D ~ dnorm( mu, sigma ), # likelihood/probability given data
    mu <- a + bM*M + bA*A, # linear model of intercept a + beta marriage rate M + beta age at marriage
    a ~ dnorm( 0 , 0.2 ), # weakly informed prior on intercept
    bM ~ dnorm( 0 , 0.5 ), # prior on beta M
    bA ~ dnorm( 0 , 0.5 ), # prior on beta A
    sigma ~ dexp( 1 ) # prior on sigma
  ), data = d )

# divorce rate is independent of marriage rate when conditioned on median age
# at marriage. divorce rate decreases as median age at marriage increases

precis( m5.3 )

# plot
plot( coeftab(m5.1, m5.2, m5.3 ), par = c("bA", "bM") )

# mar rate and a; as the median age of marriage increases, the marriage
# rate declines. States where median age of marriage is higher also have a lower
# rate of marriage
summary(lm(M ~ A, data = d))

ggplot(d, aes(x = A, y = M)) +
  geom_point()

# diagnostic plots -------------------------------------------------------------
# predictor residual plots

m5.4 <- quap(
  alist(
    M ~ dnorm( mu, sigma ),
    mu <- a + bAM * A,
    a ~ dnorm( 0, 0.2 ),
    bAM ~ dnorm( 0, 0.5 ),
    sigma ~ dexp( 1 )
  ), data = d )

# residuals subtracting the mean from predicted
mu <- link( m5.4 )
mu_mean <- apply( mu, 2, mean )
mu_resid <- d$M - mu_mean

plot(mu_resid, d$D)

# posterior prediction plots
# call link without new data; uses original
mu <- link( m5.3 )

# summarize samples across cases
mu_mean <- apply(mu, 2, mean )
mu_PI <- apply(mu, 2, PI )

# simulate observations; use original data
D_sim <- sim( m5.3, n = 1e4 )
D_PI <- apply( D_sim, 2, PI )

# posterior prediction plot
plot( mu_mean ~ d$D, col = rangi2 , ylim = range(mu_PI) ,
      xlab = 'Observed divorce', ylab = 'Predicted divorce' )
abline( a=0, b=1, lty = 2)
for (i in 1:nrow(d) ) lines( rep(d$D[i], 2) , mu_PI[,i] , col = rangi2 )
#identify( x=d$D , y=mu_mean , labels = d$Loc)


# counterfactual plots; two regressions at once for D <- A -> M -> D dag
m5.3_A <- quap(
 alist(
   # A -> D <- M
   D ~ dnorm( mu, sigma ),
   mu <- a + bA*A + bM*M,
   a ~ dnorm( 0, 0.2 ),
   bA ~ dnorm( 0 , 0.5 ),
   bM ~ dnorm( 0, 0.5 ),
   sigma ~ dexp( 1 ),
   #  A -> M
   M ~ dnorm( mu_M, sigma_M ),
   mu_M <- aM + bAM*A, 
   aM ~ dnorm( 0 , 0.2 ),
   bAM ~ dnorm( 0 , 0.5 ),
   sigma_M ~ dexp( 1 )
   ) , data = d)

# summary
precis( m5.3_A )

# range of A values to predict at
A_seq <- seq( from = -2, to = 2 , length.out = 30)

# simulate  M -> D of joint influence of A on M, and then A and M on D
sim_dat <- data.frame(A = A_seq)

s <- sim( m5.3_A, data = sim_dat, vars = c('M', 'D'))

# plot counterfactual plots of different realizations of Age at marriage
# this includes both the A -> D controlled direct effect and
# the indirect / mediated effect from A -> M -> D
plot( sim_dat$A, colMeans(s$D), ylim=c(-2,2), type = 'l', 
      xlab = 'manipulated A', ylab = 'counterfactual D')
shade( apply(s$D, 2, PI), sim_dat$A )
mtext('Total counterfactual effect of A on D')

plot(sim_dat$A, colMeans(s$M), ylim=c(-2,2), type = 'l',
     xlab='manipulated A', ylab = 'counterfactual M')
shade( apply(s$M, 2, PI), sim_dat$A)
mtext('Counterfactual effect of A -> M')


# calculate relationship of divorce rate
sim2_dat <- data.frame( A = (c(20, 30) - 26.1)/1.24)
s2 <- sim(m5.3_A, data = sim2_dat, vars = c('M', 'D'))
mean( s2$D[,2] - s2$D[,1])

# sim data of A -> D <- M; setting A = 0
sim_dat <- data.frame( M = seq( from = -2, to = 2, length.out = 30), A = 0 )
s <- sim(m5.3_A, data = sim_dat, vars = 'D')  

plot( sim_dat$M, colMeans(s), ylim = c(-2,2), type = 'l', 
      xlab = 'manipuated M', ylab = 'counterfactual D')
shade( apply(s, 2, PI), sim_dat$M )
mtext( "Total counterfactual effect of M on D" )  
  
# simulatej your own counterfactuals
A_seq <- seq( from =  -2, to = 2, length.out = 30)
post <- extract.samples( m5.3_A )
M_sim <- with(
  post, 
  sapply(1:30 , function(i) rnorm( 1e3 , aM + bAM*A_seq[i] , sigma_M )
         )
  )

D_sim <- with(
  post,
  sapply(1:30, function(i) rnorm( 1e3 , a + bA*A_seq[i] + bM*M_sim[,i], sigma))
)
plot(A_seq, colMeans(D_sim), type = 'l')

# masked relationship ----------------------------------------------------------
# primate milk relatonship
data(milk)
d <- milk
str(d)

d$K <- standardize( d$kcal.per.g )
d$N <- standardize( d$neocortex.perc )
d$M <- standardize( log(d$mass) )

detach(package:rethinking, unload = TRUE)
# going to try out brms; but may eventually use stan since it's both python and r
library(brms)

d %>% 
  select(kcal.per.g, mass, neocortex.perc) %>% 
  pairs(col = 'firebrick4')

# standardized
d %>% 
  select(K, M, N) %>% 
  pairs(col = 'firebrick4')

# model
b5.5_draft <- brm(
  data = d,
  family = 'gaussian',
  K ~ 1 + N,
  prior = c(
    prior(normal(0,1), class = Intercept),
    prior(normal(0,1), class = b),
    prior(exponential(1), class = sigma)
  ),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 5,
  sample_prior = TRUE,
  file = './r_code/fits/b05.05_draft'
)

# missing data; complete case
dcc <- d %>% 
  drop_na(K, M, N)

# update fit model
b5.5_draft <- update(b5.5_draft, newdata = dcc, seed = 5)

# simulate and plot 50 prior regression lines
set.seed(5)

prior <- prior_samples(b5.5_draft) %>% 
  sample_n(size = 50) %>% 
  rownames_to_column() %>% 
  expand(nesting(rowname, Intercept, b), N = c(-2,2)) %>% 
  mutate(kcal.per.g_s = Intercept + b * N)

# plot priors; note this is not a good fit
ggplot(prior, aes(x=N, y=kcal.per.g_s)) +
  geom_line(aes(group = rowname), color = 'lightblue', alpha = 0.8) +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "neocortex percent (std)",
       y = "kilocal per g (std)",
       subtitle = "Intercept ~ dnorm(0, 1)\nb ~ dnorm(0, 1)") +
  theme_bw() +
  theme(panel.grid = element_blank()) 

print(b5.5_draft)

# tighten up priors
b5.5 <- brm(
  data = dcc,
  family = 'gaussian',
  formula = K ~ N,
  prior = c(
    prior(normal( 0,0.2 ), class = Intercept),
    prior(normal( 0, 0.5 ), class = b),
    prior(exponential( 1 ), class = sigma)
  ), 
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 5,
  sample_prior = TRUE,
  file = './r_code/fits/b05.05'
)

print(b5.5)

# priors look better
set.seed(5)
prior_samples(b5.5) %>% 
  sample_n(size = 50) %>% 
  rownames_to_column() %>% 
  expand(nesting(rowname, Intercept, b),
         neocortex.perc_s = c(-2, 2)) %>% 
  mutate(kcal.per.g_s = Intercept + b * neocortex.perc_s) %>% 
  
  ggplot(aes(x = neocortex.perc_s, y = kcal.per.g_s)) +
  geom_line(aes(group = rowname),
            color = "firebrick", alpha = .4) +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "neocortex percent (std)",
       y = "kilocal per g (std)",
       subtitle = "Intercept ~ dnorm(0, 0.2)\nb ~ dnorm(0, 0.5)") +
  theme_bw() +
  theme(panel.grid = element_blank()) 

# posterior
nd <- tibble(N = seq(from = -2.5, to = 2, length.out = 30))

fitted(b5.5, 
       newdata = nd,
       probs = c(.025, .975, .25, .75)) %>%
  data.frame() %>%
  bind_cols(nd) %>% 
  
  ggplot(aes(x = N, y = Estimate)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5),
              fill = "firebrick", alpha = 1/5) +
  geom_smooth(aes(ymin = Q25, ymax = Q75),
              stat = "identity",
              fill = "firebrick4", color = "firebrick4", alpha = 1/5, size = 1/2) +
  geom_point(data = dcc, 
             aes(x = N, y = K),
             size = 2, color = "firebrick4") +
  coord_cartesian(xlim = range(dcc$N), 
                  ylim = range(dcc$K)) +
  labs(x = "neocortex percent (std)",
       y = "kilocal per g (std)") +
  theme_bw() +
  theme(panel.grid = element_blank())

# log mass standardized as new predictor
b5.6 <- brm(
  data = dcc,
  family = 'gaussian',
  formula = K ~ M,
  prior = c(
    prior(normal( 0 , 0.2 ), class = Intercept),
    prior(normal( 0 , 0.5 ), class = b),
    prior(exponential( 1 ), class = sigma)
  ),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 5,
  sample_prior = TRUE,
  file = './r_code/fits/b05.06'
)

print(b5.6)

nd <- tibble(M = seq(from = -2.5, to = 2.5, length.out = 30))

fitted(b5.6, 
       newdata = nd,
       probs = c(.025, .975, .25, .75)) %>%
  as_tibble() %>%
  bind_cols(nd) %>% 
  
  ggplot(aes(x = M, y = Estimate)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5),
              fill = "firebrick", alpha = 1/5) +
  geom_smooth(aes(ymin = Q25, ymax = Q75),
              stat = "identity",
              fill = "firebrick4", color = "firebrick4", alpha = 1/5, size = 1/2) +
  geom_point(data = dcc, 
             aes(x = M, y = K),
             size = 2, color = "firebrick4") +
  coord_cartesian(xlim = range(dcc$M), 
                  ylim = range(dcc$K)) +
  labs(x = "log body mass (std)",
       y = "kilocal per g (std)") +
  theme_bw() +
  theme(panel.grid = element_blank())

# multivariate model in kilocals ~ Neurocortex + log(Mass)

m5.7 <- brm(
  data = dcc,
  family = gaussian(),
  formula = K ~ 1 + M + N,
  prior = c(
    prior( normal( 0 , 0.2 ), class = Intercept ),
    prior( normal( 0, 0.5 ), class = b), # interesting, only need to define beta once?
    prior( exponential( 1 ), class = sigma)
  ),
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 5,
  file = './r_code/fits/b05.07'
)
  
print(m5.7)
  
bind_cols(
  posterior_samples(m5.7) %>% 
    transmute(`m5.7_beta[N]` = b_N,
              `m5.7_beta[M]` = b_M)
  ) %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  summarise(mean = mean(value),
            ll   = quantile(value, prob = .025),
            ul   = quantile(value, prob = .975)) %>% 
  separate(name, into = c("fit", "parameter"), sep = "_") %>% 
  # complete(fit, parameter) %>% 
  
  ggplot(aes(x = mean, y = fit, xmin = ll, xmax = ul)) +
  geom_pointrange(color = "firebrick") +
  geom_hline(yintercept = 0, color = "firebrick", alpha = 1/5) +
  ylab(NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "transparent", color = "transparent")) +
  facet_wrap(~parameter, ncol = 1, labeller = label_parsed) 
  
# categorical ------------------------------------------------------------------
data("Howell1")
d <- Howell1

# simulate directly from priors
prior <- tibble(
  mu_female = rnorm(1e4, mean = 178, sd = 20),
  mu_male = mu_female + rnorm(1e4, mean = 0, sd = 10)
)

prior %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  summarise(mean = mean(value),
            sd   = sd(value),
            ll   = quantile(value, prob = .025),
            ul   = quantile(value, prob = .975)) %>% 
  mutate_if(is.double, round, digits = 2) 

prior %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(x = value, fill = name, color = name)) +
  geom_density(size = 3/4, alpha = 2/3) +
  scale_fill_manual(NULL, values = c("firebrick4", "black")) +
  scale_color_manual(NULL, values = c("firebrick4", "black")) +
  scale_y_continuous(NULL, breaks = NULL) +
  xlab("prior predictive distribution for our dummy groups") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(.84, .84))
  
# note variance is wider for males compared to females; other approach is to 
# use index variable
# make index
d <- d %>% 
  mutate(sex = ifelse(male == 1, 2, 1))

d <- d %>% 
  mutate(sex = factor(sex))

d %>% glimpse()

b5.8 <- brm(
  data = d,
  family = 'gaussian',
  formula = height ~ 0 + sex,
  prior = c(
    prior(normal(178, 20), class = b),
    prior(exponential(1), class = sigma)
  ), 
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 5,
  file = './r_code/fits/b05.08'
)

print(b5.8)  
summary(b5.8)

# use milk data
d <- milk %>% 
  mutate(
    clade = factor(clade),
    kcal.per.g_s = (kcal.per.g - mean(kcal.per.g))/sd(kcal.per.g))



b5.9 <- brm(
  data = d,
  family = gaussian(),
  formula = kcal.per.g_s ~ 0 + clade,
  prior = c(
    prior(normal(0, 0.5), class = b),
    prior(exponential(1), class = sigma)
  ), 
  iter = 2000, warmup = 1000, chains = 4, cores = 4,
  seed = 5,
  file = './r_code/fits/b05.09'
)

print(b5.9)
# can use brms function
mcmc_plot( b5.9, pars = 'b')
  