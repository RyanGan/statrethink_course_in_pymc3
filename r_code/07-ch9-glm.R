# ------------------------------------------------------------------------------
# GLMs
# ------------------------------------------------------------------------------

# library
library(rethinking)
library(tidybayes)
library(rstan)
library(tidyverse)

# chimpanzees
data("chimpanzees")
d <- chimpanzees

# look at variables
?chimpanzees

# create index variable to model diferent conditions
d$treatment <- 1 + d$prosoc_left + 2*d$condition
unique(d$treatment)

xtabs(~ treatment + prosoc_left + condition, d)

unique(chimpanzees$actor)

# notes: quadratics don't always work as well as mcmc for glms

# flat prior for a simple binomial model using quadratic approximation
m11.1 <- quap(
  alist(
    pulled_left ~ dbinom( 1, p ) , 
    logit(p) <- a ,
    a ~ dnorm( 0 , 10)
  ), data = d)

precis( m11.1 )

# sample from the prior
set.seed(1999)
prior <- extract.prior( m11.1 , n = 1e4)

# convert to outcome scale using inverse link funciton of logit
p <- inv_logit( prior$a )

dens( p, adj = 0.1 )

# key lesson here; flat prior in the logit space is not a flat prior in the 
# outcome space!

# second model with a prior of 1.5 instead of 10
m11.1a <- quap(
  alist(
    pulled_left ~ dbinom( 1, p ) , 
    logit(p) <- a ,
    a ~ dnorm( 0 , 1.5)
  ), data = d)

precis( m11.1a )

# sample from the prior
set.seed(1999)
prior <- extract.prior( m11.1a , n = 1e4)

# convert to outcome scale using inverse link funciton of logit
p_a <- inv_logit( prior$a )

dens( p_a, adj = 0.1 )

# Find priors for beta; using a flat prior to make point on why they are bad in 
# glms
m11.2 <- quap(
  alist(
    pulled_left ~ dbinom( 1, p ) , 
    logit(p) <- a + b[treatment] , 
    a ~ dnorm( 0, 1.5 ),
    b[treatment] ~ dnorm( 0, 10 )
  ), data = d
)

precis(m11.2, depth = 2)

# extract priors
prior <- extract.prior( m11.2, n = 1e4 )

p <- sapply( 1:4, function(k) inv_logit(prior$a + prior$b[, k] ) )

dens( abs(p[,1] - p[,2]), adj = 0.1 )

# repeated with tighter prior for beta

m11.3 <- quap(
  alist(
    pulled_left ~ dbinom( 1, p ) , 
    logit(p) <- a + b[treatment] , 
    a ~ dnorm( 0, 1.5 ),
    b[treatment] ~ dnorm( 0, 0.5 )
  ), data = d
)

# extract priors
prior <- extract.prior( m11.3 , n = 1e4 )

p <- sapply( 1:4, function(k) inv_logit(prior$a + prior$b[, k] ) )

dens( abs(p[,1] - p[,2]), adj = 0.1 )
# mean difference
mean( abs( p[,1] - p[,2]))

# stan version ----
# trimmed data list
dat_list <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  treatment = d$treatment
)

m11.4 <- ulam(
  alist(
    pulled_left ~ dbinom( 1, p ),
    logit(p) <- a[actor] + b[treatment] , 
    a[actor] ~ dnorm( 0 , 1.5 ),
    b[treatment] ~ dnorm( 0 , 0.5 )
  ), data = dat_list, chains = 4, log_lik = TRUE
)

precis( m11.4, depth = 2 )

# estimate individual chimps preference for pulling left lever
post <- extract.samples( m11.4 )
p_left <- inv_logit( post$a )

plot( precis( as.data.frame( p_left ) ), xlim = c(0,1))

# estimate treatment
labs <- c('r/n', 'l/n', 'r/p', 'l/p')
plot(precis(m11.4, depth = 2, pars = 'b' ), labels = labs)
# no real preference for pulling left lever based on these treatment types

# contrast by partner types
diffs <- list(
  db13 = post$b[,1] - post$b[,3],
  db24 = post$b[,2] - post$b[,4]
)

plot(precis(diffs))
# look at stan code
stancode( m11.4 )

"
data{
  int pulled_left[504];
  int treatment[504];
  int actor[504];
}
parameters{
  vector[7] a;
  vector[4] b;
}
model{
  vector[504] p;
  b ~ normal( 0 , 0.5 );
  a ~ normal( 0 , 1.5 );
  for ( i in 1:504 ) {
    p[i] = a[actor[i]] + b[treatment[i]];
    p[i] = inv_logit(p[i]);
  }
  pulled_left ~ binomial( 1 , p );
}
generated quantities{
  vector[504] log_lik;
  vector[504] p;
  for ( i in 1:504 ) {
    p[i] = a[actor[i]] + b[treatment[i]];
    p[i] = inv_logit(p[i]);
  }
  for ( i in 1:504 ) log_lik[i] = binomial_lpmf( pulled_left[i] | 1 , p[i] );
}
"

# posterior prediction check for each treatment
pl <- by( d$pulled_left, list( d$actor , d$treatment ), mean)
pl[1,]

# see plot in book. too much code for me to write

# simpler model without interaction
d$side <- d$ prosoc_left + 1 # right 1, left 2
d$cond <- d$condition + 1 # no partner 1 , partner 2

dat_list2 <- list(
  pulled_left = d$pulled_left,
  actor = d$actor, 
  side = d$side,
  cond = d$cond
)

m11.5 <- ulam(
  alist(
    pulled_left ~ dbinom( 1 , p ),
    logit(p) <- a[actor] + bs[side] + bc[cond] ,
    a[actor] ~ dnorm( 0 , 1.5 ),
    bs[side] ~ dnorm( 0 , 0.5 ),
    bc[cond] ~ dnorm( 0, 0.5 )
  ), data = dat_list2, chains = 4, log_lik = TRUE
)

# compare two models
compare( m11.5, m11.4, func = PSIS)

# going to practice writing stan code; more helpful if i want to write in python
# aggregated count model of the same chimp data

agg_d <- d %>% 
  mutate(
    treatment = factor(1 + prosoc_left + 2 * condition),
    actor = factor(actor), 
    side = prosoc_left + 1,
    cond = condition + 1
  ) %>% 
  # aggregate
  group_by(treatment, actor, side, cond) %>% 
  summarize(left_pulls = sum(pulled_left))

head(agg_d)

# tidy bayes set up list of stan data
stan_data <- compose_data(agg_d)

# aggregate stan program
stan_program <- "
// data block where data elements are defined
data {
  int n; // n group
  int n_actor; // number of unique chimpanzee ids
  int n_treatment; // number of unique treatments groups
  int treatment[n]; // define treatments
  int actor[n]; // define actors
  int left_pulls[n]; // define events of left pulls
}
// define parameters
parameters {
  real a[n_actor];
  real b[n_treatment];
}
// model
model {
  vector[n] p; // define p of length of n
  a ~ normal( 0 , 1.5 ); // alpha parameter for intercept
  b ~ normal( 0 , 0.5 ); // beta parameter
  // random intercept model
  for ( i in 1:n ) {
    p[i] = a[actor[i]] + b[treatment[i]];
    p[i] = inv_logit(p[i]);
  }
  left_pulls ~ binomial( 18 , p ); // likelihood binomial over 18 trials with p based on i
}
"
# run stat model
m11.6 <- stan(model_code = stan_program, data = stan_data)

summary(m11.6, pars = 'b')
traceplot(m11.6)

datplot <- m11.6 %>% 
  spread_draws(b[treatment]) %>% 
  mean_qi() %>% 
  mutate(treatment = c('R/N', 'L/N', 'R/P', 'L/P')[treatment])

ggplot() +
  geom_pointrange(data = datplot, aes(b, treatment, xmin = .lower, xmax = .upper)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme_minimal()

# aggregate graduate admin admissions ----
# load berkley data
data("UCBadmit") 

ucb_admit <- UCBadmit %>% 
    mutate(
      applicant_gender = ifelse(applicant.gender == 'male', 1, 2),
      ratio = admit / applications
      ) %>% 
  select(-applicant.gender)


# to stan data
stan_data <- compose_data(ucb_admit)


stan_program_admit <- "
data {
  int n; // number
  int admit[n]; // vector of admissions
  int applications[n]; // vector of number of applications
  int applicant_gender[n]; // vector of applicant gender
}
// model parameter of alpha only
parameters {
  real a[2]; // i think the a of vector length 2 is for n and p
}
// calculated parameter
transformed parameters {
  vector[n] p; // p 
  for ( i in 1:n ) {
  // caculate probability of admission
    p[i] = a[applicant_gender[i]]; 
    p[i] = inv_logit( p[i] );
  }
}
// model
model {
  a ~ normal( 0 , 1.5 ); // prior for alpha intercept
  for (i in 1:n ) {
    admit[i] ~ binomial(applications[i] , p[i] );
  }
}
"

# intercept only model
m11.7 <- stan( model_code = stan_program_admit , data = stan_data )

# parameters for male 1 female 2 all chains
summary(m11.7, pars = c('a'))$summary
?spread_draws # seems like the similar function to extract draws

# extract posteriods
# using rstand to extract columns
post <- rstan::extract(m11.7, pars = 'a')
diff_a <- post$a[,1] - post$a[,2]
diff_p <- inv_logit(post$a[,1]) - inv_logit(post$a[,2])

# males are 60% to 100% more likely to be admitted relative to females
# absolute difference is 12 to 16% higher
precis(list(diff_a = diff_a , diff_p = diff_p))

# differences in departments likely explain overall difference; new model

head(ucb_admit)

unique(ucb_admit$dept)
stan_program_admit_dept <- "
data {
  int n; // number
  int admit[n]; // vector of admissions
  int applications[n]; // vector of number of applications
  int applicant_gender[n]; // vector of applicant gender
  // adding in dept
  int n_dept; // number of departments
  int dept[n]; // vector of departments 
}
// model parameter of alpha only
parameters {
  real a[2]; // i think the a of vector length 2 is for n and p
  real b[n_dept]; // beta for each department
}
// calculated parameter
transformed parameters {
  vector[n] p; // p 
  for ( i in 1:n ) {
  // caculate probability of admission
    p[i] = a[applicant_gender[i]] + b[dept[i]]; 
    p[i] = inv_logit( p[i] );
  }
}
// model
model {
  a ~ normal( 0 , 1.5 ); // prior for alpha intercept gender
  b ~ normal( 0 , 0.5 ); // prior for beta department
  for (i in 1:n ) {
    admit[i] ~ binomial(applications[i] , p[i] );
  }
}
"

# intercept only model
m11.8 <- stan( model_code = stan_program_admit_dept , data = stan_data )

summary(m11.8, pars = c('a', 'b'))$summary

# calculate new contrast in gender, accounting for department

# extract posteriods
# using rstand to extract columns
post <- rstan::extract(m11.8, pars = 'a')
diff_a <- post$a[,1] - post$a[,2]
diff_p <- inv_logit(post$a[,1]) - inv_logit(post$a[,2])

# males are 60% to 100% more likely to be admitted relative to females
# absolute difference is 12 to 16% higher
precis(list(diff_a = diff_a , diff_p = diff_p))


# Poisson models ------
# poisson is special case of binomial when p is very small and N is very large
# then mean Np and variance Np(1-p) are ~ the same

# Oceanic tools example
data("Kline")
d <- Kline
# log population
d$P <- as.vector(scale( log(d$population) ))
d$contact_id <- ifelse( d$contact == 'high', 2, 1)

# stan data
stan_data <- compose_data(d)

# intercept only model
stan_program <- "
data {
  int n;
  int total_tools[n];
}
parameters {
  real a; // declare alpha
}
model {
  vector[n] lambda; // make lambda vector of length n
  a ~ normal(3, 0.5); // intercept only
  for (i in 1:n) {
    lambda[i] = exp(a);
  }
  total_tools ~ poisson(lambda);  
}
"
m11.9 <- stan(model_code = stan_program, data = stan_data)

summary(m11.9)$summary


# interaction model
stan_program <- "
data {
  int n;
  int n_contact; // number of unique contact
  int contact[n]; // contact observation
  real P[n]; // log population observation
  real Pnew[100]; // not sure what this is? I think simulated populatio?
  int total_tools[n]; // observation of tool counts
}
parameters {
  real a[n_contact]; // declare alpha
  real b[n_contact]; // declare beta
}
model {
  vector[n] lambda; // make lambda vector of length n
  // relationship of lambda
  for (i in 1:n) {
    lambda[i] = exp(a[contact[i]] + b[contact[i]] * P[i]);
  }
  a ~ normal(3, 0.5); // intercept
  b ~ normal(0, 0.2); // beta
  total_tools ~ poisson(lambda);  
}
generated quantities {
  real yhat[100, 2]; 
  for ( i in 1:100 ) {
    for ( j in 1:2 ) {
      yhat[i,j] = exp(a[j] + b[j] * Pnew[i]);
    }
  }
}
"
# add Pnew to stan data
stan_data$Pnew <- seq(min(stan_data$P), max(stan_data$P), length.out = 100)

m11.10 <- stan(model_code = stan_program, data = stan_data)

summary(m11.10, pars = c('a', 'b'))$summary


datplot <- m11.10 %>%
  gather_draws(yhat[idx, contact]) %>%
  median_qi() %>%
  left_join(tibble(idx = 1:100, 'Population' = stan_data$Pnew), by = 'idx') %>%
  mutate(contact = ifelse(contact == 1, 'Low', 'High')) %>%
  rename(Tools = .value)

ggplot(datplot, aes(Population, Tools, linetype = contact, ymax = .upper, ymin = .lower, fill = contact)) +
  geom_ribbon(alpha = .2) +
  geom_line() +
  xlab('log population (std)')

# scientific model

stan_program <- '
data {
  int n;
  int n_contact;
  int contact[n];
  real population[n];
  int total_tools[n];
}
parameters {
  real a[n_contact];
  real<lower=0> b[n_contact];
  real<lower=0> g;
}
model {
  vector[n] lambda;
  for (i in 1:n) {
    lambda[i] = exp(a[contact[i]]) * population[i]^b[contact[i]] / g;
  }
  a ~ normal(1, 1);
  b ~ exponential(1);
  g ~ exponential(1);
  total_tools ~ poisson(lambda);
}
'

m11.11 <- stan(model_code = stan_program, data = stan_data)

m11.11

# negative binomial / gamma poisson -----------
num_days <- 30
y <- rpois( num_days, 1.5 )

num_weeks <- 4
y_new <- rpois( num_weeks, 0.5 * 7 )

y_all <- c( y, y_new )
exposure <- c( rep( 1, 30 ), rep( 7, 4) )
monastery <- c(rep( 0, 30 ) , rep( 1, 4 ))
d <- data.frame( y = y_all, days = exposure, monastery = monastery)
# log days offset
d$log_days <- log( d$days )

# compile stan data
stan_data <- compose_data( d )

stan_data

stan_program <- "
data {
  int <lower=1> n; // lower bound at least 1
  real <lower=0> log_days[n]; // offset
  int <lower=0, upper=1> monastery[n]; // monastery indicator
  int y[n];
}
parameters {
  real a;
  real b;
}
model {
  real lambda[n]; // define lambda/rate
  // functional relationship
  for ( i in 1:n ) {
    lambda[i] = exp( a + b * monastery[i] + log_days[i] );
  }
  // priors
  y ~ poisson( lambda ); 
  a ~ normal( 0 , 1 );
  b ~ normal( 0, 1 );
}
generated quantities {
  real Old;
  real New;
  Old = exp( a );
  New = exp( a +  b );
}
"
m11.12 <- stan(model_code = stan_program, data = stan_data)

m11.12

datplot <- m11.12 %>% gather_draws(Old, New)

ggplot(datplot, aes(.value, .variable)) + 
  stat_halfeye() +
  labs(x = 'Posterior distribution of monastic productivity', y = '')

# cool

# multinomial logistic/multiclass problem ---- 
# simulate career choices among 500 individuals
N <- 500             # number of individuals
income <- c(1,2,5)   # expected income of each career
score <- 0.5*income  # scores for each career, based on income

# next line converts scores to probabilities
p <- softmax(score[1],score[2],score[3])

# now simulate choice
# outcome career holds event type values, not counts
career <- rep(NA,N)  # empty vector of choices for each individual

# sample chosen career for each individual
set.seed(34302)
for ( i in 1:N ) career[i] <- sample( 1:3 , size=1 , prob=p )
stan_data <- list( N=N , K=3 , career=career , career_income=income )

stan_program <- "
data{
    int N; // number of individuals
    int K; // number of possible careers
    int career[N]; // outcome
    vector[K] career_income;
}
parameters{
    vector[K-1] a; // intercepts
    real<lower=0> b; // association of income with choice
}
model{
    vector[K] p;
    vector[K] s;
    a ~ normal( 0 , 1 );
    b ~ normal( 0 , 0.5 );
    s[1] = a[1] + b*career_income[1];
    s[2] = a[2] + b*career_income[2];
    s[3] = 0; // pivot
    p = softmax( s );
    career ~ categorical( p );
}
"

m11.13 <- stan(model_code = stan_program, data = stan_data)

traceplot(m11.13)
m11.13

post <- extract.samples( m11.13 )

# set up logit scores
s1 <- with( post, a[, 1] + b*income[1] )
s2_orig <- with( post , a[,2] * b*income[2] )
s2_new <- with( post, a[,2])

# compute probability for ogirinal and counterfactual
p_orig <- sapply( 1:length(post$b) , function(i) softmax( c(s1[i], s2_orig[i], 0)))
p_new <- sapply( 1:length(post$b) , function(i) softmax( c(s1[i], s2_new[i], 0)))

p_diff <- p_new[2,] - p_orig[2,]

precis(p_diff)



## R code 11.59
N <- 500
# simulate family incomes for each individual
family_income <- runif(N)
# assign a unique coefficient for each type of event
b <- c(-2,0,2)
career <- rep(NA,N)  # empty vector of choices for each individual
for ( i in 1:N ) {
  score <- 0.5*(1:3) + b*family_income[i]
  p <- softmax(score[1],score[2],score[3])
  career[i] <- sample( 1:3 , size=1 , prob=p )
}
stan_data <- list( N=N , K=3 , career=career , family_income=family_income )

stan_program <- "
data{
    int N; // number of observations
    int K; // number of outcome values
    int career[N]; // outcome
    real family_income[N];
}
parameters{
    vector[K-1] a; // intercepts
    vector[K-1] b; // coefficients on family income
}
model{
    vector[K] p;
    vector[K] s;
    a ~ normal(0,1.5);
    b ~ normal(0,1);
    for ( i in 1:N ) {
        for ( j in 1:(K-1) ) s[j] = a[j] + b[j]*family_income[i];
        s[K] = 0; // the pivot
        p = softmax( s );
        career[i] ~ categorical( p );
    }
}
"

m11.14 <- stan(model_code = stan_program, data = stan_data)

precis(m11.14, 2)

# multinomial poisson ------
data(UCBadmit)

UCBadmit <- UCBadmit %>% 
  rename(rejection = reject) %>%
  janitor::clean_names() %>%
  mutate(ratio = admit / applications)


stan_data <- compose_data(UCBadmit)
stan_program <- '
data {
  int n;
  int admit[n];
  int applications[n];
}
parameters {
  real a;
}
transformed parameters {
  real p;
  p = inv_logit(a);
}
model {
  for (i in 1:n) {
    admit[i] ~ binomial(applications[i], p);
  }
}
'
m_binom <- stan(model_code = stan_program, data = stan_data)

stan_program <- '
data {
  int n;
  int admit[n];
  int rejection[n];
  int applications[n];
}
parameters {
  real lambda[2];
  real a[2];
}
model {
  for (i in 1:n) {
    admit[i] ~ poisson(exp(a[1]));
    rejection[i] ~ poisson(exp(a[2]));
  }
}
generated quantities {
  real p;
  p = exp(a[1]) / (exp(a[1]) + exp(a[2]));
}
'
m_pois <- stan(model_code = stan_program, data = stan_data)

summary(m_binom, pars = 'p')$summary

summary(m_pois, pars = 'p')$summary
