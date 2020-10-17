# ------------------------------------------------------------------------------
# Title: Statistical Rethinking Chaper 7 Ulysses compass
# ------------------------------------------------------------------------------

library(tidyverse)
library(rethinking)
library(rcartocolor)
library(patchwork)

# probelm with parameters ------------------------------------------------------
# point of this exercise is that more parameters almost always improves fit
# and fit alone will not tell us if a model is useful

d <- tibble(
  species = c("afarensis", "africanus", "habilis", "boisei", 
              "rudolfensis", "ergaster", "sapiens"), 
  brain   = c(438, 452, 612, 521, 752, 871, 1350), 
  mass    = c(37.0, 35.5, 34.5, 41.5, 55.5, 61.0, 53.5)
  ) %>% 
  mutate(
    mass_std = (mass - mean(mass)) / sd(mass),
    brain_std = brain / max(brain)
  )

summary(d)

m7.1 <- quap( 
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)), 
    mu <- a + b * mass_std,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10),
    log_sigma ~ dnorm(0, 1)
  ), 
  data = d)

precis( m7.1 )

# R squared examle
R2_is_bad <- function(quap_fit, seed = 7, ...) {
  set.seed(seed)
  s <- sim(quap_fit, refresh = 0, ...)
  r <- apply(s, 2, mean) - d$brain_std
  1 - var2(r) / var2(d$brain_std)
}

R2_is_bad(m7.1)

# increasingly complex polynomial models 
# quadratic
m7.2 <- quap( 
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)), 
    mu <- a + b[1] * mass_std + b[2] * mass_std^2,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10),
    log_sigma ~ dnorm(0, 1)
  ), 
  data = d, start = list(b = rep(0, 2)))

# cubic
m7.3 <- quap( 
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)), 
    mu <- a + b[1] * mass_std + b[2] * mass_std^2 + b[3] * mass_std^3,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10),
    log_sigma ~ dnorm(0, 1)
  ), 
  data = d, start = list(b = rep(0, 3)))

# fourth-order
m7.4 <- quap( 
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)), 
    mu <- a + b[1] * mass_std + b[2] * mass_std^2 + b[3] * mass_std^3 + b[4] * mass_std^4,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10),
    log_sigma ~ dnorm(0, 1)
  ), 
  data = d, start = list(b = rep(0, 4)))

# fifth-order
m7.5 <- quap( 
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)), 
    mu <- a + b[1] * mass_std + b[2] * mass_std^2 + b[3] * mass_std^3 + b[4] * mass_std^4 + b[5] * mass_std^5,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10),
    log_sigma ~ dnorm(0, 1)
  ), 
  data = d, start = list(b = rep(0, 5)))

# sixth-order
m7.6 <- quap( 
  alist(
    brain_std ~ dnorm(mu, 0.001), 
    mu <- a + b[1] * mass_std + b[2] * mass_std^2 + b[3] * mass_std^3 + b[4] * mass_std^4 + b[5] * mass_std^5 + b[6] * mass_std^6,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10)
  ), 
  data = d, start = list(b = rep(0, 6)))

# plot models
make_figure7.3 <- function(quap_fit, ylim = range(d$brain_std)) {
  
  # compute the R2
  r2 <- R2_is_bad(quap_fit)
  
  # define the mass_seq 
  mass_seq <- seq(from = -2, to = 2, length.out = 100) 
  
  # simulate and wrangle
  link(quap_fit, data = list(mass_std = mass_seq)) %>% 
    data.frame() %>% 
    pivot_longer(everything()) %>% 
    mutate(mass_std = rep(mass_seq, times = 1000)) %>% 
    group_by(mass_std) %>% 
    mean_qi(value, .width = .89) %>% 
    
    # plot!  
    ggplot(aes(x = mass_std, y = value)) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper),
                    color = carto_pal(7, "BurgYl")[7], size = 1/2, 
                    fill = alpha(carto_pal(7, "BurgYl")[6], 1/3)) +
    geom_point(data = d,
               aes(y = brain_std),
               color = carto_pal(7, "BurgYl")[7]) +
    labs(subtitle = bquote(italic(R)^2==.(round(r2, digits = 2))),
         x = "body mass (std)",
         y = "brain volume (std)") +
    coord_cartesian(xlim = c(-1.2, 1.5),
                    ylim = ylim)
}

p1 <- make_figure7.3(m7.1)
p2 <- make_figure7.3(m7.2)
p3 <- make_figure7.3(m7.3)
p4 <- make_figure7.3(m7.4, ylim = c(.25, 1.1))
p5 <- make_figure7.3(m7.5, ylim = c(.1, 1.4))
p6 <- make_figure7.3(m7.6, ylim = c(-0.25, 1.5)) +
  geom_hline(yintercept = 0, color = carto_pal(7, "BurgYl")[2], linetype = 2) 


((p1 | p2) / (p3 | p4) / (p5 | p6)) +
  plot_annotation(title = "Figure7.3. Polynomial linear models of increasing\ndegree for the hominin data.")

