## Code prepared by Ullrika and Dmytro

library(dplyr)
library(ggplot2)
library(ggstance)
library(purrr)
library(tidyr)
library(tidybayes)
library(cowplot)
library(RColorBrewer)
library(gganimate)

## Plotting a subjective probablity distribution
df.par <- data.frame(ef = rnorm(1000))

df.par %>%
  ggplot(aes(y = "",  x = ef)) +
  stat_intervalh()+
  scale_color_brewer()+
  labs(
    title = "Precise probability",
    x = "Exceedence frequency",
    y = ""
  )

df.par %>%
  ggplot(aes(y = "",  x = ef)) +
  stat_pointintervalh() + 
  scale_color_brewer()+
  labs(
    title = "Precise probability",
    x = "Exceedence frequency",
    y = ""
  )

df.par %>%
  ggplot(aes(y = "",  x = ef)) +
  stat_gradientintervalh() +
  scale_color_brewer() +
  labs(
    title = "Precise probability",
    x = "Exceedence frequency",
    y = ""
  )

df.par %>%
  ggplot(aes(y = "",  x = ef)) +
  geom_halfeyeh() +
  coord_cartesian(expand = FALSE) +
  scale_color_brewer() +
  labs(
    title = "Precise probability",
    x = "Exceedence frequency",
    y = "pdf for epistemic uncertainty"
  )
  
## Plotting a two-dimensional probability distribution 
## made up data and the function I use for sampling from the posterior
df = data.frame(x = rnorm(10,5,2))
sample_param <- function(data, mu0 = x_mu0, v0 = 1, alpha0 = 1, beta0 = 1){
  #quantify uncertainty in parameters with normal gamma conjugate model
  n = length(data)
  s2 = var(data)
  m = mean(data) 
  mu_n = ( v0 * mu0 + n * m ) / (v0 + n)
  v_n =  v0 + n
  alpha_n = alpha0 + n/2
  beta_n = beta0 + 1/2*(s2*(n - 1) + (n*v0/(v0 + n) * (m - mu0)^2))
  precision = rgamma(1, shape = alpha_n, rate = beta_n)  # variance
  sigma_n = (1/sqrt(precision))
  return(list(mu = rnorm(1,mu_n,sigma_n/v_n), sigma = sigma_n))
}

#####
#Generate data for plotting
ndraws <- 20 ## number of spaghetti straws
sample_param_df <- map_dfr(seq(ndraws), 
                           ~sample_param(data = df$x,
                                         mu0 = 0, v0 = 1, alpha0 = 1, beta0 = 1)) %>% 
  tibble::rownames_to_column(var=".iter")

# 2 d plots
df_ale <- sample_param_df  %>% 
  mutate(q=map2(mu, sigma, ~tibble(qs=c(0.001, seq(0.01, 0.99, by=0.01), 0.999),
                                   vals=qnorm(qs, mean=.x, sd=.y),
                                   ds=dnorm(vals, mean=.x, sd=.y)))
  ) %>% unnest(q)

p_pdf <- df_ale %>%
  mutate(grp = .iter) %>% 
  ggplot(aes(group=grp, x=vals, y=ds)) +
  geom_line(data=. %>% select(-.iter), size=1, alpha=0.2, color="grey50")+
  coord_cartesian(expand = FALSE) +
  labs(
    title = "The spaghetti plot",
    x = "X",
    y = "pdf for aleatory uncertainty"
  )
p_pdf

p_cdf <- df_ale %>% 
  mutate(grp = .iter) %>% 
  ggplot(aes(group=grp, x=vals, y=qs)) +
  geom_line(data=. %>% select(-.iter), size=1, alpha=0.2, color="grey50")+
  stat_ecdf(data = df, aes(x = x, y = NULL, group=NULL), pad = TRUE) + ## adds the empirical cdf for data
  coord_cartesian(expand = FALSE) +
  labs(
    title = "The spaghetti plot",
    x = "X",
    y = "cdf for aleatory uncertainty"
  )
p_cdf

#####
## Plot 2dim with animations
p_pdf_anime <- df_ale %>% 
  mutate(grp = .iter) %>% 
  ggplot(aes(group=grp, x=vals, y=ds)) +
  geom_line(data=. %>% select(-.iter), size=1, alpha=0.2, color="grey50")+
  geom_line(size=1, color="firebrick")+
  coord_cartesian(expand = FALSE) +
  transition_manual(.iter) +
  labs(
    title = "The spaghetti plot",
    x = "X",
    y = "pdf for aleatory uncertainty"
  )

p_pdf_anime

p_cdf_anime <- df_ale %>% 
  mutate(grp = .iter) %>% 
  ggplot(aes(group=grp, x=vals, y=qs)) +
  geom_line(data=. %>% select(-.iter), size=1, alpha=0.2, color="grey50")+
  geom_line(size=1, color="firebrick")+
  coord_cartesian(expand = FALSE) +
  transition_manual(.iter) +
  labs(
    title = "The spaghetti plot",
    x = "X",
    y = "pdf for aleatory uncertainty"
  )

p_cdf_anime