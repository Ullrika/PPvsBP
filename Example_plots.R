## Different ways to represent uncertainty at the variable and parameter levels
## This is the code for the plot in the intro of the paper
## Example Body weight

library(dplyr)
library(ggplot2)
library(ggstance)
library(purrr)
library(tidyr)
library(tidybayes)
library(cowplot)
#library(RColorBrewer)
#library(gganimate)

######################################################
## Case 1 - 2D
## Variable level

## Description
# A normal-normal conjugate model
# X_i  ~ N(\mu, \simga2)   (Likelihood)
# \mu|X  ~ N(\mu_0, \sigma2_0)  (Prior distribution)

df <- data.frame(x = rnorm(n = 100, mean = 25, sd = sqrt(15)))
sample_param <- function(data, sigma2, mu_0, sigma2_0){
  n <- length(data)
  
  sigma2_n <- 1/((1 / sigma2_0) + (n/sigma2))
  #mu_n <- (1 / ((1/sigma2_0) + (n/sigma2)))*((mu_0 / sigma2_0) + sum(data)/sigma2)
  mu_n <- sigma2_n *((mu_0 / sigma2_0) + sum(data)/sigma2)
   
  mu <- rnorm(1,mean = mu_n, sd = sqrt(sigma2_n))
 #postpred_data <- rnorm(1, mean = mu_n, sd = sqrt(sigma2 + sigma2_n))
  return(list(mu = mu, sigma = sqrt(sigma2_n), sigma2 = sigma2, mu_n = mu_n, sigma2_n = sigma2_n))
}


#Update the parameter (posterior distribution)  
ndraws <- 20 ## number of spaghetti straws
sample_param_df <- map_dfr(seq(ndraws), 
                           ~sample_param(data = df$x,
                                         sigma2 = 10, mu_0 = 30, sigma2_0 = 5)) %>% 
  tibble::rownames_to_column(var=".iter")


# 2 d plots 
df_ale <- sample_param_df  %>% 
  mutate(q=map2(mu, sigma2, ~tibble(qs=c(0.001, seq(0.01, 0.99, by=0.01), 0.999),
                                   vals=qnorm(qs, mean=.x, sd=.y),
                                   ds=dnorm(vals, mean=.x, sd=.y)))
  ) %>% unnest(q)

p_cdf <- df_ale %>% 
  mutate(grp = .iter) %>% 
  ggplot(aes(group=grp, x=vals, y=qs)) +
  geom_line(data=. %>% select(-.iter), size=1, alpha=0.2, color="grey50")+
  stat_ecdf(data = df, aes(x = x, y = NULL, group=NULL), pad = TRUE) + ## adds the empirical cdf for data
  coord_cartesian(expand = FALSE) +
  labs(
    title = "",
    x = "Body weight",
    y = "cdf"
  ) + 
  theme(axis.title = element_text(size = 30), axis.text = element_text(size = 15))
  p_cdf 

 #####################################################
  ### Different dataframes
  
  df_ale <- sample_param_df  %>% 
    mutate(q=map2(mu, sigma2, ~tibble(qs=c(0.001, seq(0.01, 0.99, by=0.01), 0.999),
                                      vals=qnorm(qs, mean=.x, sd=.y),
                                      ds=dnorm(vals, mean=.x, sd=.y)))
    ) %>% unnest(q)
  
  
  ## Posterior predictive
  
  sample_param_df_2 <- map_dfr(1, 
                               ~sample_param(data = df$x,
                                             sigma2 = 10, mu_0 = 25, sigma2_0 = 5)) 
  
  
  df_ale_2 <- sample_param_df_2  %>% 
    mutate(q=map2(mu_n, sigma2 + sigma2_n, ~tibble(qs=c(0.001, seq(0.01, 0.99, by=0.01), 0.999),
                                                   vals=qnorm(qs, mean=.x, sd=.y),
                                                   ds=dnorm(vals, mean=.x, sd=.y)))
    ) %>% unnest(q)
  
  
  ## Plot 2D and posterior predictive 
  
  p_cdf_together <- df_ale %>% 
    mutate(grp = .iter) %>% 
    ggplot(aes(group=grp, x=vals, y=qs)) +
    geom_line(data=. %>% select(-.iter), size=1, alpha=0.2, color="grey50")+
    stat_ecdf(data = df, aes(x = x, y = NULL, group=NULL, color = "empirical cdf data"), pad = TRUE ) + ## adds the empirical cdf for data
    geom_line(data = df_ale_2, aes(x = vals, y = qs, group = NULL, color = "posterior predictive"), size=1.5, alpha=0.5, ) +  ## adds the cdf of the posterior predictive distribution
    coord_cartesian(expand = FALSE) +
    scale_colour_manual(name = '', values = c("black", 'brown'), labels = c('empirical cdf data', 'posterior predictive')) +
    labs(
      title = "",
      x = "Body weight",
      y = "cdf"
    ) + 
    theme(axis.title = element_text(size = 30), axis.text = element_text(size = 15), 
          legend.text = element_text(size = 15),
          legend.justification =  'bottom', legend.position = c(0.8,0))
  p_cdf_together 
  
  
######################
## Parameter levels (mu)

# 2 d plots 
df_epi <- sample_param_df  %>% 
  mutate(q=map2(mu_n, sigma2_n, ~tibble(qs=c(0.001, seq(0.01, 0.99, by=0.01), 0.999),
                                    vals=qnorm(qs, mean=.x, sd=.y),
                                    ds=dnorm(vals, mean=.x, sd=.y)))
  ) %>% unnest(q)

## pdf plot
p_pdf <- df_epi %>%
  mutate(grp = .iter) %>% 
  ggplot(aes(group=grp, x=vals, y=ds)) +
  geom_line(data=. %>% select(-.iter), size=1, alpha=0.2, color="blue")+
  coord_cartesian(expand = FALSE) +
  labs(
    title = "",
    x = "mu",
    y = "pdf"
  )
p_pdf + theme(axis.title = element_text(size = 30), axis.text = element_text(size = 15))


#####################################################
## Case 2 - Posterior predictive data_pred ~ N(mu_n, sigma2_n + sigma2)
## CDF plot

sample_param_df_2 <- map_dfr(1, 
                           ~sample_param(data = df$x,
                                         sigma2 = 10, mu_0 = 25, sigma2_0 = 5)) 
  
df_ale_2 <- sample_param_df_2  %>% 
  mutate(q=map2(mu_n, sigma2 + sigma2_n, ~tibble(qs=c(0.001, seq(0.01, 0.99, by=0.01), 0.999),
                                   vals=qnorm(qs, mean=.x, sd=.y),
                                   ds=dnorm(vals, mean=.x, sd=.y)))
  ) %>% unnest(q)


p_cdf_2 <- df_ale_2 %>% 
  # mutate(grp = .iter) %>% 
  ggplot(aes(x=vals, y=qs)) +
  geom_line(size=1.2, color="brown") +
  coord_cartesian(expand = FALSE) +
  labs(
    title = "",
    x = "Body weight",
    y = "cdf"
  ) + 
  theme(axis.title = element_text(size = 30), axis.text = element_text(size = 15))
p_cdf_2 


## Parameter level (single point)
# sample_param_df_2 %>% 
#  ggplot(aes(y = '', x = mu)) +
#  geom_point(size = 3, col = 'blue') + 
#  labs(
#    title = "",
#    x = "mu",
#    y = ""
#  ) +
# theme(axis.title = element_text(size = 30), axis.text = element_text(size = 15))

####################################33
### Case 3 - P-boxes
## Variable level 
## Normal([mu_1, mu_2],sigma2)
mu_1 <- 20
mu_2 <- 25

sample_param_df_3_mu0_1 <- data.frame(mu = 20, sigma = sqrt(50))

sample_param_df_3_mu0_2 <- data.frame(mu = 25, sigma = sqrt(50))

df_ale_3_mu0_1 <- sample_param_df_3_mu0_1  %>% 
  mutate(q=map2(mu, sigma, ~tibble(qs=c(0.001, seq(0.01, 0.99, by=0.01), 0.999),
                                   vals=qnorm(qs, mean=.x, sd=.y),
                                   ds=dnorm(vals, mean=.x, sd=.y)))
  ) %>% unnest(q)

df_ale_3_mu0_2 <- sample_param_df_3_mu0_2  %>% 
  mutate(q=map2(mu, sigma, ~tibble(qs=c(0.001, seq(0.01, 0.99, by=0.01), 0.999),
                                   vals=qnorm(qs, mean=.x, sd=.y),
                                   ds=dnorm(vals, mean=.x, sd=.y)))
  ) %>% unnest(q)

graph_bp <- function(lower_points, upper_points){
  
  n_values <- length(lower_points)
  l_points <- rep(0,n_values)
  u_points <- rep(0,n_values)
  
  for(i in 1:n_values){
    l_points[i] <- min(lower_points[i],upper_points[i]) 
    u_points[i] <- max(lower_points[i],upper_points[i])
  }
  
  # data wide format
  data_plot_wide <- data.frame(l_points = sort(l_points), 
                               u_points = sort(u_points), 
                               cdf = c(1:n_values / n_values))
  
  app_data_plot <- data.frame(l_points= c(min(data_plot_wide$l_points),max(data_plot_wide$u_points)),
                              u_points = c(min(data_plot_wide$l_points),max(data_plot_wide$u_points)),
                              cdf=c(0,1))
  
  data_plot <- rbind(data_plot_wide, app_data_plot)
  
  
  # data long format
  data_plot <- gather(data_plot, bound, values, l_points, u_points, factor_key = TRUE)
  
  
  p <- data_plot %>% 
    ggplot(aes(x = values, y = cdf, group = bound, col = bound)) +
    geom_line() +
    scale_color_manual(labels = c('Upper', 'Lower'), values = c('red', 'blue')) +
    guides(color = guide_legend("Bounds")) +
    labs(
      title = "",
      x = "Body weight",
      y = "cdf")
  p + theme(axis.title = element_text(size = 30), axis.text = element_text(size = 15), 
            legend.title = element_text(size = 15), legend.text = element_text(size = 15),
            legend.justification =  'bottom', legend.position = c(0.9,0))
}

graph_bp(lower_points = df_ale_3_mu0_1$vals, upper_points = df_ale_3_mu0_2$vals)

# Parameter level (interval)

sample_param_df_3_mu0 <- data.frame(x = c(sample_param_df_3_mu0_1$mu, sample_param_df_3_mu0_2$mu))

sample_param_df_3_mu0 %>%
  ggplot(aes(y = '', x = x)) +
  geom_line(size = 1.5, col = 'blue') + 
  labs(
    title = "",
    x = "mu",
    y = ""
  ) +
theme(axis.title = element_text(size = 30), axis.text = element_text(size = 15))



 