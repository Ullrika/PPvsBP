##########
# Code for paper comparing precise to bounded probability 
# Ivette Raices Cruz, Matthias C M Troffaes, Ullrika Sahlin
# 17 april 2020

########################
library(SHELF)
library(HDInterval)
library(dfoptim)
library(ggplot2)
library(dplyr)
library(tidyr)

####################################################
## Load data 
load('data_assessment.Rdata')

###################################################
## Parametric inference quantifying epistemic uncertainty by bounded probability

## Functions 
{
  ## https://en.wikipedia.org/wiki/Conjugate_prior
  
  ## Normal-Gamma conjugate model (Bayesian inference). 
  ## sufficient statistics (sample mean and sample variance) for the normal likelihood are known for log-aluminium concentration and log-consumption 
  
  update_normal_gamma <- function(data, mu0 = 0, v0 = 5, alpha0 = 1, beta0 = 1, suff_stat = TRUE){
    
    if(suff_stat == TRUE){
      n <- data[['sample_size']]
      m <- data[['sample_mean']]
      S2 <- data[['sample_sd']]^2  # sample variance
    }
    else{
     n <- length(data) 
     m <- mean(data)
     S2 <- var(data)
    }
    
    v_n <-  v0 + n
    alpha_n <- alpha0 + n/2
    
    mu_n <- (v0 * mu0 + n * m ) / (v0 + n)
    beta_n <- beta0 + 1/2*S2*(n - 1) + (n*v0/(v0 + n) * (m - mu0)^2)
      
    output <- list(data = data, param = list(prior = list(mu = mu0, v = v0, alpha = alpha0, beta = beta0),
                                            posterior = list(mu = mu_n, v = v_n , alpha = alpha_n, beta = beta_n)))
      
    return(output)
  }
  
  ## Bernoulli-Beta conjugate model (Bayesian inference). 
  
  update_bernoulli_beta <- function(data, alpha0 = 1, beta0 = 1){
    n_non_consumer <- data[['non_consumer_sample_size']]
    n_consumer <- data[['consumer_sample_size']]
    
    alpha_n <- alpha0 + n_consumer
    beta_n <- beta0 + n_non_consumer
    
    output <- list(data = data, param = list(prior = list(alpha = alpha0, beta = beta0),
                                            posterior = list(alpha = alpha_n, beta = beta_n)))
    
    return(output)
  }

  ################################################################################
  ## post refers to the parameters of the posterior distribution -- normal-gamma 
  ## niter_ale refers to the number of generated samples  
  ## percentile_ale indicates if the assessment is done in a high consumer child by 95 or on all children by 0
  
  propagate_normal_gamma <- function(niter_ale, post, percentile_ale){
  
    precision <- rgamma(1, shape = post$param$posterior$alpha, rate = post$param$posterior$beta)  # precision
    sigma_n <- (1/sqrt(precision))  # standard deviation
    mu <- rnorm(1,post$param$posterior$mu, sigma_n / sqrt(post$param$posterior$v))
    if(percentile_ale == 0){
    gen_sample <- rlnorm(niter_ale, meanlog = mu, sdlog = sigma_n)
    }
    else{
    gen_sample <- rep(qlnorm(percentile_ale/100, meanlog = mu, sdlog = sigma_n), niter_ale)
    }
    
    output <- list(gen_sample = gen_sample)
  
    return(output)
  }  

  ################################################################################
  ## Precise Probability EKE 
  ## data from EKE are taken from the report
  ## fitdist is a function in the SHELF R-package
  
  # change_cocoa_dist <- SHELF::fitdist(vals = change.cons$vals, probs = change.cons$probs/100, lower = -31.5, upper = 21.5)
  # plotfit(change_cocoa_dist)
  
  ## The fitted distribution is beta(1.39, 1.33)
  quantify_uncertainty_pp_change_eke <- function(vals, probs){
    change_cocoa_dist <- SHELF::fitdist(vals = vals, probs = probs, lower = min(vals) - 0.5 , upper = max(vals) + 0.5 )
    best_fit <- toString(change_cocoa_dist$best.fitting$best.fit)
    if(best_fit == "normal"){
      mean <-  change_cocoa_dist$Normal$mean
      sd <-  change_cocoa_dist$Normal$sd
      param <-  c(mean = mean, sd = sd)
    }
    if(best_fit == "t"){
      location <-  change_cocoa_dist$Student.t$location 
      scale <-  change_cocoa_dist$Student.t$scale
      df <-  change_cocoa_dist$Student.t$df
      param <-  c(location = location, scale = scale, df = df)
    }
    if(best_fit == "gamma"){
      shape <-  change_cocoa_dist$Gamma$shape
      rate  <-  change_cocoa_dist$Gamma$rate
      param <-  c(shape = shape, rate = rate)
    }
    if(best_fit == "lognormal"){
      mean.log.X  <-  change_cocoa_dist$Log.normal$mean.log.X 
      sd.log.X <-  change_cocoa_dist$Log.normal$sd.log.X
      param <-  c(mean.log.X = mean.log.X, sd.log.X = sd.log.X)
    }
    if(best_fit == "logt"){
      location.log.X <-  change_cocoa_dist$Log.Student.t$location.log.X
      scale.log.X <-  change_cocoa_dist$Log.Student.t$scale.log.X
      df.log.X <-  change_cocoa_dist$Log.Student.t$df.log.X
      param <-  c(location.log.X = location.log.X, scale.log.X = scale.log.X,  df.log.X =  df.log.X)
    }
    if(best_fit == "beta"){
      shape1 <-  change_cocoa_dist$Beta$shape1
      shape2 <-  change_cocoa_dist$Beta$shape2
      param <-  c(shape1 = shape1, shape2 = shape2)
    }
    output <-  list(vals = vals, best_fit = best_fit, param = param)
    return(output)
  }
  
  ################################################################################
  ## Combination of uncertainty

  combine_uncertainty <- function(gen_data_concentration, gen_data_consumption, 
                                  gen_data_EKE, threshold, niter_ale){
    
    gen_data_concentration <- matrix(unlist(gen_data_concentration), ncol = 7, nrow = niter_ale)
    gen_data_consumption <- matrix(unlist(gen_data_consumption), ncol = 7, nrow = niter_ale)
    
    weekly_intake <- rowSums(gen_data_EKE * (gen_data_consumption * 0.007) * gen_data_concentration)
    
    prob_exceed_wi <- mean(weekly_intake > threshold)
    
    return(prob_exceed_wi = prob_exceed_wi)
    
  }
  
  ###############################################################################
  ## Final assessment model -- Precise Probability 
  #**US to IR: This is the uncertainty analysis on the assessment, I 
  #suggest to put everything that is fixed inside the function, and clarify what is changeable
  #I.e. niter_ale, niter_epi, threshold, percentile_ale

  unc_analysis_assessment <- function(niter_ale = 1000, niter_epi = 1000, 
                                      threshold = 1, percentile_ale = 0,
                                      data_concentration, data_consumption, gen_data_EKE, consumers_info_sample_size,
                                      concentration_mu0 = 3.5, concentration_v0 = 5, concentration_alpha0 = 1, 
                                      concentration_beta0 = 1, suff_stat_concentration = TRUE,
                                      consumption_mu0 = 0, consumption_v0 = 5, consumption_alpha0 = 1, 
                                      consumption_beta0 = 1, suff_stat_consumption = TRUE,
                                      consumption_event_alpha0 = 1, consumption_event_beta0 = 1){
    
    nr_products <-  length(data_concentration)
    prob_consumption <- parameters_consumption <- parameters_concentration <- vector('list', nr_products)
    
    # Probability of a child i consumes chocolate product k
    param_consumption <-  lapply(consumers_info_sample_size, update_bernoulli_beta, alpha0 = consumption_event_alpha0, beta0 = consumption_event_beta0)
    
    for(k in 1:7){
      prob_consumption[[k]] <-  rbeta(1,shape1 = param_consumption[[k]]$param$posterior$alpha, shape2 = param_consumption[[k]]$param$posterior$beta)
    }
    
    prob_exceed <- rep(0, niter_epi)
    
    post_concentration <-  lapply(data_concentration, update_normal_gamma, mu0 = concentration_mu0,
                                  v0 = concentration_v0, alpha0 = concentration_alpha0, 
                                  beta0 = concentration_beta0, suff_stat = suff_stat_concentration)
    
    post_consumption <-  lapply(data_consumption, update_normal_gamma, mu0 = consumption_mu0, 
                                v0 = consumption_v0, alpha0 = consumption_alpha0, 
                                beta0 = consumption_beta0, suff_stat = suff_stat_consumption)
    for(j in 1:nr_products){
      parameters_concentration[[j]] <-  post_concentration[[j]]$param
      
      parameters_consumption[[j]] <-  post_consumption[[j]]$param
    }
    
    for(i in 1:niter_epi){
      
      gen_data_concentration <-  lapply(post_concentration, propagate_normal_gamma, niter_ale = niter_ale, percentile_ale = 0)
      
      gen_data_consumption <-  lapply(post_consumption, propagate_normal_gamma, niter_ale = niter_ale, percentile_ale = percentile_ale)
      
      prob_exceed[[i]] <- combine_uncertainty(gen_data_concentration =  gen_data_concentration, gen_data_consumption = gen_data_consumption, 
                                            gen_data_EKE = gen_data_EKE, threshold =  threshold, niter_ale = niter_ale)
      
    }
    
    expected_prob_exceed <- mean(prob_exceed)
    hdi_prob_exceed <- hdi(prob_exceed, credMass = 0.95) # Highest (Posterior) Density Interval
    
    return(list(prob_consumption_event = prob_consumption,
                parameters_concentration = parameters_concentration,
                parameters_consumption = parameters_consumption,
                prob_exceed = prob_exceed, 
                expected_prob_exceed = expected_prob_exceed,
                hdi_prob_exceed = hdi_prob_exceed))
  }
  
  
 ###################################################################################
  ## Figures  
  
  graph_plot_pp <- function(assessment_output){
    
    prob_cdf <- c(1:length(assessment_output$prob_exceed))/(length(assessment_output$prob_exceed))
    data_plot <- data.frame(prob_exceed = sort(assessment_output$prob_exceed), prob_cdf = prob_cdf)
    
    data_plot %>% 
      ggplot(mapping = aes(y = prob_cdf,  x = prob_exceed)) +
      geom_line() +
      labs(
        title = "Uncertainty",
        x = "Frequency of exceeding TWI",
        y = "cdf")
   }

}
#########################################################################################################################################
#########################################################################################################################################
## Precise probability

## Generate one eke sample 
  fit <- quantify_uncertainty_pp_change_eke(vals = data_assessment$change_cons$vals, probs = data_assessment$change_cons$probs/100)
  gen_eke <- rbeta(1,shape1 = fit$param[1], shape2 = fit$param[2])

## Final assessment
  
  TWI_pp <-  unc_analysis_assessment(niter_ale = 10000, niter_epi = 10000, threshold = 1, percentile_ale = 0,
                                     data_concentration = data_assessment$log_concentration_ss_data, 
                                     data_consumption = data_assessment$log_consumption_ss_data, gen_data_EKE = gen_eke,
                                     consumers_info_sample_size = data_assessment$consumers_info_sample_size,
                                     concentration_mu0 = 3.5, concentration_v0 = 5, concentration_alpha0 = 1, 
                                     concentration_beta0 = 1, suff_stat_concentration = TRUE,
                                     consumption_mu0 = 0, consumption_v0 = 5, consumption_alpha0 = 1, 
                                     consumption_beta0 = 1, suff_stat_consumption = TRUE,
                                     consumption_event_alpha0 = 1, consumption_event_beta0 = 1)
  
  save(TWI_pp, file = 'TWI_pp.Rdata')

  TWI_pp_high_consumer <-  unc_analysis_assessment(niter_ale = 10000, niter_epi = 10000, threshold = 1, percentile_ale = 95,
                                     data_concentration = data_assessment$log_concentration_ss_data, 
                                     data_consumption = data_assessment$log_consumption_ss_data, gen_data_EKE = gen_eke,
                                     consumers_info_sample_size = data_assessment$consumers_info_sample_size,
                                     concentration_mu0 = 3.5, concentration_v0 = 5, concentration_alpha0 = 1, 
                                     concentration_beta0 = 1, suff_stat_concentration = TRUE,
                                     consumption_mu0 = 0, consumption_v0 = 5, consumption_alpha0 = 1, 
                                     consumption_beta0 = 1, suff_stat_consumption = TRUE,
                                     consumption_event_alpha0 = 1, consumption_event_beta0 = 1)
  
  save(TWI_pp_high_consumer, file = 'TWI_pp_high_consumer.Rdata')
  
  
##########################################################################################
  ## Bounded probability 
{
  ### initial_mu0 = c(concentration_mu0, consumption_mu0)

  obj_func <- function(parameters, niter_ale = 1000, niter_epi = 1000,
                       threshold = 0.5, percentile_ale = 0,
                       data_concentration = data_assessment$log_concentration_ss_data, 
                       data_consumption = data_assessment$log_consumption_ss_data,
                       gen_data_EKE = gen_eke, consumers_info_sample_size = data_assessment$consumers_info_sample_size, 
                       concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, suff_stat_concentration = TRUE,
                       consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, suff_stat_consumption = TRUE,
                       consumption_event_alpha0 = 1, consumption_event_beta0 = 1){
    
    concentration_mu0 <- parameters[1] 
    consumption_mu0 <- parameters[2] 
    
    out <- unc_analysis_assessment(niter_ale = niter_ale, niter_epi= niter_epi, data_concentration = data_concentration, 
                                   threshold = threshold, percentile_ale = percentile_ale,
                                   data_consumption = data_consumption,  gen_data_EKE = gen_data_EKE,
                                   consumers_info_sample_size = consumers_info_sample_size, 
                                   concentration_mu0 = concentration_mu0, concentration_v0 = concentration_v0, 
                                   concentration_alpha0 = concentration_alpha0, concentration_beta0 = concentration_beta0,
                                   suff_stat_concentration = suff_stat_concentration,
                                   consumption_mu0 = consumption_mu0, consumption_v0 =  consumption_v0, 
                                   consumption_alpha0 = consumption_alpha0, consumption_beta0 = consumption_beta0, 
                                   suff_stat_consumption = suff_stat_consumption,
                                   consumption_event_alpha0 = consumption_event_alpha0, consumption_event_beta0 = consumption_event_beta0)
    
    output <- out$expected_prob_exceed
   
    return(output)
  }
 
  bound_prob  <- function(obj_func, maximize = 'TRUE', lower = c(1, -5), upper = c(6, 1),
                   niter_ale = 1000, niter_epi = 1000, threshold = 1, percentile_ale = 0,
                   data_concentration = data_assessment$log_concentration_ss_data,
                   data_consumption = data_assessment$log_consumption_ss_data,
                   gen_data_EKE = gen_eke, 
                   consumers_info_sample_size = data_assessment$consumers_info_sample_size,
                   concentration_mu0 = 5.75,
                   concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                   suff_stat_concentration = TRUE,
                   consumption_mu0 = 0.75,
                   consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                   suff_stat_consumption = TRUE,
                   consumption_event_alpha0 = 1, consumption_event_beta0 = 1){
    
    initial_mu0 = c(concentration_mu0, consumption_mu0)
    
    opt_value <- nmkb(par = initial_mu0, fn = obj_func, lower = lower, upper = upper,
                             control = list(maximize =  maximize),
                             niter_ale = niter_ale, niter_epi = niter_epi, threshold = threshold, percentile_ale = percentile_ale,
                             data_concentration = data_concentration, data_consumption = data_consumption, 
                             gen_data_EKE = gen_data_EKE, 
                             consumers_info_sample_size = data_assessment$consumers_info_sample_size, 
                             concentration_v0 = concentration_v0, concentration_alpha0 = concentration_alpha0,
                             concentration_beta0 = concentration_beta0, suff_stat_concentration = suff_stat_concentration,
                             consumption_v0 = consumption_v0, consumption_alpha0 = consumption_alpha0, 
                             consumption_beta0 = consumption_beta0, suff_stat_consumption = suff_stat_consumption,
                             consumption_event_alpha0 = consumption_event_alpha0, 
                             consumption_event_beta0 = consumption_event_beta0)
    
    out <- unc_analysis_assessment(niter_ale = niter_ale, niter_epi= niter_epi, data_concentration = data_concentration, 
                                   threshold = threshold, percentile_ale = percentile_ale,
                                   data_consumption = data_consumption,  gen_data_EKE = gen_data_EKE,
                                   consumers_info_sample_size = consumers_info_sample_size, 
                                   concentration_mu0 = opt_value$par[1], concentration_v0 = concentration_v0, 
                                   concentration_alpha0 = concentration_alpha0, concentration_beta0 = concentration_beta0,
                                   suff_stat_concentration = suff_stat_concentration,
                                   consumption_mu0 =  opt_value$par[2], consumption_v0 =  consumption_v0, 
                                   consumption_alpha0 = consumption_alpha0, consumption_beta0 = consumption_beta0, 
                                   suff_stat_consumption = suff_stat_consumption,
                                   consumption_event_alpha0 = consumption_event_alpha0, consumption_event_beta0 = consumption_event_beta0)
    
    return(list(opt_value = opt_value, opt_prob_exceed = out$prob_exceed))
  }
  
  graph_bp <- function(lower_points, upper_points){
     
    n_values <- length(lower_points)
    l_points <- rep(0,n_values)
    u_points <- rep(0,n_values)
    
    for(i in 1:n_values){
      l_points[i] <- min(lower_points[i],upper_points[i]) 
      u_points[i] <- max(lower_points[i],upper_points[i])
    }
    
    # data wide format
    data_plot <- data.frame(l_points = sort(l_points), u_points = sort(u_points), cdf = c(1:length(l_points)/ length(l_points)))
    x1 <-  data_plot[1,1]
    x1_end <-   data_plot[1,2]
    
    x2 <-  tail(data_plot[,1], 1)
    x2_end <- tail(data_plot[,2], 1)
    
    # data long format
    data_plot <- gather(data_plot, bound, values, l_points, u_points, factor_key = TRUE)
    
     
  p <- data_plot %>% 
     ggplot(aes(x = values, y = cdf, group = bound, col = bound)) +
      geom_line() +
      scale_color_manual(labels = c('Upper', 'Lower'), values = c('red', 'blue')) +
      guides(color = guide_legend("Bounds")) +
      labs(
        title = "Uncertainty",
        x = "Frequency of exceeding TWI",
        y = "cdf")
  
  p = p + geom_segment(x = x1, y = 0, xend = x1_end, yend = 0, col = 'blue') 
  p + geom_segment(x = x2, y = 1, xend = x2_end, yend = 1, col = 'red') 
    
  }
  
}
  #################################################
  
  upper_bound <- bound_prob(obj_func = obj_func, maximize = TRUE, lower = c(1, -5), upper = c(6, 1),
                          niter_ale = 1000, niter_epi = 1000, threshold = 1, percentile_ale = 0,
                          data_concentration = data_assessment$log_concentration_ss_data,
                          data_consumption = data_assessment$log_consumption_ss_data,
                          gen_data_EKE = gen_eke, 
                          consumers_info_sample_size = data_assessment$consumers_info_sample_size,
                          concentration_mu0 = 5.75,
                          concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                          suff_stat_concentration = TRUE,
                          consumption_mu0 = 0.75,
                          consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                          suff_stat_consumption = TRUE,
                          consumption_event_alpha0 = 1, consumption_event_beta0 = 1)
  
  save(upper_bound, file = 'upper_bound.Rdata')
  
  lower_bound <- bound_prob(obj_func = obj_func, maximize = FALSE, lower = c(1, -5), upper = c(6, 1),
                           niter_ale = 1000, niter_epi = 1000, threshold = 1, percentile_ale = 0,
                           data_concentration = data_assessment$log_concentration_ss_data,
                           data_consumption = data_assessment$log_consumption_ss_data,
                           gen_data_EKE = gen_eke, 
                           consumers_info_sample_size = data_assessment$consumers_info_sample_size,
                           concentration_mu0 = 2.75,
                           concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                           suff_stat_concentration = TRUE,
                           consumption_mu0 = -2.5,
                           consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                           suff_stat_consumption = TRUE,
                           consumption_event_alpha0 = 1, consumption_event_beta0 = 1)
  
  save(lower_bound, file = 'lower_bound.Rdata')
  
  upper_bound_high_consumer <- bound_prob(obj_func = obj_func, maximize = TRUE, lower = c(1, -5), upper = c(6, 1),
                           niter_ale = 1000, niter_epi = 1000, threshold = 1, percentile_ale = 95,
                           data_concentration = data_assessment$log_concentration_ss_data,
                           data_consumption = data_assessment$log_consumption_ss_data,
                           gen_data_EKE = gen_eke, 
                           consumers_info_sample_size = data_assessment$consumers_info_sample_size,
                           concentration_mu0 = 5.75,
                           concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                           suff_stat_concentration = TRUE,
                           consumption_mu0 = 0.75,
                           consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                           suff_stat_consumption = TRUE,
                           consumption_event_alpha0 = 1, consumption_event_beta0 = 1)
  
  save(upper_bound_high_consumer, file = 'upper_bound_high_consumer.Rdata')
  
  lower_bound_high_consumer <- bound_prob(obj_func = obj_func, maximize = FALSE, lower = c(1, -5), upper = c(6, 1),
                           niter_ale = 1000, niter_epi = 1000, threshold = 1, percentile_ale = 95,
                           data_concentration = data_assessment$log_concentration_ss_data,
                           data_consumption = data_assessment$log_consumption_ss_data,
                           gen_data_EKE = gen_eke, 
                           consumers_info_sample_size = data_assessment$consumers_info_sample_size,
                           concentration_mu0 = 2.75,
                           concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                           suff_stat_concentration = TRUE,
                           consumption_mu0 = -2.5,
                           consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                           suff_stat_consumption = TRUE,
                           consumption_event_alpha0 = 1, consumption_event_beta0 = 1)
  
  save(lower_bound_high_consumer, file = 'lower_bound_high_consumer.Rdata')
  
#############################################  
  ## Visualizations
  ## Precise probability
  load('TWI_pp.Rdata')
  graph_plot_pp(assessment_output = TWI_pp)
  
  load('TWI_pp_high_consumer.Rdata')
  graph_plot_pp(assessment_output = TWI_pp_high_consumer)
  
  ############################
  ## Bounded probability
  load('lower_bound.Rdata')
  load('upper_bound.Rdata') 
  graph_bp(lower_points = lower_bound$opt_prob_exceed, upper_points = upper_bound$opt_prob_exceed)
  
  # high consumer
  load('lower_bound_high_consumer.Rdata')
  load('upper_bound_high_consumer.Rdata') 
  graph_bp(lower_points = lower_bound_high_consumer$opt_prob_exceed, upper_points = upper_bound_high_consumer$opt_prob_exceed)
  
