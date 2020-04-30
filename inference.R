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
 
  ## Normal Normal-Gamma conjugate model (Bayesian inference). 
  ## sufficient statistics (sample mean and sample variance) for the normal likelihood are known for log-aluminium concentration and log-consumption 
  
  ### Description
  ## This function implements the Normal Normal-Gamma conjugate model using Bayesian inference.
  
 update_normal_gamma <- function(suff_stat_data, mu0 = 0, v0 = 5, alpha0 = 1, beta0 = 1, sufficient_statistics = TRUE){
    
    if(sufficient_statistics == TRUE){
      n <- suff_stat_data[['sample_size']]
      m <- suff_stat_data[['sample_mean']]
      S2 <- suff_stat_data[['sample_sd']]^2  # sample variance
    }
    else{
     n <- length(suff_stat_data) 
     m <- mean(suff_stat_data)
     S2 <- var(suff_stat_data)
    }
    
    v_n <-  v0 + n
    alpha_n <- alpha0 + n/2
    
    mu_n <- (v0 * mu0 + n * m ) / (v0 + n)
    beta_n <- beta0 + 1/2*S2*(n - 1) + (n*v0/(v0 + n) * (m - mu0)^2)
      
    output <- list(param = list(prior = list(mu = mu0, v = v0, alpha = alpha0, beta = beta0),
                                             posterior = list(mu = mu_n, v = v_n , alpha = alpha_n, beta = beta_n)))
    
    return(output)
  }
 
 ### Usage
 ## update_normal_gamma(suff_stat_data, mu0, v0, alpha0, beta0, sufficient_statistics)
 
 ### Arguments
 ## suff_stat_data        a vector of sufficient statistics: sample_size, sample_mean and sample_sd. If sufficient_statistics = FALSE, then it is vector of observed data 
 ## mu0                   prior hyperparameter mu0 for the normal-gamma distribution
 ## v0                    prior hyperparameter v0 for the normal-gamma distribution
 ## alpha0                prior hyperparameter alpha0 for the normal-gamma distribution
 ## beta0                 prior hyperparameter beta0 for the normal-gamma distribution
 ## sufficient_statistics logical; if TRUE, sufficient statistics, suff_stat_data, are given as input data, 
 ##                       otherwise suff_stat_data is given as observed data. Default is TRUE
 
 ### Details
 ## The conjugate model is (x| mu,tau) ~ Normal(mu, tau); (mu, tau) ~ Normal-Gamma(mu0,v0,alpha0,beta0)
 ## See wikipedia for details https://en.wikipedia.org/wiki/Conjugate_prior
 
 ### Value
 ## param      A list with the prior and posterior hyperparameters. It constains the following components
 ## prior      A list with the prior hyperparameters mu, v, alpha, beta
 ## posterior  A list with the posterior hyperparameters mu, v, alpha, beta

 ############################################################################
 ############################################################################
 ############################################################################
 
 ## Bernoulli-Beta conjugate model (Bayesian inference). 

 ### Description
 ## This function implements the Bernoulli Beta conjugate model using Bayesian inference.
 
  update_bernoulli_beta <- function(suff_stat_data, alpha0 = 1, beta0 = 1){
    n_non_consumer <- suff_stat_data[['non_consumer_sample_size']]
    n_consumer <- suff_stat_data[['consumer_sample_size']]
    
    alpha_n <- alpha0 + n_consumer
    beta_n <- beta0 + n_non_consumer
    
    output <- list(param = list(prior = list(alpha = alpha0, beta = beta0),
                                            posterior = list(alpha = alpha_n, beta = beta_n)))
    
    return(output)
  }
  
  ### Usage
  ## update_bernoulli_beta(suff_stat_data, alpha0 = 1, beta0 = 1)
  
  ## Arguments
  ## alpha0             prior hyperparameter alpha0 for the beta distribution
  ## beta0              prior hyperparameter beta0 for the beta distribution
  ## suff_stat_data     a vector of sufficient statistics: non_consumer_sample_size and consumer_sample_size
  
  ### Details
  ## The conjugate model is (B | pi) ~ Bernoulli(pi); pi ~ Beta(alpha0,beta0)
  ## See wikipedia for details https://en.wikipedia.org/wiki/Conjugate_prior
  
  ### Value
  ## param      A list with the prior and posterior hyperparameters. It constains the following components
  ## prior      A list with the prior hyperparameters alpha, beta
  ## posterior  A list with the posterior hyperparameters alpha, beta
  

  ############################################################################
  ############################################################################
  ############################################################################
  
  ### Description
  ## This function generates samples from a log-normal distribution 
  ## given the posterior hyperparameters of the normal-gamma distribution
  
  generate_samples_normal_gamma <- function(niter_ale, post, percentile_ale){
  
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
  
  
  ###Usage
  ## generate_samples_normal_gamma(niter_ale, post, percentile_ale)
  
  
  ## Arguments
  ## niter_ale      number of generated samples  
  ## post           the output of update_normal_gamma function. Post is a list with the prior and posterior hyperparameters of the Normal-Gamma distribution. Prior and posterior are also list with hyperparameters mu, v, alpha and beta. 
  ## percentile_ale a value that indicates if the assessment is done on all population by 0 or on a high consumer child by 95. Default is 0
  
  
  ### Value
  ## gen_sample      A vector with the generated samples
   
  
  ############################################################################
  ############################################################################
  ############################################################################
  
  ## Precise Probability EKE 
  ## The fitted distribution is beta(1.39, 1.33)
  
  ### Description
  ## Takes elicited probabilities as inputs, and fits parametric distributions using least squares on 
  ## the cumulative distribution function. If separate judgements from multiple experts are specified,
  ## the function will fit one set of distributions per expert.
  
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
    output <-  list(best_fit = best_fit, param = param)
    return(output)
  }
  
  ### Usage 
  ## quantify_uncertainty_pp_change_eke(vals, probs)
  
    
  ### Arguments
  ## vals     A vector of elicited values for one expert, or a matrix of elicited values for multiple 
  ##          experts (one column per expert). Note that the an elicited judgement about X should be of the form P(X<= vals[i,j]) = probs[i,j]
  ## probs    A vector of elicited probabilies for one expert, or a matrix of elicited values for multiple 
  ##          experts (one column per expert). A single vector can be used if the probabilities are the same 
  ##          for each expert. For each expert, the smallest elicited probability must be less than 0.4, 
  ##          and the largest elicited probability must be greater than 0.6.
  
  
  ### Details
  ## The parametric distributions are: Normal, T-student, Gamma, Log-Normal, log-T student and Beta distributions
  
  ### Value
  ## A list with the following components 
  ## best_fit    the best fitting distribution for each expert, determined by the smallest sum of squared errors.
  ## param       a vector with the parameters of the best fitted distribution.
  
  ############################################################################
  ############################################################################
  ############################################################################
  
  ## Combination of uncertainty
  
  ### Description
  ## This function estimates the probability of exceeding the weekly intake
  
  combine_uncertainty <- function(gen_data_concentration, gen_data_consumption, 
                                  gen_data_EKE, threshold, niter_ale){
    
    gen_data_concentration <- matrix(unlist(gen_data_concentration), ncol = 7, nrow = niter_ale)
    gen_data_consumption <- matrix(unlist(gen_data_consumption), ncol = 7, nrow = niter_ale)
    
    weekly_intake <- rowSums(gen_data_EKE * (gen_data_consumption * 0.007) * gen_data_concentration)
    
    prob_exceed_wi <- mean(weekly_intake > threshold)
    
    return(prob_exceed_wi = prob_exceed_wi)
    
  }
  
  ### Usage 
  ## combine_uncertainty(gen_data_concentration, gen_data_consumption, gen_data_EKE, threshold, niter_ale)
  
  ## Arguments
  ## gen_data_concentration   generated concentration samples after Bayesian inference (output of function generate_samples_normal_gamma)
  ## gen_data_consumption     generated consumption samples after Bayesian inference (output of function generate_samples_normal_gamma)
  ## gen_data_EKE             a generated sample from the fitted distribution to the expert's elicited values
  ## threshold                safety threshold
  ## niter_ale                number of generated samples
  
  ### Value
  ## prob_exceed_wi           the probability of exceeding the weekly intake
  
  ############################################################################
  ############################################################################
  ############################################################################
  
  ## Final assessment model -- Precise Probability 
  
  ### Description
  ## This function does the aluminium exposure assessment. It estimates the expected value and 
  ## the highest posterior density of the probability of exceeding the threshold

  unc_analysis_assessment <- function(niter_ale = 1000, niter_epi = 1000, 
                                      threshold = 1, percentile_ale = 0,
                                      suff_stat_concentration, suff_stat_consumption, gen_data_EKE, consumers_info_sample_size,
                                      concentration_mu0 = 3.5, concentration_v0 = 5, concentration_alpha0 = 1, 
                                      concentration_beta0 = 1, sufficient_statistics_concentration = TRUE,
                                      consumption_mu0 = 0, consumption_v0 = 5, consumption_alpha0 = 1, 
                                      consumption_beta0 = 1, sufficient_statistics_consumption = TRUE,
                                      consumption_event_alpha0 = 1, consumption_event_beta0 = 1){
    
    nr_products <-  length(suff_stat_concentration)
    prob_consumption <- parameters_consumption <- parameters_concentration <- vector('list', nr_products)
    
    # Probability of a child i consumes chocolate product k
    param_consumption <-  lapply(consumers_info_sample_size, update_bernoulli_beta, alpha0 = consumption_event_alpha0, beta0 = consumption_event_beta0)
    
    for(k in 1:nr_products){
      prob_consumption[[k]] <-  rbeta(1,shape1 = param_consumption[[k]]$param$posterior$alpha, shape2 = param_consumption[[k]]$param$posterior$beta)
    }
    
    prob_exceed <- rep(0, niter_epi)
    
    post_concentration <-  lapply(suff_stat_concentration, update_normal_gamma, mu0 = concentration_mu0,
                                  v0 = concentration_v0, alpha0 = concentration_alpha0, 
                                  beta0 = concentration_beta0, sufficient_statistics = sufficient_statistics_concentration)
    
    post_consumption <-  lapply(suff_stat_consumption, update_normal_gamma, mu0 = consumption_mu0, 
                                v0 = consumption_v0, alpha0 = consumption_alpha0, 
                                beta0 = consumption_beta0, sufficient_statistics = sufficient_statistics_consumption)
    for(j in 1:nr_products){
      parameters_concentration[[j]] <-  post_concentration[[j]]$param
      
      parameters_consumption[[j]] <-  post_consumption[[j]]$param
    }
    
    for(i in 1:niter_epi){
      
      gen_data_concentration <-  lapply(post_concentration, generate_samples_normal_gamma, niter_ale = niter_ale, percentile_ale = 0)
      
      gen_data_consumption <-  lapply(post_consumption, generate_samples_normal_gamma, niter_ale = niter_ale, percentile_ale = percentile_ale)
      
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
  
  ### Usage
  ## unc_analysis_assessment(niter_ale, niter_epi, threshold, percentile_ale,
  ##                         suff_stat_concentration, suff_stat_consumption, gen_data_EKE, consumers_info_sample_size,
  ##                         concentration_mu0, concentration_v0 = 5, concentration_alpha0, 
  ##                         concentration_beta0, sufficient_statistics_concentration,
  ##                         consumption_mu0, consumption_v0 = 5, consumption_alpha0, 
  ##                         consumption_beta0, sufficient_statistics_consumption,
  ##                         consumption_event_alpha0, consumption_event_beta0)

  ### Arguments
  
  ## niter_ale                            number of generated samples
  ## niter_epi                            number of generated parameters from the posterior distrbutions 
  ##                                      (it indicates the number of repetitions the assessment will be done)
  ## threshold                            safety threshold
  ## percentile_ale                       a value that indicates if the assessment is done on all population by 0 or on a high consumer child by 95. Default is 0
  ## suff_stat_concentration              a vector of sufficient statistics: sample_size, sample_mean and sample_sd 
  ##                                      corresponding to concentration. If sufficient_statistics_concentration = FALSE, 
  ##                                      then it is vector of observed data
  ## suff_stat_consumption                a vector of sufficient statistics: sample_size, sample_mean and sample_sd 
  ##                                      corresponding to consumption. If sufficient_statistics_consumption = FALSE, 
  ##                                      then it is vector of observed data 
  ## gen_data_EKE                         a generated sample from the fitted distribution to the expert's elicited values
  ## consumers_info_sample_size           a vector with the sample size of non_consumer_sample_size and consumer_sample_size
  ## concentration_mu0                    prior hyperparameter mu0 for the normal-gamma distribution corresponding to concentration 
  ## concentration_v0                     prior hyperparameter v0 for the normal-gamma distribution corresponding to concentration 
  ## concentration_alpha0                 prior hyperparameter alpha0 for the normal-gamma distribution  corresponding to concentration 
  ## concentration_beta0                  prior hyperparameter beta0 for the normal-gamma distribution corresponding to concentration 
  ## sufficient_statistics_concentration  logical; if TRUE, sufficient statistics (sample_size, sample_mean, sample_variance) 
  ##                                      corresponding to concentration are given as input data, otherwise 
  ##                                      sufficient_statistics_concentration is given as observed data. Default is TRUE
  ## consumption_mu0                      prior hyperparameter mu0 for the normal-gamma distribution corresponding to consumption
  ## consumption_v0                       prior hyperparameter v0 for the normal-gamma distribution corresponding to consumption
  ## consumption_alpha0                   prior hyperparameter alpha0 for the normal-gamma distribution corresponding to consumption 
  ## consumption_beta0                    prior hyperparameter beta0 for the normal-gamma distribution corresponding to consumption 
  ## sufficient_statistics_consumption    logical; if TRUE, sufficient statistics (sample_size, sample_mean, sample_variance) 
  ##                                      corresponding to consumption are given as input data, otherwise 
  ##                                      sufficient_statistics_consumption is given as observed data. Default is TRUE
  ## consumption_event_alpha0             prior hyperparameter alpha0 for the beta distribution corresponding to consumption event
  ## consumption_event_beta0              prior hyperparameter beta0 for the beta distribution corresponding to consumption event
  
  ### value
  
  ## prob_consumption_event    the estimated probability of consumption events
  ## parameters_concentration  a list with the values of the prior and posterior parameters of concentration
  ## parameters_consumption    a list with the prior and posterior parameters of consumption
  ## prob_exceed               a vector with the estimated probabilities of exceeding the threshold (the lenght is niter_epi)
  ## expected_prob_exceed      the expected value of the probability of exceeding the threshold
  ## hdi_prob_exceed           the highest posterior density interval of the probability of exceeding the threshold
  
 ###################################################################################
  ## Figures  
  
  ### Description
  ## Plot a cumulative distribution function.
  
  graph_pp <- function(prob_exceed){
    
    prob_cdf <- c(1:length(prob_exceed))/(length(prob_exceed))
    data_plot <- data.frame(prob_exceed = sort(prob_exceed), prob_cdf = prob_cdf)
    
    data_plot %>% 
      ggplot(mapping = aes(y = prob_cdf,  x = prob_exceed)) +
      geom_line() +
      labs(
        title = "Uncertainty",
        x = "Frequency of exceeding TWI",
        y = "cdf")
  }
  
  ### Usage
  ## graph_plot_pp(prob_exceed)
  
  ### Arguments
  ## prob_exceed   a vector of probabilities of exceedance
  
}

#########################################################################################################################################
#########################################################################################################################################
## Precise probability

## Generate one eke sample 
  fit <- quantify_uncertainty_pp_change_eke(vals = data_assessment$change_cons$vals, probs = data_assessment$change_cons$probs/100)
  gen_eke <- rbeta(1,shape1 = fit$param[1], shape2 = fit$param[2])

## Final assessment
  
  ## all population
  TWI_pp <-  unc_analysis_assessment(niter_ale = 1000, niter_epi = 1000, threshold = 1, percentile_ale = 0,
                                     suff_stat_concentration = data_assessment$log_concentration_ss_data, 
                                     suff_stat_consumption = data_assessment$log_consumption_ss_data, gen_data_EKE = gen_eke,
                                     consumers_info_sample_size = data_assessment$consumers_info_sample_size,
                                     concentration_mu0 = 3.5, concentration_v0 = 5, concentration_alpha0 = 1, 
                                     concentration_beta0 = 1, sufficient_statistics_concentration = TRUE,
                                     consumption_mu0 = 0, consumption_v0 = 5, consumption_alpha0 = 1, 
                                     consumption_beta0 = 1, sufficient_statistics_consumption = TRUE,
                                     consumption_event_alpha0 = 1, consumption_event_beta0 = 1)
  
  save(TWI_pp, file = 'TWI_pp.Rdata')

  ## A high consumer
  TWI_pp_high_consumer <-  unc_analysis_assessment(niter_ale = 1000, niter_epi = 1000, threshold = 1, percentile_ale = 95,
                                                   suff_stat_concentration = data_assessment$log_concentration_ss_data, 
                                                   suff_stat_consumption = data_assessment$log_consumption_ss_data, gen_data_EKE = gen_eke,
                                                   consumers_info_sample_size = data_assessment$consumers_info_sample_size,
                                                   concentration_mu0 = 3.5, concentration_v0 = 5, concentration_alpha0 = 1, 
                                                   concentration_beta0 = 1, sufficient_statistics_concentration = TRUE,
                                                   consumption_mu0 = 0, consumption_v0 = 5, consumption_alpha0 = 1, 
                                                   consumption_beta0 = 1, sufficient_statistics_consumption = TRUE,
                                                   consumption_event_alpha0 = 1, consumption_event_beta0 = 1)
  
  save(TWI_pp_high_consumer, file = 'TWI_pp_high_consumer.Rdata')
  
  
#################################################################################################
#################################################################################################
#################################################################################################

  ## Bounded probability 
{
  ### This is the objective function (a function to be optimized) 
  ### The objective function is the unc_analysis_assessment function where 
  ### concentration_mu0 and consumption_mu0 are the parameters and the rest of the inputs arguments are fixed.   
  
  obj_func <- function(parameters, niter_ale = 1000, niter_epi = 1000,
                       threshold = 0.5, percentile_ale = 0,
                       suff_stat_concentration = data_assessment$log_concentration_ss_data, 
                       suff_stat_consumption = data_assessment$log_consumption_ss_data,
                       gen_data_EKE = gen_eke, consumers_info_sample_size = data_assessment$consumers_info_sample_size, 
                       concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, sufficient_statistics_concentration = TRUE,
                       consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, sufficient_statistics_consumption = TRUE,
                       consumption_event_alpha0 = 1, consumption_event_beta0 = 1){
    
    concentration_mu0 <- parameters[1] 
    consumption_mu0 <- parameters[2] 
    
    out <- unc_analysis_assessment(niter_ale = niter_ale, niter_epi= niter_epi, 
                                   threshold = threshold, percentile_ale = percentile_ale,
                                   suff_stat_concentration = suff_stat_concentration, 
                                   suff_stat_consumption = suff_stat_consumption,  gen_data_EKE = gen_data_EKE,
                                   consumers_info_sample_size = consumers_info_sample_size, 
                                   concentration_mu0 = concentration_mu0, concentration_v0 = concentration_v0, 
                                   concentration_alpha0 = concentration_alpha0, concentration_beta0 = concentration_beta0,
                                   sufficient_statistics_concentration = sufficient_statistics_concentration,
                                   consumption_mu0 = consumption_mu0, consumption_v0 =  consumption_v0, 
                                   consumption_alpha0 = consumption_alpha0, consumption_beta0 = consumption_beta0, 
                                   sufficient_statistics_consumption = sufficient_statistics_consumption,
                                   consumption_event_alpha0 = consumption_event_alpha0, consumption_event_beta0 = consumption_event_beta0)
    
    output <- out$expected_prob_exceed
   
    return(output)
  }
  ### Usage
  ## obj_func(parameters, niter_ale, niter_epi, threshold, percentile_ale,
  ##          suff_stat_concentration, suff_stat_consumption, gen_data_EKE, consumers_info_sample_size,
  ##          concentration_v0, concentration_alpha0, concentration_beta0, sufficient_statistics_concentration,
  ##          consumption_v0, consumption_alpha0, consumption_beta0, sufficient_statistics_consumption,
  ##          consumption_event_alpha0, consumption_event_beta0)
  
  ### Arguments
  
  ## parameters                           parameters of the objective function (concentration_mu0 and consumption_mu0) 
  ##                                      which are the hyperparameters mu0 of the normal-gamma distribution corresponding to 
  ##                                      concentration and consumption respectively
  ## niter_ale                            number of generated samples
  ## niter_epi                            number of generated parameters from the posterior distrbutions 
  ##                                      (it indicates the number of repetitions the assessment will be done)
  ## threshold                            safety threshold
  ## percentile_ale                       a value that indicates if the assessment is done on all population by 0 or on a high consumer child by 95. Default is 0
  ## suff_stat_concentration              a vector of sufficient statistics: sample_size, sample_mean and sample_sd 
  ##                                      corresponding to concentration. If sufficient_statistics_concentration = FALSE, 
  ##                                      then it is vector of observed data
  ## suff_stat_consumption                a vector of sufficient statistics: sample_size, sample_mean and sample_sd 
  ##                                      corresponding to consumption. If sufficient_statistics_consumption = FALSE, 
  ##                                      then it is vector of observed data 
  ## gen_data_EKE                         a generated sample from the fitted distribution to the expert's elicited values
  ## consumers_info_sample_size           a vector with the sample size of non_consumer_sample_size and consumer_sample_size
  ## concentration_v0                     prior hyperparameter v0 for the normal-gamma distribution corresponding to concentration 
  ## concentration_alpha0                 prior hyperparameter alpha0 for the normal-gamma distribution  corresponding to concentration 
  ## concentration_beta0                  prior hyperparameter beta0 for the normal-gamma distribution corresponding to concentration 
  ## sufficient_statistics_concentration  logical; if TRUE, sufficient statistics (sample_size, sample_mean, sample_variance) 
  ##                                      corresponding to concentration are given as input data, otherwise 
  ##                                      sufficient_statistics_concentration is given as observed data. Default is TRUE
  ## consumption_v0                       prior hyperparameter v0 for the normal-gamma distribution corresponding to consumption
  ## consumption_alpha0                   prior hyperparameter alpha0 for the normal-gamma distribution corresponding to consumption 
  ## consumption_beta0                    prior hyperparameter beta0 for the normal-gamma distribution corresponding to consumption 
  ## sufficient_statistics_consumption    logical; if TRUE, sufficient statistics (sample_size, sample_mean, sample_variance) 
  ##                                      corresponding to consumption are given as input data, otherwise 
  ##                                      sufficient_statistics_consumption is given as observed data. Default is TRUE
  ## consumption_event_alpha0             prior hyperparameter alpha0 for the beta distribution corresponding to consumption event
  ## consumption_event_beta0              prior hyperparameter beta0 for the beta distribution corresponding to consumption event
  
  
  ### Value
  ## expected_prob_exceed                 the expected value of the probability of exceeding the threshold 
  
  
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  
  ### Description
  ## This function uses the Nelder-Mead algorithm  for the optimization procedure. This algorithm 
  ## is implemented in the 'nmkb' function from the {dfoptim} package. 
  ##  
  
  bound_prob_exceed  <- function(obj_func, maximize = FALSE, lower_mu0 = c(1, -5), upper_mu0 = c(6, 1),
                          niter_ale = 1000, niter_epi = 1000, threshold = 1, percentile_ale = 0,
                          suff_stat_concentration = data_assessment$log_concentration_ss_data,
                          suff_stat_consumption = data_assessment$log_consumption_ss_data,
                          gen_data_EKE = gen_eke, 
                          consumers_info_sample_size = data_assessment$consumers_info_sample_size,
                          concentration_mu0 = 2.75,
                          concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                          sufficient_statistics_concentration = TRUE,
                          consumption_mu0 = -2.5,
                          consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                          sufficient_statistics_consumption = TRUE,
                          consumption_event_alpha0 = 1, consumption_event_beta0 = 1){
    
    initial_mu0 = c(concentration_mu0, consumption_mu0)
    
    opt_value <- nmkb(par = initial_mu0, fn = obj_func, lower = lower_mu0, upper = upper_mu0,
                      control = list(maximize =  maximize),
                      niter_ale = niter_ale, niter_epi = niter_epi, threshold = threshold, percentile_ale = percentile_ale,
                      suff_stat_concentration = suff_stat_concentration, suff_stat_consumption = suff_stat_consumption, 
                      gen_data_EKE = gen_data_EKE, 
                      consumers_info_sample_size = data_assessment$consumers_info_sample_size, 
                      concentration_v0 = concentration_v0, concentration_alpha0 = concentration_alpha0,
                      concentration_beta0 = concentration_beta0, sufficient_statistics_concentration = sufficient_statistics_concentration,
                      consumption_v0 = consumption_v0, consumption_alpha0 = consumption_alpha0, 
                      consumption_beta0 = consumption_beta0, sufficient_statistics_consumption = sufficient_statistics_consumption,
                      consumption_event_alpha0 = consumption_event_alpha0, 
                      consumption_event_beta0 = consumption_event_beta0)
    
    out <- unc_analysis_assessment(niter_ale = niter_ale, niter_epi= niter_epi,  
                                   threshold = threshold, percentile_ale = percentile_ale,
                                   suff_stat_concentration = suff_stat_concentration,
                                   suff_stat_consumption = suff_stat_consumption,  gen_data_EKE = gen_data_EKE,
                                   consumers_info_sample_size = consumers_info_sample_size, 
                                   concentration_mu0 = opt_value$par[1], concentration_v0 = concentration_v0, 
                                   concentration_alpha0 = concentration_alpha0, concentration_beta0 = concentration_beta0,
                                   sufficient_statistics_concentration = sufficient_statistics_concentration,
                                   consumption_mu0 =  opt_value$par[2], consumption_v0 =  consumption_v0, 
                                   consumption_alpha0 = consumption_alpha0, consumption_beta0 = consumption_beta0, 
                                   sufficient_statistics_consumption = sufficient_statistics_consumption,
                                   consumption_event_alpha0 = consumption_event_alpha0, consumption_event_beta0 = consumption_event_beta0)
    
    return(list(opt_value = opt_value, opt_prob_exceed = out$prob_exceed))
  }
  
  ### Usage
  ##  bound_prob_exceed(obj_func, maximize, lower_mu0, upper_mu0, niter_ale, niter_epi, threshold = 1, percentile_ale = 0,
  ##                    suff_stat_concentration, suff_stat_consumption, gen_data_EKE, consumers_info_sample_size,
  ##                    concentration_mu0, concentration_v0, concentration_alpha0, concentration_beta0, 
  ##                    sufficient_statistics_concentration,
  ##                    consumption_mu0, consumption_v0, consumption_alpha0, consumption_beta0, 
  ##                    sufficient_statistics_consumption,
  ##                    consumption_event_alpha0, consumption_event_beta0)
  
  
  ### Arguments
  
  ## obj_function                         the function to be optimized
  ## maximize                             a logical variable indicating whether the objective function should be maximized. Default is FALSE.
  ## lower_mu0                            lower bounds on the parameters mu0
  ## upper_mu0                            upper bounds on the parameters mu0
  ## niter_ale                            number of generated samples
  ## niter_epi                            number of generated parameters from the posterior distrbutions 
  ##                                      (it indicates the number of repetitions the assessment will be done)
  ## threshold                            safety threshold
  ## percentile_ale                       a value that indicates if the assessment is done on all population by 0 or on a high consumer child by 95. Default is 0
  ## suff_stat_concentration              a vector of sufficient statistics: sample_size, sample_mean and sample_sd 
  ##                                      corresponding to concentration. If sufficient_statistics_concentration = FALSE, 
  ##                                      then it is vector of observed data
  ## suff_stat_consumption                a vector of sufficient statistics: sample_size, sample_mean and sample_sd 
  ##                                      corresponding to consumption. If sufficient_statistics_consumption = FALSE, 
  ##                                      then it is vector of observed data 
  ## gen_data_EKE                         a generated sample from the fitted distribution to the expert's elicited values
  ## consumers_info_sample_size           a vector with the sample size of non_consumer_sample_size and consumer_sample_size
  ## concentration_mu0                    prior hyperparameter mu0 for the normal-gamma distribution corresponding to concentration 
  ## concentration_v0                     prior hyperparameter v0 for the normal-gamma distribution corresponding to concentration 
  ## concentration_alpha0                 prior hyperparameter alpha0 for the normal-gamma distribution  corresponding to concentration 
  ## concentration_beta0                  prior hyperparameter beta0 for the normal-gamma distribution corresponding to concentration 
  ## sufficient_statistics_concentration  logical; if TRUE, sufficient statistics (sample_size, sample_mean, sample_variance) 
  ##                                      corresponding to concentration are given as input data, otherwise 
  ##                                      sufficient_statistics_concentration is given as observed data. Default is TRUE
  ## consumption_mu0                      prior hyperparameter mu0 for the normal-gamma distribution corresponding to consumption
  ## consumption_v0                       prior hyperparameter v0 for the normal-gamma distribution corresponding to consumption
  ## consumption_alpha0                   prior hyperparameter alpha0 for the normal-gamma distribution corresponding to consumption 
  ## consumption_beta0                    prior hyperparameter beta0 for the normal-gamma distribution corresponding to consumption 
  ## sufficient_statistics_consumption    logical; if TRUE, sufficient statistics (sample_size, sample_mean, sample_variance) 
  ##                                      corresponding to consumption are given as input data, otherwise 
  ##                                      sufficient_statistics_consumption is given as observed data. Default is TRUE
  ## consumption_event_alpha0             prior hyperparameter alpha0 for the beta distribution corresponding to consumption event
  ## consumption_event_beta0              prior hyperparameter beta0 for the beta distribution corresponding to consumption event
  
  
  ## value (output)
  
  ## opt_value is a list with the following components:
  ## par                                  Best estimate of the parameter vector found by the algorithm
  ## value                                The value of the objective function at termination
  ## feval	                              The number of times the objective fn was evaluated
  ## restarts	                            The number of times the algorithm had to be restarted when it stagnated.
  ## convergence	                        An integer code indicating type of convergence, 0 indicates successful convergence. Positive integer codes indicate failure to converge.
  ## message	                            Text message indicating the type of convergence or failure
  ## opt_prob_exceed                      A vector with the estimated probabilities of exceeding the threshold corresponding to the solution
 
  
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  
  ### Description
  ## Plot a p-box (a p-box is represented by two cdf sequences)
   
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
  
  ### Usage
  ## graph_bp(lower_points, upper_points) 
  
  ### Arguments
  ## lower_points               A vector with the lower cdf sequence
  ## upper_points               A vector with the upper cdf sequence
  
  
  ###################################################################################################
  ###################################################################################################
  ###################################################################################################
  
}
  ####################################################################################################
  
  ## maximum - all population 
  upper_bound_prob_exceeding_TWI <- bound_prob_exceed(obj_func = obj_func, maximize = TRUE, lower_mu0 = c(1, -5), upper_mu0 = c(6, 1),
                          niter_ale = 1000, niter_epi = 1000, threshold = 1, percentile_ale = 0,
                          suff_stat_concentration = data_assessment$log_concentration_ss_data,
                          suff_stat_consumption = data_assessment$log_consumption_ss_data,
                          gen_data_EKE = gen_eke, 
                          consumers_info_sample_size = data_assessment$consumers_info_sample_size,
                          concentration_mu0 = 5.75,
                          concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                          sufficient_statistics_concentration = TRUE,
                          consumption_mu0 = 0.75,
                          consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                          sufficient_statistics_consumption = TRUE,
                          consumption_event_alpha0 = 1, consumption_event_beta0 = 1)
  
  save(upper_bound_prob_exceeding_TWI, file = 'upper_bound_prob_exceeding_TWI.Rdata')
  
  
  ## minimum - all population 
  lower_bound_prob_exceeding_TWI <- bound_prob_exceed(obj_func = obj_func, maximize = FALSE, lower_mu0 = c(1, -5), upper_mu0 = c(6, 1),
                           niter_ale = 1000, niter_epi = 1000, threshold = 1, percentile_ale = 0,
                           suff_stat_concentration = data_assessment$log_concentration_ss_data,
                           suff_stat_consumption = data_assessment$log_consumption_ss_data,
                           gen_data_EKE = gen_eke, 
                           consumers_info_sample_size = data_assessment$consumers_info_sample_size,
                           concentration_mu0 = 2.75,
                           concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                           sufficient_statistics_concentration = TRUE,
                           consumption_mu0 = -2.5,
                           consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                           sufficient_statistics_consumption = TRUE,
                           consumption_event_alpha0 = 1, consumption_event_beta0 = 1)
  
  save(lower_bound_prob_exceeding_TWI, file = 'lower_bound_prob_exceeding_TWI.Rdata')
  
  # maximum - high consumer
  
  upper_bound_prob_exceeding_TWI_high_consumer <- bound_prob_exceed(obj_func = obj_func, maximize = TRUE, lower_mu0 = c(1, -5), upper_mu0 = c(6, 1),
                           niter_ale = 1000, niter_epi = 1000, threshold = 1, percentile_ale = 95,
                           suff_stat_concentration = data_assessment$log_concentration_ss_data,
                           suff_stat_consumption = data_assessment$log_consumption_ss_data,
                           gen_data_EKE = gen_eke, 
                           consumers_info_sample_size = data_assessment$consumers_info_sample_size,
                           concentration_mu0 = 5.75,
                           concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                           sufficient_statistics_concentration = TRUE,
                           consumption_mu0 = 0.75,
                           consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                           sufficient_statistics_consumption = TRUE,
                           consumption_event_alpha0 = 1, consumption_event_beta0 = 1)
  
  save(upper_bound_prob_exceeding_TWI_high_consumer, file = 'upper_bound_prob_exceeding_TWI_high_consumer.Rdata')
  
  # minimum - high consumer
  
  lower_bound_prob_exceeding_TWI_high_consumer <- bound_prob_exceed(obj_func = obj_func, maximize = FALSE, lower_mu0 = c(1, -5), upper_mu0 = c(6, 1),
                           niter_ale = 1000, niter_epi = 1000, threshold = 1, percentile_ale = 95,
                           suff_stat_concentration = data_assessment$log_concentration_ss_data,
                           suff_stat_consumption = data_assessment$log_consumption_ss_data,
                           gen_data_EKE = gen_eke, 
                           consumers_info_sample_size = data_assessment$consumers_info_sample_size,
                           concentration_mu0 = 2.75,
                           concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                           sufficient_statistics_concentration = TRUE,
                           consumption_mu0 = -2.5,
                           consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                           sufficient_statistics_consumption = TRUE,
                           consumption_event_alpha0 = 1, consumption_event_beta0 = 1)
  
  save(lower_bound_prob_exceeding_TWI_high_consumer, file = 'lower_bound_prob_exceeding_TWI_high_consumer.Rdata')
  
  #########################################################################################################################
  #########################################################################################################################
  #########################################################################################################################
  
  ## Visualizations
  ## Precise probability
  
  ## all population
  load('TWI_pp.Rdata')
  graph_pp(prob_exceed = TWI_pp$prob_exceed)
  
  # high consumer
  load('TWI_pp_high_consumer.Rdata')
  graph_pp(prob_exceed = TWI_pp_high_consumer$prob_exceed)
  
  ############################
  ## Bounded probability
  
  ## all population
  load('lower_bound.Rdata')
  load('upper_bound.Rdata') 
  graph_bp(lower_points = lower_bound$opt_prob_exceed, upper_points = upper_bound$opt_prob_exceed)
  
  # high consumer
  load('lower_bound_high_consumer.Rdata')
  load('upper_bound_high_consumer.Rdata') 
  graph_bp(lower_points = lower_bound_high_consumer$opt_prob_exceed, upper_points = upper_bound_high_consumer$opt_prob_exceed)
  
  #########################################################################################################################
  #########################################################################################################################
  #########################################################################################################################
  ## Now with 3 parameters (concentration_mu0, consumption_mu0, EKE_alpha) 
  
  ## ## Bounded probability 
  { 
    
    unc_analysis_assessment_3_param <- function(niter_ale = 1000, niter_epi = 1000, 
                                        threshold = 1, percentile_ale = 0,
                                        suff_stat_concentration, suff_stat_consumption, consumers_info_sample_size,
                                        concentration_mu0 = 3.5, concentration_v0 = 5, concentration_alpha0 = 1, 
                                        concentration_beta0 = 1, sufficient_statistics_concentration = TRUE,
                                        consumption_mu0 = 0, consumption_v0 = 5, consumption_alpha0 = 1, 
                                        consumption_beta0 = 1, sufficient_statistics_consumption = TRUE,
                                        consumption_event_alpha0 = 1, consumption_event_beta0 = 1,
                                        EKE_alpha = 1.391, EKE_beta = 1.177){
      
      nr_products <-  length(suff_stat_concentration)
      prob_consumption <- parameters_consumption <- parameters_concentration <- vector('list', nr_products)
      
      # Probability of a child i consumes chocolate product k
      param_consumption <-  lapply(consumers_info_sample_size, update_bernoulli_beta, alpha0 = consumption_event_alpha0, beta0 = consumption_event_beta0)
      
      for(k in 1:nr_products){
        prob_consumption[[k]] <-  rbeta(1,shape1 = param_consumption[[k]]$param$posterior$alpha, shape2 = param_consumption[[k]]$param$posterior$beta)
      }
      
      prob_exceed <- rep(0, niter_epi)
      
      post_concentration <-  lapply(suff_stat_concentration, update_normal_gamma, mu0 = concentration_mu0,
                                    v0 = concentration_v0, alpha0 = concentration_alpha0, 
                                    beta0 = concentration_beta0, sufficient_statistics = sufficient_statistics_concentration)
      
      post_consumption <-  lapply(suff_stat_consumption, update_normal_gamma, mu0 = consumption_mu0, 
                                  v0 = consumption_v0, alpha0 = consumption_alpha0, 
                                  beta0 = consumption_beta0, sufficient_statistics = sufficient_statistics_consumption)
      for(j in 1:nr_products){
        parameters_concentration[[j]] <-  post_concentration[[j]]$param
        
        parameters_consumption[[j]] <-  post_consumption[[j]]$param
      }
      
      gen_data_EKE <- rbeta(1, shape1 = EKE_alpha, shape2 = EKE_beta)
      
      for(i in 1:niter_epi){
        
        gen_data_concentration <-  lapply(post_concentration, generate_samples_normal_gamma, niter_ale = niter_ale, percentile_ale = 0)
        
        gen_data_consumption <-  lapply(post_consumption, generate_samples_normal_gamma, niter_ale = niter_ale, percentile_ale = percentile_ale)
        
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
    
    
    #########################################################################################################################
    #########################################################################################################################
    #########################################################################################################################
    
    ### Description
    ### This is the objective function (a function to be optimized) 
    ### The objective function is the unc_analysis_assessment_3_param function where 
    ### concentration_mu0, consumption_mu0 and EKE_alpha are the parameters and the rest of the inputs arguments are fixed.   
    
    obj_func_3_param <- function(parameters, niter_ale = 1000, niter_epi = 1000,
                         threshold = 0.5, percentile_ale = 0,
                         suff_stat_concentration = data_assessment$log_concentration_ss_data, 
                         suff_stat_consumption = data_assessment$log_consumption_ss_data,
                         consumers_info_sample_size = data_assessment$consumers_info_sample_size, 
                         concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, sufficient_statistics_concentration = TRUE,
                         consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, sufficient_statistics_consumption = TRUE,
                         consumption_event_alpha0 = 1, consumption_event_beta0 = 1,
                         EKE_beta = 1.177){
      
      concentration_mu0 <- parameters[1] 
      consumption_mu0 <- parameters[2] 
      EKE_alpha <- parameters[3] 
      
      
      out <- unc_analysis_assessment_3_param(niter_ale = niter_ale, niter_epi= niter_epi, 
                                     threshold = threshold, percentile_ale = percentile_ale,
                                     suff_stat_concentration = suff_stat_concentration, 
                                     suff_stat_consumption = suff_stat_consumption,
                                     consumers_info_sample_size = consumers_info_sample_size, 
                                     concentration_mu0 = concentration_mu0, concentration_v0 = concentration_v0, 
                                     concentration_alpha0 = concentration_alpha0, concentration_beta0 = concentration_beta0,
                                     sufficient_statistics_concentration = sufficient_statistics_concentration,
                                     consumption_mu0 = consumption_mu0, consumption_v0 =  consumption_v0, 
                                     consumption_alpha0 = consumption_alpha0, consumption_beta0 = consumption_beta0, 
                                     sufficient_statistics_consumption = sufficient_statistics_consumption,
                                     consumption_event_alpha0 = consumption_event_alpha0, consumption_event_beta0 = consumption_event_beta0,
                                     EKE_alpha = EKE_alpha, EKE_beta = EKE_beta)
      
      output <- out$expected_prob_exceed
      
      return(output)
    }
    
    ### Usage
    ##  obj_func_3_param(parameters, niter_ale, niter_epi, threshold, percentile_ale,
    ##                    suff_stat_concentration, suff_stat_consumption, consumers_info_sample_size,
    ##                    concentration_v0, concentration_alpha0, concentration_beta0, sufficient_statistics_concentration,
    ##                    consumption_v0, consumption_alpha0, consumption_beta0, sufficient_statistics_consumption,
    ##                    consumption_event_alpha0, consumption_event_beta0, EKE_beta)
    
    
    ### Arguments
    
    ## parameters                           parameters of the objective function (concentration_mu0 and consumption_mu0) 
    ##                                      which are the hyperparameters mu0 of the normal-gamma distribution corresponding to 
    ##                                      concentration and consumption respectively
    ## niter_ale                            number of generated samples
    ## niter_epi                            number of generated parameters from the posterior distrbutions 
    ##                                      (it indicates the number of repetitions the assessment will be done)
    ## threshold                            safety threshold
    ## percentile_ale                       a value that indicates if the assessment is done on all population by 0 or on a high consumer child by 95. Default is 0
    ## suff_stat_concentration              a vector of sufficient statistics: sample_size, sample_mean and sample_sd 
    ##                                      corresponding to concentration. If sufficient_statistics_concentration = FALSE, 
    ##                                      then it is vector of observed data
    ## suff_stat_consumption                a vector of sufficient statistics: sample_size, sample_mean and sample_sd 
    ##                                      corresponding to consumption. If sufficient_statistics_consumption = FALSE, 
    ##                                      then it is vector of observed data 
    ## consumers_info_sample_size           a vector with the sample size of non_consumer_sample_size and consumer_sample_size
    ## concentration_v0                     prior hyperparameter v0 for the normal-gamma distribution corresponding to concentration 
    ## concentration_alpha0                 prior hyperparameter alpha0 for the normal-gamma distribution  corresponding to concentration 
    ## concentration_beta0                  prior hyperparameter beta0 for the normal-gamma distribution corresponding to concentration 
    ## sufficient_statistics_concentration  logical; if TRUE, sufficient statistics (sample_size, sample_mean, sample_variance) 
    ##                                      corresponding to concentration are given as input data, otherwise 
    ##                                      sufficient_statistics_concentration is given as observed data. Default is TRUE
    ## consumption_v0                       prior hyperparameter v0 for the normal-gamma distribution corresponding to consumption
    ## consumption_alpha0                   prior hyperparameter alpha0 for the normal-gamma distribution corresponding to consumption 
    ## consumption_beta0                    prior hyperparameter beta0 for the normal-gamma distribution corresponding to consumption 
    ## sufficient_statistics_consumption    logical; if TRUE, sufficient statistics (sample_size, sample_mean, sample_variance) 
    ##                                      corresponding to consumption are given as input data, otherwise 
    ##                                      sufficient_statistics_consumption is given as observed data. Default is TRUE
    ## consumption_event_alpha0             prior hyperparameter alpha0 for the beta distribution corresponding to consumption event
    ## consumption_event_beta0              prior hyperparameter beta0 for the beta distribution corresponding to consumption event
    ## EKE_alpha0                           prior hyperparameter alpha0 for the beta distribution corresponding to EKE
    
    
    ### Value
    ## expected_prob_exceed                 the expected value of the probability of exceeding the threshold 
    
    
    ###################################################################################################
    ###################################################################################################
    ###################################################################################################
    
    ### Description
    ## This function uses the Nelder-Mead algorithm  for the optimization procedure. This algorithm 
    ## is implemented in the 'nmkb' function from the {dfoptim} package. 
    ##  
    
    bound_prob_exceed_3_param  <- function(obj_func_3_param, maximize = FALSE, lower_parameters  = c(1, -5, 1), upper_parameters  = c(6, 1, 2),
                                   niter_ale = 1000, niter_epi = 1000, threshold = 1, percentile_ale = 0,
                                   suff_stat_concentration = data_assessment$log_concentration_ss_data,
                                   suff_stat_consumption = data_assessment$log_consumption_ss_data,
                                   consumers_info_sample_size = data_assessment$consumers_info_sample_size,
                                   concentration_mu0 = 2.75,
                                   concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                                   sufficient_statistics_concentration = TRUE,
                                   consumption_mu0 = -2.5,
                                   consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                                   sufficient_statistics_consumption = TRUE,
                                   consumption_event_alpha0 = 1, consumption_event_beta0 = 1,
                                   EKE_alpha = 1.391, EKE_beta = 1.177){
      
      initial_parameters <- c(concentration_mu0, consumption_mu0, EKE_alpha)
      
      opt_value <- nmkb(par = initial_parameters, fn = obj_func_3_param, lower = lower_parameters, upper = upper_parameters,
                        control = list(maximize =  maximize),
                        niter_ale = niter_ale, niter_epi = niter_epi, threshold = threshold, percentile_ale = percentile_ale,
                        suff_stat_concentration = suff_stat_concentration, suff_stat_consumption = suff_stat_consumption, 
                        consumers_info_sample_size = consumers_info_sample_size, 
                        concentration_v0 = concentration_v0, concentration_alpha0 = concentration_alpha0,
                        concentration_beta0 = concentration_beta0, sufficient_statistics_concentration = sufficient_statistics_concentration,
                        consumption_v0 = consumption_v0, consumption_alpha0 = consumption_alpha0, 
                        consumption_beta0 = consumption_beta0, sufficient_statistics_consumption = sufficient_statistics_consumption,
                        consumption_event_alpha0 = consumption_event_alpha0, 
                        consumption_event_beta0 = consumption_event_beta0, EKE_beta = EKE_beta)
      
      out_prob <- unc_analysis_assessment_3_param(niter_ale = niter_ale, niter_epi= niter_epi,  
                                     threshold = threshold, percentile_ale = percentile_ale,
                                     suff_stat_concentration = suff_stat_concentration,
                                     suff_stat_consumption = suff_stat_consumption,
                                     consumers_info_sample_size = consumers_info_sample_size, 
                                     concentration_mu0 = opt_value$par[1], concentration_v0 = concentration_v0, 
                                     concentration_alpha0 = concentration_alpha0, concentration_beta0 = concentration_beta0,
                                     sufficient_statistics_concentration = sufficient_statistics_concentration,
                                     consumption_mu0 =  opt_value$par[2], consumption_v0 =  consumption_v0, 
                                     consumption_alpha0 = consumption_alpha0, consumption_beta0 = consumption_beta0, 
                                     sufficient_statistics_consumption = sufficient_statistics_consumption,
                                     consumption_event_alpha0 = consumption_event_alpha0, consumption_event_beta0 = consumption_event_beta0,
                                     EKE_alpha = opt_value$par[3], EKE_beta = EKE_beta)
      
      return(list(opt_value = opt_value, opt_prob_exceed = out_prob$prob_exceed))
    }
    
    ###
    
    min = bound_prob_exceed_3_param (obj_func_3_param = obj_func_3_param, maximize = FALSE, lower_parameters  = c(1, -5, 1), upper_parameters  = c(6, 1, 2),
                                     niter_ale = 1000, niter_epi = 1000, threshold = 1, percentile_ale = 0,
                                     suff_stat_concentration = data_assessment$log_concentration_ss_data,
                                     suff_stat_consumption = data_assessment$log_consumption_ss_data,
                                     consumers_info_sample_size = data_assessment$consumers_info_sample_size,
                                     concentration_mu0 = 2.75,
                                     concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                                     sufficient_statistics_concentration = TRUE,
                                     consumption_mu0 = -2.5,
                                     consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                                     sufficient_statistics_consumption = TRUE,
                                     consumption_event_alpha0 = 1, consumption_event_beta0 = 1,
                                     EKE_alpha = 1.391, EKE_beta = 1.177)
    
    max = bound_prob_exceed_3_param(obj_func_3_param = obj_func_3_param, maximize = TRUE, lower_parameters  = c(1, -5, 1), upper_parameters  = c(6, 1, 2),
                                     niter_ale = 1000, niter_epi = 1000, threshold = 1, percentile_ale = 0,
                                     suff_stat_concentration = data_assessment$log_concentration_ss_data,
                                     suff_stat_consumption = data_assessment$log_consumption_ss_data,
                                     consumers_info_sample_size = data_assessment$consumers_info_sample_size,
                                     concentration_mu0 = 5.75,
                                     concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                                     sufficient_statistics_concentration = TRUE,
                                     consumption_mu0 = 0.75,
                                     consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                                     sufficient_statistics_consumption = TRUE,
                                     consumption_event_alpha0 = 1, consumption_event_beta0 = 1,
                                     EKE_alpha = 1.391, EKE_beta = 1.177)
    
    ### Usage
    ##  bound_prob_exceed_3_param(obj_func, maximize, lower_mu0, upper_mu0, niter_ale, niter_epi, threshold = 1, percentile_ale = 0,
    ##                    suff_stat_concentration, suff_stat_consumption, gen_data_EKE, consumers_info_sample_size,
    ##                    concentration_mu0, concentration_v0, concentration_alpha0, concentration_beta0, 
    ##                    sufficient_statistics_concentration,
    ##                    consumption_mu0, consumption_v0, consumption_alpha0, consumption_beta0, 
    ##                    sufficient_statistics_consumption,
    ##                    consumption_event_alpha0, consumption_event_beta0)
    
    
    ### Arguments
    
    ## obj_function                         the function to be optimized
    ## maximize                             a logical variable indicating whether the objective function should be maximized. Default is FALSE.
    ## lower_parameters                     lower bounds on the parameters 
    ## upper_parameters                     upper bounds on the parameters 
    ## niter_ale                            number of generated samples
    ## niter_epi                            number of generated parameters from the posterior distrbutions 
    ##                                      (it indicates the number of repetitions the assessment will be done)
    ## threshold                            safety threshold
    ## percentile_ale                       a value that indicates if the assessment is done on all population by 0 or on a high consumer child by 95. Default is 0
    ## suff_stat_concentration              a vector of sufficient statistics: sample_size, sample_mean and sample_sd 
    ##                                      corresponding to concentration. If sufficient_statistics_concentration = FALSE, 
    ##                                      then it is vector of observed data
    ## suff_stat_consumption                a vector of sufficient statistics: sample_size, sample_mean and sample_sd 
    ##                                      corresponding to consumption. If sufficient_statistics_consumption = FALSE, 
    ##                                      then it is vector of observed data 
    ## consumers_info_sample_size           a vector with the sample size of non_consumer_sample_size and consumer_sample_size
    ## concentration_mu0                    prior hyperparameter mu0 for the normal-gamma distribution corresponding to concentration 
    ## concentration_v0                     prior hyperparameter v0 for the normal-gamma distribution corresponding to concentration 
    ## concentration_alpha0                 prior hyperparameter alpha0 for the normal-gamma distribution  corresponding to concentration 
    ## concentration_beta0                  prior hyperparameter beta0 for the normal-gamma distribution corresponding to concentration 
    ## sufficient_statistics_concentration  logical; if TRUE, sufficient statistics (sample_size, sample_mean, sample_variance) 
    ##                                      corresponding to concentration are given as input data, otherwise 
    ##                                      sufficient_statistics_concentration is given as observed data. Default is TRUE
    ## consumption_mu0                      prior hyperparameter mu0 for the normal-gamma distribution corresponding to consumption
    ## consumption_v0                       prior hyperparameter v0 for the normal-gamma distribution corresponding to consumption
    ## consumption_alpha0                   prior hyperparameter alpha0 for the normal-gamma distribution corresponding to consumption 
    ## consumption_beta0                    prior hyperparameter beta0 for the normal-gamma distribution corresponding to consumption 
    ## sufficient_statistics_consumption    logical; if TRUE, sufficient statistics (sample_size, sample_mean, sample_variance) 
    ##                                      corresponding to consumption are given as input data, otherwise 
    ##                                      sufficient_statistics_consumption is given as observed data. Default is TRUE
    ## consumption_event_alpha0             prior hyperparameter alpha0 for the beta distribution corresponding to consumption event
    ## consumption_event_beta0              prior hyperparameter beta0 for the beta distribution corresponding to consumption event
    ## EKE_alpha0                           prior hyperparameter alpha0 for the beta distribution corresponding to EKE
    
    
    ## value (output)
    
    ## opt_value is a list with the following components:
    ## par                                  Best estimate of the parameter vector found by the algorithm
    ## value                                The value of the objective function at termination
    ## feval	                              The number of times the objective fn was evaluated
    ## restarts	                            The number of times the algorithm had to be restarted when it stagnated.
    ## convergence	                        An integer code indicating type of convergence, 0 indicates successful convergence. Positive integer codes indicate failure to converge.
    ## message	                            Text message indicating the type of convergence or failure
    ## opt_prob_exceed                      A vector with the estimated probabilities of exceeding the threshold corresponding to the solution
    
    
    ###################################################################################################
    ###################################################################################################
    ###################################################################################################
    
    
  }
  
  
