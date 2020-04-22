##########
# Code for paper comparing precise to bounded probability 
# Ivette Raices Cruz, Matthias C M Troffaes, Ullrika Sahlin
# 17 april 2020

########################
library('SHELF')
library('HDInterval')
library('dfoptim')

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
      
    output = list(data = data, param = list(prior = list(mu = mu0, v = v0, alpha = alpha0, beta = beta0),
                                            posterior = list(mu = mu_n, v = v_n , alpha = alpha_n, beta = beta_n)))
      
    return(output)
  }
  
  ## Bernoulli-Beta conjugate model (Bayesian inference). 
  
  update_bernoulli_beta <- function(data, alpha0 = 1, beta0 = 1){
    n_non_consumer = data[['non_consumer_sample_size']]
    n_consumer = data[['consumer_sample_size']]
    
    alpha_n = alpha0 + n_consumer
    beta_n = beta0 + n_non_consumer
    
    output = list(data = data, param = list(prior = list(alpha = alpha0, beta = beta0),
                                            posterior = list(alpha = alpha_n, beta = beta_n)))
    
    return(output)
  }

  ################################################################################
  ## post refers to the parameters of the posterior distribution -- normal-gamma 
  ## niter_ale refers to the number of generated samples  
  
  propagate_normal_gamma <- function(niter_ale, post){
  
    precision <- rgamma(1, shape = post$param$posterior$alpha, rate = post$param$posterior$beta)  # precision
    sigma_n <- (1/sqrt(precision))  # standard deviation
    mu <- rnorm(1,post$param$posterior$mu, sigma_n / sqrt(post$param$posterior$v))
    gen_sample <- rlnorm(niter_ale, meanlog = mu, sdlog = sigma_n)
  
    output <- list(gen_sample = gen_sample)
  
    return(output)
  }  

# post_update_concentration = lapply(data_assessment$log_concentration_ss_data, update_normal_gamma, mu0 = 1, v0 = 5, alpha0 = 1, beta0 = 1, suff_stat = TRUE)
# gen_data_concentration = lapply(post_update_concentration, propagate_normal_gamma, niter_ale = 1000)

# post_update_consumption = lapply(data_assessment$log_consumption_ss_data, update_normal_gamma, mu0 = 1, v0 = 5, alpha0 = 1, beta0 = 1, suff_stat = TRUE)
# gen_data_consumption = lapply(post_update_consumption, propagate_normal_gamma, niter_ale = 1000)


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
      mean = change_cocoa_dist$Normal$mean
      sd = change_cocoa_dist$Normal$sd
      param = c(mean = mean, sd = sd)
    }
    if(best_fit == "t"){
      location = change_cocoa_dist$Student.t$location 
      scale = change_cocoa_dist$Student.t$scale
      df = change_cocoa_dist$Student.t$df
      param = c(location = location, scale = scale, df = df)
    }
    if(best_fit == "gamma"){
      shape = change_cocoa_dist$Gamma$shape
      rate  = change_cocoa_dist$Gamma$rate
      param = c(shape = shape, rate = rate)
    }
    if(best_fit == "lognormal"){
      mean.log.X  = change_cocoa_dist$Log.normal$mean.log.X 
      sd.log.X = change_cocoa_dist$Log.normal$sd.log.X
      param = c(mean.log.X = mean.log.X, sd.log.X = sd.log.X)
    }
    if(best_fit == "logt"){
      location.log.X = change_cocoa_dist$Log.Student.t$location.log.X
      scale.log.X = change_cocoa_dist$Log.Student.t$scale.log.X
      df.log.X = change_cocoa_dist$Log.Student.t$df.log.X
      param = c(location.log.X = location.log.X, scale.log.X = scale.log.X,  df.log.X =  df.log.X)
    }
    if(best_fit == "beta"){
      shape1 = change_cocoa_dist$Beta$shape1
      shape2 = change_cocoa_dist$Beta$shape2
      param = c(shape1 = shape1, shape2 = shape2)
    }
    output = list(vals = vals, best_fit = best_fit, param = param)
    return(output)
  }
  
  ################################################################################
  ## Combination of uncertainty

  combine_uncertainty <- function(gen_data_concentration, gen_data_consumption, 
                                       gen_data_EKE, threshold = 0.5, niter_ale){
    
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
                                      threshold = 0.5, percentile_ale,
                                      data_concentration, 
                              concentration_mu0 = 3.5, concentration_v0 = 5, concentration_alpha0 = 1, 
                              concentration_beta0 = 1, suff_stat_concentration = TRUE,
                              data_consumption,
                              consumption_mu0 = 0, consumption_v0 = 5, consumption_alpha0 = 1, 
                              consumption_beta0 = 1, suff_stat_consumption = TRUE,
                              consumers_info_sample_size, consumption_event_alpha0 = 1, 
                              consumption_event_beta0 = 1,
                              gen_data_EKE){
    
    nr_products <-  length(data_concentration)
    prob_consumption <- parameters_consumption <- parameters_concentration <- vector('list', nr_products)
    
    # Probability of a child i consumes chocolate product k
    param_consumption = lapply(consumers_info_sample_size, update_bernoulli_beta, alpha0 = consumption_event_alpha0, beta0 = consumption_event_beta0)
    
    for(k in 1:7){
      prob_consumption[[k]] = rbeta(1,shape1 = param_consumption[[k]]$param$posterior$alpha, shape2 = param_consumption[[k]]$param$posterior$beta)
    }
    
    prob_exceed <- rep(0, niter_epi)
    
    post_concentration = lapply(data_concentration, update_normal_gamma, mu0 = concentration_mu0,
                                       v0 = concentration_v0, alpha0 = concentration_alpha0, 
                                       beta0 = concentration_beta0, suff_stat = suff_stat_concentration)
    
    post_consumption = lapply(data_consumption, update_normal_gamma, mu0 = consumption_mu0, 
                                     v0 = consumption_v0, alpha0 = consumption_alpha0, 
                                     beta0 = consumption_beta0, suff_stat = suff_stat_consumption)
    for(j in 1:nr_products){
      parameters_concentration[[j]] = post_concentration[[j]]$param
      
      parameters_consumption[[j]] = post_consumption[[j]]$param
    }
   
    for(i in 1:niter_epi){
      
      gen_data_concentration = lapply(post_concentration, propagate_normal_gamma, niter_ale = niter_ale)
      
      gen_data_consumption = lapply(post_consumption, propagate_normal_gamma, niter_ale = niter_ale)
      
      prob_exceed[i] <- combine_uncertainty(gen_data_concentration =  gen_data_concentration, gen_data_consumption = gen_data_consumption, 
                                 gen_data_EKE = gen_data_EKE, threshold =  threshold, niter_ale = niter_ale)
          
    }
    
    expected_prob_exceed <- mean(prob_exceed)
    hdi_prob_exceed <- hdi(prob_exceed, credMass = 0.95) # Highest (Posterior) Density Interval
    #**US to IR: is this correct? 
    percentile_prob_exceed = quantile(prob_exceed, probs = (percentile_ale/100)) 
    
    return(list(prob_consumption_event = prob_consumption,
                parameters_concentration = parameters_concentration,
                parameters_consumption = parameters_consumption,
                prob_exceed = prob_exceed, 
                expected_prob_exceed = expected_prob_exceed,
                hdi_prob_exceed = hdi_prob_exceed,
                percentile_prob_exceed = percentile_prob_exceed))
  }
  
 ###################################################################################
  ## Figures  
  graph_plot_pp = function(assessment_output){
    w <- 0
    plot(sort(assessment_output$prob_exceed),(1:length(assessment_output$prob_exceed))/(length(assessment_output$prob_exceed)), 
         xlab = 'Frequency of exceeding the TWI', main = 'Uncertainty', ylab = 'cdf', type = 'l')
    expected_prob_exceed <- assessment_output$expected_prob_exceed
    hdi_prob_exceed <-  assessment_output$hdi_prob_exceed # Highest (Posterior) Density Interval
    #segments(x0 = hdi_prob_exceed[1], x1 = hdi_prob_exceed[2], y0 = w, col = 'blue', lwd = 10)
    #points(expected_prob_exceed, w, pch = 20, cex = 2, col = 'darkred')
    #legend('bottomright',c('Probability interval','Expected value'), col = c('blue','darkred'), pch = c(1,1), bty = 'n')
  }
  
}
#########################################################################################################################################
#########################################################################################################################################
## Precise probability

## Generate one eke sample 
  fit <- quantify_uncertainty_pp_change_eke(vals = data_assessment$change_cons$vals, probs = data_assessment$change_cons$probs/100)
  gen_eke <- rbeta(1,shape1 = fit$param[1], shape2 = fit$param[2])

## Final assessment
  
  TWI_pp <-  unc_analysis_assessment (niter_ale = 10000, niter_epi = 10000, data_concentration = data_assessment$log_concentration_ss_data, 
                                      concentration_mu0 = 3.5, concentration_v0 = 5, concentration_alpha0 = 1,
                                      concentration_beta0 = 1, suff_stat_concentration = TRUE,
                                      data_consumption = data_assessment$log_consumption_ss_data,
                                      consumption_mu0 = 0, consumption_v0 = 5, consumption_alpha0 = 1, 
                                      consumption_beta0 = 1, suff_stat_consumption = TRUE,
                                      consumers_info_sample_size = data_assessment$consumers_info_sample_size, consumption_event_alpha0 = 1, consumption_event_beta0 = 1,
                                      gen_data_EKE = gen_eke, threshold = 0.5, percentile_ale = 95)
  
  
  save(TWI_pp, file='TWI_pp.Rdata')
  
##########################################################################################
  ## Bounded probability 
 
  ### initial_mu0 = c(concentration_mu0, consumption_mu0)

  obj_func <- function(parameters, niter_ale = 1000, niter_epi = 1000, 
                       data_concentration = data_assessment$log_concentration_ss_data, 
                       concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                       suff_stat_concentration = TRUE,
                       data_consumption = data_assessment$log_consumption_ss_data,
                       consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                       suff_stat_consumption = TRUE,
                       consumers_info_sample_size = data_assessment$consumers_info_sample_size, 
                       consumption_event_alpha0 = 1, consumption_event_beta0 = 1,
                       gen_data_EKE = gen_eke, threshold = 0.5, percentile_ale, output_exp = 'TRUE'){
    
    concentration_mu0 <- parameters[1] 
    consumption_mu0 <- parameters[2] 
    
    out <- unc_analysis_assessment(niter_ale = niter_ale, niter_epi= niter_epi, data_concentration = data_concentration, 
                              concentration_mu0 = concentration_mu0, concentration_v0 = concentration_v0, 
                              concentration_alpha0 = concentration_alpha0, concentration_beta0 = concentration_beta0,
                              suff_stat_concentration = suff_stat_concentration,
                              data_consumption = data_consumption,
                              consumption_mu0 = consumption_mu0, consumption_v0 =  consumption_v0, 
                              consumption_alpha0 = consumption_alpha0, consumption_beta0 = consumption_beta0, 
                              suff_stat_consumption = suff_stat_consumption,
                              consumers_info_sample_size = consumers_info_sample_size, 
                              consumption_event_alpha0 = consumption_event_alpha0, consumption_event_beta0 = consumption_event_beta0,
                              gen_data_EKE = gen_data_EKE, threshold = threshold, percentile_ale = percentile_ale)
    
    if(output_exp == 'TRUE'){
      output <- out$expected_prob_exceed
    }
    else{
      output <- out$percentile_prob_exceed[[1]]
    }
    return(output)
  }
 
  
  graph_bp = function(lower_points, upper_points, min_exp, max_exp){
    
    l = length(lower_points)
    l_points = rep(0,l)
    u_points = rep(0,l)
    
    for(i in 1:l){
      l_points[i] = min(lower_points[i],upper_points[i]) 
      u_points[i] = max(lower_points[i],upper_points[i])
    }
    
    par(mfrow = c(1,1))
    plot(sort(l_points),(1:l)/(l), xlim = c(0,max(u_points)),
         xlab = 'Frequency of exceeding the TWI', main = 'Uncertainty', ylab = 'cdf', type = 'l', col = 'red')
    par(new=TRUE)
    
    plot(sort(u_points),(1:l)/(l), xlim = c(0,max(u_points)),
         xlab = 'Frequency of exceeding the TWI', main = 'Uncertainty', ylab = 'cdf', type = 'l', col = 'blue')
    par(new=TRUE)
    segments(x0 = sort(l_points)[1], x1 = sort(u_points)[1], y0 = 0, col = 'blue', lwd = 1) 
    segments(x0 = tail(sort(l_points),1), x1 = tail(sort(u_points),1), y0 = 1, col = 'red', lwd = 1) 
    
    #segments(x0 = min_exp, x1 = max_exp, y0 = 0, col = 'black', lwd = 2) 
    
    legend('bottomright',c('Lower cdf', 'Upper cdf'), col = c('blue','red'), pch = c(1,1), bty = 'n')
    
  }
  
#####################################################################
  ### initial_mu0 = c(concentration_mu0, consumption_mu0)
  initial_mu0 = c(5.75, 0.75)
  
  argmin_upper_exp = nmkb(par = initial_mu0, fn = obj_func,  lower = c(1, -5), upper = c(6, 1),
                    control = list(maximize = TRUE),
                    niter_ale = 1000, niter_epi = 1000, 
                    data_concentration = data_assessment$log_concentration_ss_data, 
                    concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                    suff_stat_concentration = TRUE,
                    data_consumption = data_assessment$log_consumption_ss_data,
                    consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                    suff_stat_consumption = TRUE,
                    consumers_info_sample_size = data_assessment$consumers_info_sample_size, 
                    consumption_event_alpha0 = 1, consumption_event_beta0 = 1,
                    gen_data_EKE = gen_eke, threshold = 0.5, percentile_ale = 95, output_exp = 'TRUE')
   
   argmin_upper_perc = nmkb(par = initial_mu0, fn = obj_func,  lower = c(1, -5), upper = c(6, 1),
                    control = list(maximize = TRUE),
                    niter_ale = 1000, niter_epi = 1000, 
                    data_concentration = data_assessment$log_concentration_ss_data, 
                    concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                    suff_stat_concentration = TRUE,
                    data_consumption = data_assessment$log_consumption_ss_data,
                    consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                    suff_stat_consumption = TRUE,
                    consumers_info_sample_size = data_assessment$consumers_info_sample_size, 
                    consumption_event_alpha0 = 1, consumption_event_beta0 = 1,
                    gen_data_EKE = gen_eke, threshold = 0.5, percentile_ale = 95, output_exp = 'FALSE')
   
   ############
   ### initial_mu0 = c(concentration_mu0, consumption_mu0)
   initial_mu0 = c(2.75, -2.5)
   
   argmin_lower_exp = nmkb(par = initial_mu0, fn = obj_func,  lower = c(1, -5), upper = c(6, 1),
                           control = list(maximize = FALSE),
                           niter_ale = 1000, niter_epi = 1000, 
                           data_concentration = data_assessment$log_concentration_ss_data, 
                           concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                           suff_stat_concentration = TRUE,
                           data_consumption = data_assessment$log_consumption_ss_data,
                           consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                           suff_stat_consumption = TRUE,
                           consumers_info_sample_size = data_assessment$consumers_info_sample_size, 
                           consumption_event_alpha0 = 1, consumption_event_beta0 = 1,
                           gen_data_EKE = gen_eke, threshold = 0.5, percentile_ale = 95, output_exp = 'TRUE')
   

   argmin_lower_perc = nmkb(par = initial_mu0, fn = obj_func,  lower = c(1, -5), upper = c(6, 1),
                            control = list(maximize = FALSE),
                            niter_ale = 1000, niter_epi = 1000, 
                            data_concentration = data_assessment$log_concentration_ss_data, 
                            concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                            suff_stat_concentration = TRUE,
                            data_consumption = data_assessment$log_consumption_ss_data,
                            consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                            suff_stat_consumption = TRUE,
                            consumers_info_sample_size = data_assessment$consumers_info_sample_size, 
                            consumption_event_alpha0 = 1, consumption_event_beta0 = 1,
                            gen_data_EKE = gen_eke, threshold = 0.5, percentile_ale = 95, output_exp = 'FALSE')
   
   #################
   percentiles = seq(1,99, by = 1)

   initial_mu0 = c(5.75, 0.75)
   
   upper_bound_perc <- lapply(c(1,5,25,50,75,95,99), nmkb, par = initial_mu0, fn = obj_func,  lower = c(1, -5), upper = c(6, 1),
                        control = list(maximize = TRUE),
                        niter_ale = 1000, niter_epi = 1000, 
                        data_concentration = data_assessment$log_concentration_ss_data, 
                        concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                        suff_stat_concentration = TRUE,
                        data_consumption = data_assessment$log_consumption_ss_data,
                        consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                        suff_stat_consumption = TRUE,
                        consumers_info_sample_size = data_assessment$consumers_info_sample_size, 
                        consumption_event_alpha0 = 1, consumption_event_beta0 = 1,
                        gen_data_EKE = gen_eke, threshold = 0.5, output_exp = 'FALSE')
   
   save(upper_bound_perc, file = 'upper_bound_perc_bp.Rdata')
   
   percentiles_1 <- c(1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99)
   
   upper_bound_perc1 <- lapply(percentiles_1, nmkb, par = initial_mu0, fn = obj_func,  lower = c(1, -5), upper = c(6, 1),
                            control = list(maximize = TRUE),
                            niter_ale = 1000, niter_epi = 1000, 
                            data_concentration = data_assessment$log_concentration_ss_data, 
                            concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                            suff_stat_concentration = TRUE,
                            data_consumption = data_assessment$log_consumption_ss_data,
                            consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                            suff_stat_consumption = TRUE,
                            consumers_info_sample_size = data_assessment$consumers_info_sample_size, 
                            consumption_event_alpha0 = 1, consumption_event_beta0 = 1,
                            gen_data_EKE = gen_eke, threshold = 0.5, output_exp = 'FALSE')
   
   save(upper_bound_perc1, file = 'upper_bound_perc1_bp.Rdata')
   ####
   
   initial_mu0 = c(2.75, -2.5)
   
   lower_bound_perc <- lapply(c(1,5,25,50,75,95,99), nmkb, par = initial_mu0 , fn = obj_func,  lower = c(1, -5), upper = c(6, 1),
                        control = list(maximize = FALSE),
                        niter_ale = 1000, niter_epi = 1000, 
                        data_concentration = data_assessment$log_concentration_ss_data, 
                        concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                        suff_stat_concentration = TRUE,
                        data_consumption = data_assessment$log_consumption_ss_data,
                        consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                        suff_stat_consumption = TRUE,
                        consumers_info_sample_size = data_assessment$consumers_info_sample_size, 
                        consumption_event_alpha0 = 1, consumption_event_beta0 = 1,
                        gen_data_EKE = gen_eke, threshold = 0.5, output_exp = 'FALSE')
  
   save(lower_bound_perc , file = 'lower_bound_perc_bp.Rdata')
   
   percentiles_1 <- c(1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99)
   
   lower_bound_perc1 <- lapply(percentiles_1, nmkb, par = initial_mu0, fn = obj_func,  lower = c(1, -5), upper = c(6, 1),
                            control = list(maximize = FALSE),
                            niter_ale = 1000, niter_epi = 1000, 
                            data_concentration = data_assessment$log_concentration_ss_data, 
                            concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                            suff_stat_concentration = TRUE,
                            data_consumption = data_assessment$log_consumption_ss_data,
                            consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                            suff_stat_consumption = TRUE,
                            consumers_info_sample_size = data_assessment$consumers_info_sample_size, 
                            consumption_event_alpha0 = 1, consumption_event_beta0 = 1,
                            gen_data_EKE = gen_eke, threshold = 0.5, output_exp = 'FALSE')
   
   save(lower_bound_perc1, file = 'lower_bound_perc1_bp.Rdata')
   
   ## extract the values corresponding to percentiles_1
   
   upper_vals <- unlist(lapply(upper_bound_perc1, function(x){x$value}))
   lower_vals <- unlist(lapply(lower_bound_perc1, function(x){x$value}))
    
  ## Visualizations
load('TWI_pp.Rdata')
graph_plot_pp(assessment_output = TWI_pp)
  
   
load('upper_bound_perc_bp.Rdata')
load('upper_bound_perc1_bp.Rdata')
load('lower_bound_perc_bp.Rdata')
load('lower_bound_perc1_bp.Rdata')
graph_bp...