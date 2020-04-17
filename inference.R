######
# Code for paper comparing precise to bounded probability 
# Ivette Raices Cruz, Matthias C M Troffaes, Ullrika Sahlin
# 17 april 2020

########################
library('SHELF')
library('HDInterval')
library('dfoptim')

####################################################
## Load data 
load('EKE_change_cons.Rdata')
load('log_concentration_ss_data.Rdata')
load('log_consumption_ss_data.Rdata')
load('consumers_info_sample_size.Rdata')

###################################################
## Parametric inference quantifying epistemic uncertainty by bounded probability

## Functions 
{
  ## https://en.wikipedia.org/wiki/Conjugate_prior
  
  ## Normal-Gamma conjugate model (Bayesian inference). 
  ## Bounded probability takes into account a lower and upper bound for the prior mean.
  
  update_normal_gamma_suff_stat <- function(data, mu0 = 0, v0 = 5, alpha0 = 1, beta0 = 1){
    n <- data[['sample_size']]
    m <- data[['sample_mean']]
    S2 <- data[['sample_sd']]^2  # sample variance
    
    v_n <-  v0 + n
    alpha_n <- alpha0 + n/2
    
    mu_n <- (v0 * mu0 + n * m ) / (v0 + n)
    beta_n <- beta0 + 1/2*S2*(n - 1) + (n*v0/(v0 + n) * (m - mu0)^2)
      
    output = list(data = data, param = list(prior = list(mu = mu0, v = v0, alpha = alpha0, beta = beta0),
                                            posterior = list(mu = mu_n, v = v_n , alpha = alpha_n, beta = beta_n)))
      
    return(output)
  }
  
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
  ## sufficient statistics (sample mean and sample variance) for the normal likelihood are known for log-aluminium concentration and log-consumption 
  
  quantify_uncertainty_normal <- function(niter_ale, data, mu0 = 1, v0 = 5, alpha0 = 1, beta0 = 1){
       fun <- function(k){
         post <- update_normal_gamma_suff_stat(data = data[[k]], mu0 = mu0, v0 = v0, alpha0 = alpha0, beta0 = beta0)
          
         precision <- rgamma(1, shape = post$param$posterior$alpha, rate = post$param$posterior$beta)  # precision
         sigma_n <- (1/sqrt(precision))  # standard deviation
         mu <- rnorm(1,post$param$posterior$mu, sigma_n / sqrt(post$param$posterior$v))
         gen_sample <- rlnorm(niter_ale, meanlog = mu, sdlog = sigma_n)
         
         list(post_mu = post$param$posterior$mu, post_v = post$param$posterior$v, post_alpha = post$param$posterior$alpha, post_beta = post$param$posterior$beta,  
                 gen_sample)
    }
    gen <- lapply(1:length(data),fun)
    return(gen)
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
  ## Propagation of uncertainty
  
  propagate_uncertainty <- function(gen_data_concentration, gen_data_consumption, 
                                       gen_data_EKE, th = 0.5, niter_ale){
    
    l_concentration <- length(gen_data_concentration[[1]]) #5 
    l_consumption <- length(gen_data_consumption[[1]])  #5
     
    l_conc <- length(gen_data_concentration) #7
    l_cons <- length(gen_data_consumption)  #7
    
    gen_conc  <- vector('list',l_conc) 
    gen_cons  <- vector('list',l_cons)
    parameters_concentration <- vector('list',l_conc)
    parameters_consumption <- vector('list',l_cons)
    
    for(i in 1:l_conc){
      gen_conc[[i]] <- gen_data_concentration[[i]][5]
      gen_cons[[i]] <- gen_data_consumption[[i]][5]
      parameters_concentration[[i]] <- gen_data_concentration[[i]][1:4]
      parameters_consumption[[i]] <- gen_data_consumption[[i]][1:4]
    }
    
    gen_concentration <- matrix(unlist(gen_conc), ncol = 7, nrow = niter_ale)
    gen_consumption <- matrix(unlist(gen_cons), ncol = 7, nrow = niter_ale)
        
    weekly_intake <- rowSums(gen_data_EKE * (gen_consumption * 0.007) * gen_concentration)
    
    prob_exceed_wi <- mean(weekly_intake > th)
   
    return(list(parameters_concentration = parameters_concentration , parameters_consumption = parameters_consumption, 
           prob_exceed_wi = prob_exceed_wi))
  }
  
  
  ###############################################################################
  ## Final assessment model -- Precise Probability 
  
  do_assessment <- function(niter_ale = 1000, niter_epi = 1000, data_concentration, 
                              concentration_mu0 = 3.5, concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                              data_consumption,
                              consumption_mu0 = 0, consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                              gen_data_EKE, th = 0.5, percentile){
    
    
    l_data <- length(data_concentration)  #7
    parameters_concentration <- vector('list', l_data)
    parameters_consumption <- vector('list', l_data)
  
    prob_exceed <- rep(0, niter_epi)
       
    for(i in 1:niter_epi){
        gen_concentration <- quantify_uncertainty_normal(niter_ale = niter_ale, data = data_concentration, 
                                                          mu0 = concentration_mu0, v0 = concentration_v0, 
                                                          alpha0 = concentration_alpha0, beta0 = concentration_beta0)
          
        gen_consumption <- quantify_uncertainty_normal(niter_ale = niter_ale, data =  data_consumption, 
                                                            mu0 = consumption_mu0, v0 = consumption_v0, 
                                                            alpha0 = consumption_alpha0, beta0 = consumption_beta0)
          
        aux <- propagate_uncertainty(gen_data_concentration =  gen_concentration, gen_data_consumption = gen_consumption, 
                                     gen_data_EKE = gen_data_EKE, th = th, niter_ale = niter_ale)
          
        prob_exceed[i] <- aux[[3]]
          
    }
    
    expected_prob_exceed <- mean(prob_exceed)
    hdi_prob_exceed <- hdi(prob_exceed, credMass = 0.95) # Highest (Posterior) Density Interval
    percentile_prob_exceed = quantile (prob_exceed, probs = (percentile/100)) 
    
    
    parameters_concentration = aux[[1]][1:l_data]
    parameters_consumption = aux[[2]][1:l_data]

    return(list(parameters_concentration = parameters_concentration,
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
         xlab = 'Likelihood to exceed TWI', main = 'Uncertainty', ylab = 'cdf', type = 'l')
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

## Probability of a child i consumes chocolate product k
param_consumption = lapply(consumers_info_sample_size, update_bernoulli_beta, alpha0 = 1, beta0 = 1)

prob_consumption = list(vector, 7)
for(k in 1:7){
  prob_consumption[[k]] = rbeta(1,shape1 = param_consumption[[k]]$param$posterior$alpha, shape2 = param_consumption[[k]]$param$posterior$beta)
}
prob_consumption  

## Generate one eke sample 
  fit <- quantify_uncertainty_pp_change_eke(vals = change.cons$vals, probs = change.cons$probs/100)
  gen_eke <- rbeta(1,shape1 = fit$param[1], shape2 = fit$param[2])

## Final assessment
  
  TWI_pp <- do_assessment(niter_ale = 1000, niter_epi= 1000, data_concentration = log_concentration_ss_data, 
                             concentration_mu0 = 3.5, concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                             data_consumption = log_consumption_ss_data,
                             consumption_mu0 = 0, consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                             gen_data_EKE = gen_eke, th = 0.5, percentile = 95)
  
  
  graph_plot_pp(assessment_output = TWI_pp)
  
##########################################################################################
  ## Bounded probability 
 
  ### parameters = c(concentration_mu0, consumption_mu0)

  obj_func_exp <- function(parameters, niter_ale = 1000, niter_epi = 1000, 
                       data_concentration = log_concentration_ss_data, 
                       concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                       data_consumption = log_consumption_ss_data,
                       consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                       gen_data_EKE = gen_eke, th = 0.5, percentile){
    
    concentration_mu0 <- parameters[1] 
    consumption_mu0 <- parameters[2] 
    
    out <- final_assessment(niter_ale = 1000, niter_epi= 1000, data_concentration = log_concentration_ss_data, 
                              concentration_mu0 = concentration_mu0, concentration_v0 = concentration_v0, 
                              concentration_alpha0 = concentration_alpha0, concentration_beta0 = concentration_beta0, 
                              data_consumption = log_consumption_ss_data,
                              consumption_mu0 = consumption_mu0, consumption_v0 =  consumption_v0, 
                              consumption_alpha0 = consumption_alpha0, consumption_beta0 = consumption_beta0, 
                              gen_data_EKE = gen_data_EKE, th = th, percentile = percentile)
      
    expected_prob_exceed <- out$expected_prob_exceed
    #percentile_prob_exceed <- out$percentile_prob_exceed[[1]]
    
    return(expected_prob_exceed)
  }
 
  
  obj_func_perc <- function(parameters, niter_ale = 1000, niter_epi = 1000, 
                       data_concentration = log_concentration_ss_data, 
                       concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                       data_consumption = log_consumption_ss_data,
                       consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                       gen_data_EKE = gen_eke, th = 0.5, percentile){
    
    concentration_mu0 <- parameters[1] 
    consumption_mu0 <- parameters[2] 
    
    out <- final_assessment(niter_ale = 1000, niter_epi= 1000, data_concentration = log_concentration_ss_data, 
                              concentration_mu0 = concentration_mu0, concentration_v0 = concentration_v0, 
                              concentration_alpha0 = concentration_alpha0, concentration_beta0 = concentration_beta0, 
                              data_consumption = log_consumption_ss_data,
                              consumption_mu0 = consumption_mu0, consumption_v0 =  consumption_v0, 
                              consumption_alpha0 = consumption_alpha0, consumption_beta0 = consumption_beta0, 
                              gen_data_EKE = gen_data_EKE, th = th, percentile = percentile)
      
    #expected_prob_exceed <- out$expected_prob_exceed
    percentile_prob_exceed <- out$percentile_prob_exceed[[1]]
    
    return(percentile_prob_exceed)
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
         xlab = 'Likelihood of exceed TWI', main = 'Uncertainty', ylab = 'cdf', type = 'l', col = 'red')
    par(new=TRUE)
    
    plot(sort(u_points),(1:l)/(l), xlim = c(0,max(u_points)),
         xlab = 'Likelihood of exceed TWI', main = 'Uncertainty', ylab = 'cdf', type = 'l', col = 'blue')
    par(new=TRUE)
    segments(x0 = sort(l_points)[1], x1 = sort(u_points)[1], y0 = 0, col = 'blue', lwd = 1) 
    segments(x0 = tail(sort(l_points),1), x1 = tail(sort(u_points),1), y0 = 1, col = 'red', lwd = 1) 
    
    #segments(x0 = min_exp, x1 = max_exp, y0 = 0, col = 'black', lwd = 2) 
    
    legend('bottomright',c('Lower cdf', 'Upper cdf'), col = c('blue','red'), pch = c(1,1), bty = 'n')
    
  }
  
#####################################################################
   
  initial_mu0 = c(5.75, 0.75)
  
  argmin_max_mean = nmkb(par = initial_mu0, fn = obj_func_mean,  lower = c(1, -5), upper = c(6, 1),
                    control = list(maximize = TRUE),
                    niter_ale = 1000, niter_epi = 1000, 
                    data_concentration = log_concentration_ss_data, 
                    concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                    data_consumption = log_consumption_ss_data,
                    consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                    gen_data_EKE = gen_eke, th = 0.5, percentile = 95)
   
   argmin_max_perc = nmkb(par = initial_mu0, fn = obj_func_perc,  lower = c(1, -5), upper = c(6, 1),
                    control = list(maximize = TRUE),
                    niter_ale = 1000, niter_epi = 1000, 
                    data_concentration = log_concentration_ss_data, 
                    concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                    data_consumption = log_consumption_ss_data,
                    consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                    gen_data_EKE = gen_eke, th = 0.5, percentile = 95)
   
   initial_mu0 = c(2.75, -2.5)
   
   argmin_min_mean = nmkb(par = parameters, fn = obj_func_mean,  lower = c(1, -5), upper = c(6, 1),
                    niter_ale = 1000, niter_epi = 1000, 
                    data_concentration = log_concentration_ss_data, 
                    concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                    data_consumption = log_consumption_ss_data,
                    consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                    gen_data_EKE = gen_eke, th = 0.5, percentile = 95)
   

   argmin_min_perc = nmkb(par = parameters, fn = obj_func_perc,  lower = c(1, -5), upper = c(6, 1),
                    niter_ale = 1000, niter_epi = 1000, 
                    data_concentration = log_concentration_ss_data, 
                    concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                    data_consumption = log_consumption_ss_data,
                    consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                    gen_data_EKE = gen_eke, th = 0.5, percentile = 95)
   
   #################
   percentiles = seq(1,99, by = 1)

   parameters = c(5.75, 0.75)
   
   upper_bound_perc <- lapply(c(1,5,25,50,75,95,99), nmkb, par = parameters, fn = obj_func_perc,  lower = c(1, -5), upper = c(6, 1),
                        control = list(maximize = TRUE),
                        niter_ale = 1000, niter_epi = 1000, 
                        data_concentration = log_concentration_ss_data, 
                        concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                        data_consumption = log_consumption_ss_data,
                        consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                        gen_data_EKE = gen_eke, th = 0.5)
   
   save(max_values, file = 'max_values_bp.Rdata')
   
   upper_bound_perc <- lapply(c(1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99), nmkb, par = parameters, fn = obj_func_perc,  lower = c(1, -5), upper = c(6, 1),
                            control = list(maximize = TRUE),
                            niter_ale = 1000, niter_epi = 1000, 
                            data_concentration = log_concentration_ss_data, 
                            concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                            data_consumption = log_consumption_ss_data,
                            consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                            gen_data_EKE = gen_eke, th = 0.5)
   
   save(max_values_per, file = 'max_values_percentiles_bp.Rdata')
   
   ####
   
   parameters = c(2.75, -2.5)
   
   min_values <- lapply(c(1,5,25,50,75,95,99), nmkb, par = parameters, fn = obj_func,  lower = c(1, -5), upper = c(6, 1),
                        control = list(maximize = FALSE),
                        niter_ale = 1000, niter_epi = 1000, 
                        data_concentration = log_concentration_ss_data, 
                        concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                        data_consumption = log_consumption_ss_data,
                        consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                        gen_data_EKE = gen_eke, th = 0.5)
  
   save(min_values, file = 'min_values_bp.Rdata')
   
   
   min_values_per <- lapply(c(1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,99), nmkb, par = parameters, fn = obj_func,  lower = c(1, -5), upper = c(6, 1),
                            control = list(maximize = FALSE),
                            niter_ale = 1000, niter_epi = 1000, 
                            data_concentration = log_concentration_ss_data, 
                            concentration_v0 = 5, concentration_alpha0 = 1, concentration_beta0 = 1, 
                            data_consumption = log_consumption_ss_data,
                            consumption_v0 = 5, consumption_alpha0 = 1, consumption_beta0 = 1, 
                            gen_data_EKE = gen_eke, th = 0.5)
   
   save(min_values_per, file = 'min_values_percentiles_bp.Rdata')
   
   load('min_values_percentiles_bp.Rdata')
   
  min_valss <- unlist(lapply(min_values_per, function(x){x$value}))
   