######################################################################################################
#  Key Functions that are Used in the Report
######################################################################################################

# results = tibble(code = codes, 
#                  rmse = rep(NA, N_tests), # Root mean square of ultimate cumulative claims
#                  rmse_loss_ratio  = NA, # Root mean square of ultimate loss ratio
#                  rmse_loss_ratio_w  = NA # Root mean square of weighted ultimate loss ratio
# ) 


source('loss_reserve.R')
#Load Data
data = load_losses_subset()

# Parameters
# 
# iter = 500
# warmup = 400
# chains = 6
# adapt_delta = 0.9
# max_treedepth = 10
# business_lines = c("comauto", "medmal", "othliab","ppauto","prodliab", "wkcomp")



#' Quantile Function for Kolmogorov Distribution
#' 
#' @param p p value 
#' @param N Number of samples
#' @return Quantile for the Kolmogorov distribution
qkolm = function(p, N){
  diff = function(D, target_p, N){
    actual_p = pkolm(d=D, n = N)
    (actual_p-target_p)^2
  }
  ans = optim(par = 1.36/sqrt(N), fn = diff, target_p = p, N=N, method = "L-BFGS-B" , lower=0.00000001)
  
  if (ans$convergence !=0){
    print("Convergence Values")
    print(ans)
    warning("qkolm did not converge.  optim results:")
  }
  ans$par
}

#' Performs a Kolmogorov Smirnov Test on a set of percentiles and produces a PP Plot
#' 
#' @param p_i A vector of percentiles
#' @param p_level The p level at which to reject the Null Hypothesis for a Kolmogorov-Smirnov test of uniformity of the p_i values
#' @param title Title to include in graphs
#' @param na.rm if TRUE then remove NA values from p_i
#' @param breaks Breaks used for histogram plot
#' @return A list containing: (i) a PP Plot of the p_i values with the rejection levels at the specified p_level, (ii) K-S test results
ks_percentile = function(p_i, p_level = 0.05, title = "", na.rm = TRUE, breaks=0:5/5){
  if (na.rm){
    p_i  =p_i[!is.na(p_i)]
  }
  
  N=length(p_i)
  p_i=sort(p_i)
  f_i=1:N/N
  D = max(abs(p_i - f_i))
  p_value =  1-pkolm(d=D, n = N)
  
  data = data.frame(Theoretical = f_i, Actual = p_i)
  rejection_level = qkolm(p = 1-p_level, N = N)
  
  pp_plot = ggplot(data) + 
    theme_light()+
    geom_point(aes(Theoretical, Actual), color = "blue") + 
    geom_abline(slope = 1, intercept=0, color = "red") +
    geom_abline(slope = 1, intercept=-rejection_level, color = "blue", linetype = 2) +
    geom_abline(slope = 1, intercept=+rejection_level, color = "blue", linetype = 2) +
    labs(title = paste0(title, " PP Plot: \n p-level:", p_level))
  hist_plot = ggplot(data) + 
    theme_light()+ 
    labs(title = paste0(title, " Percentile Distribution"), x="Percentiles", y="Count")+
    geom_histogram(aes(p_i), fill="lightblue", color="black", breaks=breaks)
  hist_plot
  return(list(pp_plot= pp_plot, hist_plot = hist_plot, rejection_level = rejection_level, D= D, reject_hypothesis = D >=rejection_level, p_value = p_value))
}




#' Calculate the predicted loss ratio and supporting analytics using an incremental loss ratio model
#' 
#' Uses a vector for mu prior
#' 
#' @param loss_triangle the Loss triangle and loss square
#' @param premiums Vector of premiums
#' @param calc_ranks If TRUE then additional statistics are calculated
#' @param formula_x formula for square exponential 
#' @param formula_H formula for linear terms in kernel (not in mean function)
#' @param mu The vector for the prior mu values for the GP
#' @param prior_eta, prior_rho, prior_theta, prior_sigma, prior_alpha, prior_beta Priors for various parameters in the Stan model.  If a model does not require one of these priors, then such prior is ignored.
#' @param take_logs If TRUE, then response is logarithm is used as the predictor.  Note, if TRUE, then LHS of formula_x must be y_log_scaled, otherwise use y_scaled
#' @param plot_graphs If TRUE, then sample investigative graphs are plotted.
#' @param virtual_data if TRUE, then a training column for duration 11 is added with incremental loss = 0%
#' @return If calc_ranks = FALSE, then only the predicted loss_square is returned.  If TRUE, then the percentage ranks, sigma estimates are also returned.
stan_incremental_loss_ratio3 = function(loss_triangle, premiums=NULL, loss_square = NULL, calc_ranks = FALSE, verbose=TRUE,
                                        formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                        formula_H = y_scaled ~  AccidentYear_scaled + DevelopmentLag_scaled -1 , formula_cat = NULL,
                                        eps = 0.001, plot_graphs = FALSE,
                                        chains = 6,iter = 600, warmup =300, adapt_delta = 0.9,max_treedepth = 10,
                                        prior_eta = 1, prior_rho = c(5.566019, 5.570317), prior_theta = 1, prior_sigma = 1, 
                                        prior_beta = 1, take_logs = FALSE, 
                                        scaled_lower_bound = NULL, l_bound = 0, mu=10:1/10-.5,
                                        virtual_data = FALSE, virtual_point = -.3, virtual_accident_year = 1988:1997, virtual_development_lag = 11,
                                        stanFile = "gp_compound_mean_03.stan", use_vb = FALSE){
  
  


  if(is.null(loss_square)) 
    loss_square = matrix(NA, nrow(loss_triangle), ncol(loss_triangle))
  
  
  AccidentYear_unique = AccidentYear = as.numeric(rownames(loss_triangle))
  DevelopmentLag_unique = DevelopmentLag = as.numeric(colnames(loss_triangle))
  #DevelopmentLag = DevelopmentLag[-1]
  no_acc = length(AccidentYear); 
  no_dev = length(DevelopmentLag)
  
  AccidentYear = rep(AccidentYear, each = no_dev)
  DevelopmentLag = rep(DevelopmentLag, times = no_acc)
  
  # if (virtual_data){
  #   no_dev_train = no_dev + 1
  # } else {
  #   no_dev_train = no_dev
  # }
  no_dev_train = no_dev
  
  loss_ratio_triangle = (loss_triangle / premiums)
  loss_ratio_triangle = loss_ratio_triangle - cbind(0,loss_ratio_triangle[,-ncol(loss_ratio_triangle)])
  loss_ratio_triangle_log = log(apply((loss_triangle / premiums),2,function(x) pmax(x,eps)))
  loss_ratio_triangle_log = loss_ratio_triangle_log - cbind(0,loss_ratio_triangle_log[,-ncol(loss_ratio_triangle_log)])
  
  # if (virtual_data){
  #   loss_ratio_triangle = cbind(loss_ratio_triangle, 0)
  #   loss_ratio_triangle_log = cbind(loss_ratio_triangle_log, 0)
  #   colnames(loss_ratio_triangle) = 2:(no_dev_train+1)
  #   colnames(loss_ratio_triangle_log) = 2:(no_dev_train+1)
  # }
  
  
  loss_ratio = as_tibble(loss_ratio_triangle) %>%
    mutate(AccidentYear = as.integer(rownames(loss_ratio_triangle))) %>%
    gather(DevelopmentLag, loss_ratio, -AccidentYear) %>%
    filter(!is.na(loss_ratio)) %>%
    mutate(
      DevelopmentLag = as.integer(DevelopmentLag)
    )
  loss_ratio_log = as_tibble(loss_ratio_triangle_log) %>%
    mutate(AccidentYear = as.integer(rownames(loss_ratio_triangle_log))) %>%
    gather(DevelopmentLag, loss_ratio_log, -AccidentYear) %>%
    filter(!is.na(loss_ratio_log)) %>%
    mutate(
      DevelopmentLag = as.integer(DevelopmentLag)
    )
  loss_ratio = loss_ratio %>%
    left_join(loss_ratio_log, by = c("AccidentYear", "DevelopmentLag"))
  
  scaler_accident_year = scaler(loss_ratio$AccidentYear)
  scaler_development_lag = scaler(loss_ratio$DevelopmentLag)
  DevelopmentLag_transform = function(x) log(x-min(loss_ratio$DevelopmentLag)+1)
  scaler_loss_ratio = scaler(loss_ratio$loss_ratio, l_bound = scaled_lower_bound)
  scaler_loss_ratio_log = scaler(loss_ratio$loss_ratio_log, l_bound = scaled_lower_bound)

  if (!take_logs) {
    l_bound = scaler_loss_ratio(l_bound)[1,1]
    mu = c(scaler_loss_ratio(mu))
  } else {
    l_bound = scaler_loss_ratio_log(l_bound)[1,1]
    mu = c(scaler_loss_ratio_log(mu))
  }
  loss_ratio = loss_ratio %>%
    mutate(AccidentYear_scaled = scaler_accident_year(AccidentYear)[,1],
           DevelopmentLag_scaled = scaler_development_lag(DevelopmentLag)[,1], 
           DevelopmentLag_trans = DevelopmentLag_transform(DevelopmentLag),
           y_scaled = scaler_loss_ratio(loss_ratio)[,1],
           y_log_scaled = scaler_loss_ratio_log(loss_ratio_log)[,1]) 
  
  
  # Create inputs for out-of-sample predictions
  loss_ratio_new = tibble(AccidentYear ,DevelopmentLag) %>%
    mutate(AccidentYear_scaled = scaler_accident_year(AccidentYear)[,1],
           DevelopmentLag_scaled = scaler_development_lag(DevelopmentLag)[,1], 
           DevelopmentLag_trans = DevelopmentLag_transform(DevelopmentLag),
           y_scaled = 0,
           y_log_scaled = 0)
  
  formula_x = update(formula_x, . ~ . -1) # remove intercept term
  formula_H = update(formula_H, . ~ . -1) # remove intercept term
  
  if (virtual_data){
    virtual_loss_ratio = tibble(
      AccidentYear = virtual_accident_year,
      DevelopmentLag = virtual_development_lag,
      loss_ratio = virtual_point,
      loss_ratio_log = NA 
    ) %>%
      mutate(
        AccidentYear_scaled = scaler_accident_year(AccidentYear),
        DevelopmentLag_scaled = scaler_development_lag(DevelopmentLag),
        DevelopmentLag_trans = NA,
        y_scaled = scaler_loss_ratio(loss_ratio),
        y_log_scaled = NA 
      )
    loss_ratio = rbind(loss_ratio, virtual_loss_ratio)
    
    
  }
  
  dimSigma = length(unique(loss_ratio$DevelopmentLag))
  if (!all(loss_ratio$DevelopmentLag %in% 1:dimSigma)) {
    stop(paste0("DevelopmentLag values including virtual data points must be in range 1:", dimSigma))
  }
  
  if (length(mu) != dimSigma & length(mu) != 1) {
    stop(paste0("The length of mu must equal dimSigma, namely: ", dimSigma))
  }
  
  stan_list <- list(
    N1 = nrow(loss_ratio),
    N2 = nrow(loss_ratio_new),
    x1 =  model.matrix(formula_x, loss_ratio),
    y1 = c(loss_ratio[,all.vars(update(formula_x, .~0)) ][[1]]),
    dev_lag1 = loss_ratio$DevelopmentLag,
    x2 = model.matrix(formula_x, loss_ratio_new),
    dev_lag2 = loss_ratio_new$DevelopmentLag,
    mu = mu,
    l_bound=l_bound,
    prior_eta=prior_eta,
    prior_rho = prior_rho,
    prior_theta = prior_theta,
    prior_sigma = prior_sigma,
    prior_beta = prior_beta,
    dimSigma = dimSigma
  )
  
  if (take_logs & (all.vars(formula_x)[1] != "y_log_scaled")){
    warning("if logs are used, LHS of formula_x must be 'y_log_scaled'")
  } 
  
  stan_list$D = ncol(stan_list$x1)
  
  if(!is.null(formula_H)){
    stan_list$H1 = t(model.matrix(formula_H, loss_ratio))
    stan_list$H2 = t(model.matrix(formula_H, loss_ratio_new))
    stan_list$H1 = stan_list$H1 - apply(stan_list$H1,1,min)
    stan_list$H2 = stan_list$H2 - apply(stan_list$H1,1,min)
    stan_list$dimH <- nrow(stan_list$H1)
  }
  
  
  if(!is.null(formula_cat)){
    formula_cat = update(formula_cat,  ~ . -1) # remove intercept term
    stan_list$x_cat1 = c(model.matrix(formula_cat, loss_ratio))
    stan_list$x_cat2 = c(model.matrix(formula_cat, loss_ratio_new))
  }
  
  ## generate initial values
  init_list = list()
  for (chain in 1:chains){
    sigma=rep(0, stan_list$dimSigma)
    sigma[1] =exp(max(0.1,rnorm(1,0.05,0.05)))-1
    for (s in 2:stan_list$dimSigma){
      sigma[s] = sigma[s-1] + exp(max(0.1,rnorm(1,0.05,0.05)))-1
    }
    
    init_list[[chain]] =
      list(rho=pmax(0.1, rnorm(n = 2,mean =  c(0.8, 0.8),sd = .2)), 
           sigma  = sigma)
    
  }
  

  
  stanMod <- stan_model(file=paste0(stan_dir, stanFile))
  
  
  if (!use_vb){
    stanRun <- sampling(object = stanMod, data = stan_list, open_progress=FALSE, iter=iter, warmup=warmup, chains=chains,init=init_list,
                        control = list(adapt_delta=adapt_delta, max_treedepth =max_treedepth), verbose=verbose)
    summary_stats = summary(stanRun, pars = c("eta", "rho", "sigma", "theta"))$summary
    if (plot_graphs) {
      print(stan_trace(stanRun,pars = c("rho","sigma"))); 
      print(stan_dens(stanRun, pars = "yStar[90]"))
      print(summary_stats)
    }
  } else{
    stanRun <- vb(object = stanMod, data = stan_list, init=init_list)
  }
  
  sigma = rstan::extract(stanRun, pars= "sigma")[[1]]
  sigma = colMeans(sigma)
  
  if (!take_logs) {
    preds = apply(rstan::extract(stanRun, pars= "yStar")[[1]], 2, FUN = scaler_loss_ratio, unscale = T)
  } else {
    preds = apply(rstan::extract(stanRun, pars= "yStar")[[1]], 2, FUN = scaler_loss_ratio_log, unscale = T)
    #preds = exp(preds)
  }
  preds = cbind(loss_ratio_new[, c("AccidentYear","DevelopmentLag")],t(preds))
  #preds[1:5, 1:5]
  preds = as.data.frame(preds) %>%
    arrange(AccidentYear, DevelopmentLag)
  
  preds = preds %>%
    left_join(loss_ratio %>% select(AccidentYear, DevelopmentLag, loss_ratio, loss_ratio_log), by = c("AccidentYear", "DevelopmentLag"))
  
  if (!take_logs){
    loss_ratio_inc_train = preds$loss_ratio
    lr_cumulative = loss_square/ premiums
  } else {
    loss_ratio_inc_train = preds$loss_ratio_log
    lr_cumulative = log(apply((loss_square / premiums),2,function(x) pmax(x,eps)))
  }
  preds$loss_ratio = NULL
  preds$loss_ratio_log = NULL
  
  # Calculate the percentage ranks for incremental loss ratios
  
  lr_increment = lr_cumulative - cbind(0,lr_cumulative[,-ncol(lr_cumulative)]) 
  #lr_cumulative = lr_cumulative[,-1]
  lr_increment = c(t(lr_increment)) #actual incremental loss ratios for full square
  lr_cumulative = c(t(lr_cumulative)) #actual cumulative loss ratios for full square
  rank_val_incremental = vector_rank(v = lr_increment, m = as.matrix(preds[,-1:-2]))
  
  
  # Calculate the cumulative loss ratios
  preds_cum = as.matrix(preds[,-1:-2])

  #replaces the predicted loss ratio with the actual for the training period and calculate cumulative claims
  preds_cum[!is.na(loss_ratio_inc_train), ] = loss_ratio_inc_train[!is.na(loss_ratio_inc_train)] 

  for (i in 1:length(AccidentYear)){
    if (DevelopmentLag[i] >= 2)  {
      preds_cum[i, ] = preds_cum[i, ] + preds_cum[i-1, ]
    } 
  }
  #preds_cum[1:5,1:5]
  # if (!take_logs){
  #   preds_cum = preds_cum +  rep(loss_triangle[,1]/premiums, each= no_dev) # add back initial year claims
  # } else {
  #   preds_cum = preds_cum +  rep(log(pmax(eps,loss_triangle[,1]/premiums)), each= no_dev) # add back initial year claims
  # }  
  #round(rowMeans(preds_cum) - lr_cumulative,4)
  
  # Calculate the percentage ranks for cumulative loss ratios
  rank_val_cumulative = vector_rank(lr_cumulative, preds_cum[,-1:-2])
  #plot(rank_val_cumulative)
  
  
  # Calculation of simulated ultimate losses
  preds_cum_ult = loss_ratio_new %>%
    select(AccidentYear, DevelopmentLag) %>%
    mutate(premiums = rep(premiums, each=no_dev)) 
  
  preds_cum_ult =   cbind(preds_cum_ult, preds_cum * preds_cum_ult$premiums) %>%
    filter(DevelopmentLag == 10)
  sim_ult_losses = colSums(as.matrix(preds_cum_ult[,-1:-3]))
  
  
  
  loss_ratio_new = loss_ratio_new %>%
    mutate(
      predicted_lr = rowMeans(preds_cum) , #colMeans(preds_cum), 
      predicted_lr_med = apply(preds_cum,1, quantile, 0.5) ,
      lower_lr=apply(preds_cum,1,quantile,probs = 0.05), 
      upper_lr=apply(preds_cum,1,quantile,probs = 0.95),
      actual = lr_cumulative,
      predicted_lr_inc = rowMeans(preds[,-1:-2]) , #colMeans(preds_cum), 
      lower_lr_inc=apply(preds[,-1:-2],1,quantile,probs = 0.05), 
      upper_lr_inc=apply(preds[,-1:-2],1,quantile,probs = 0.95),
      actual_inc = lr_increment
    )
  
  if (!take_logs){
    predicted_loss_ratios = loss_ratio_new$predicted_lr
    loss_square_pred = t(matrix(predicted_loss_ratios, no_dev,no_acc))* premiums
    # loss_square_pred_med = cbind(loss_triangle[,1],t(matrix(loss_ratio_new$predicted_lr_med, no_dev,no_acc))* premiums) 
    # loss_square_pred_2 = NA 
  } else {
    predicted_loss_ratios = rowMeans(exp(preds_cum))
    predicted_loss_ratios_med = apply(exp(preds_cum),1,median)
    predicted_loss_ratios_2 = exp(rowMeans(preds_cum))
    loss_square_pred = t(matrix(predicted_loss_ratios, no_dev,no_acc))* premiums
    # loss_square_pred_med = cbind(loss_triangle[,1],t(matrix(predicted_loss_ratios_med, no_dev,no_acc))* premiums) 
    # loss_square_pred_2 = cbind(loss_triangle[,1],t(matrix(predicted_loss_ratios_2, no_dev,no_acc))* premiums)
  }

  if (plot_graphs){
    gplot = ggplot(loss_ratio_new) + 
      geom_vline(aes(xintercept = 1998-AccidentYear), color = "green", linetype="dashed")+
      geom_line(aes(DevelopmentLag, actual,color="actual"))+ 
      geom_line(aes(DevelopmentLag, predicted_lr,color="pred"))+
      geom_point(aes(DevelopmentLag, actual,color="actual"), size=.5)+ 
      geom_point(aes(DevelopmentLag, predicted_lr,color="pred"), size=.5)+
      geom_ribbon(aes(DevelopmentLag, ymin=lower_lr, ymax=upper_lr), color= "grey", alpha=0.2)+
      facet_wrap(~AccidentYear)+ 
      labs(title = "Cumulative Loss Ratio: Actual vs. Predicted") 
    print(gplot)
    
    gplot = ggplot(loss_ratio_new) + 
      geom_vline(aes(xintercept = 1998-AccidentYear), color = "green", linetype="dashed")+
      geom_hline(aes(yintercept = 0), color = "orange", linetype="dashed") +
      geom_line(aes(DevelopmentLag, actual_inc,color="actual")) + 
      geom_line(aes(DevelopmentLag, predicted_lr_inc,color="pred")) +
      geom_point(aes(DevelopmentLag, actual_inc,color="actual"), size=.5) + 
      geom_point(aes(DevelopmentLag, predicted_lr_inc,color="pred"), size=.5) +
      geom_ribbon(aes(DevelopmentLag, ymin=lower_lr_inc, ymax=upper_lr_inc), color= "grey", alpha=0.2)+
      facet_wrap(~AccidentYear)+ 
      scale_x_continuous(breaks =DevelopmentLag_unique)+
      labs(title = "Incremental Loss Ratio: Actual vs. Predicted") 
    print(gplot)
    
    print(stan_plot(stanRun,pars= "sigma"))
    
  }
  if (!calc_ranks){
    return(loss_square = loss_square_pred)
  } else {
    
    return(list(loss_square=loss_square_pred, #predicted loss square
                sim_ult_losses = sim_ult_losses, # simulated ultimate losses
                
                #loss_square_med=loss_square_pred_med, loss_square_2=loss_square_pred_2, 
                #rank_val_incremental=rank_val_incremental,rank_val_cumulative= rank_val_cumulative, rank_ultimate= rank_ultimate, 
                #sigma = sigma,
                summary_stats = summary_stats,
                ae = sum(loss_square[,10])/mean(sim_ult_losses)))
  }
}



#' Calculate the predicted loss ratio and supporting analytics using an incremental loss ratio model
#' 
#' Uses a vector for mu prior
#' 
#' @param loss_triangle the Loss triangle and loss square
#' @param premiums Vector of premiums
#' @param calc_ranks If TRUE then additional statistics are calculated
#' @param formula_x formula for square exponential 
#' @param formula_H formula for linear terms in kernel (not in mean function)
#' @param mu The vector for the prior mu values for the GP
#' @param prior_eta, prior_rho, prior_theta, prior_sigma, prior_alpha, prior_beta Priors for various parameters in the Stan model.  If a model does not require one of these priors, then such prior is ignored.
#' @param take_logs If TRUE, then response is logarithm is used as the predictor.  Note, if TRUE, then LHS of formula_x must be y_log_scaled, otherwise use y_scaled
#' @param plot_graphs If TRUE, then sample investigative graphs are plotted.
#' @param virtual_data if TRUE, then a training column for duration 11 is added with incremental loss = 0%
#' @return If calc_ranks = FALSE, then only the predicted loss_square is returned.  If TRUE, then the percentage ranks, sigma estimates are also returned.
stan_incremental_loss_ratio_virtual3 = function(loss_triangle, premiums=NULL, loss_square = NULL, calc_ranks = FALSE, 
                                               formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                               formula_H = y_scaled ~  AccidentYear_scaled + DevelopmentLag_scaled -1 , formula_cat = NULL,
                                               eps = 0.001, plot_graphs = FALSE,
                                               chains = 6,iter = 600, warmup =300, adapt_delta = 0.9,max_treedepth = 10,
                                               prior_eta = 1, prior_rho = c(5.566019, 5.570317), prior_theta = 1, prior_sigma = 1, 
                                               prior_beta = 1, take_logs = FALSE, 
                                               scaled_lower_bound = NULL, l_bound = 0, mu=9:1/9-.5,
                                               virtual_data = FALSE,virtual_point = NA, 
                                               stanFile = "gp_compound_mean_03.stan", use_vb = FALSE, verbose = FALSE){

  if(is.null(loss_square)) 
    loss_square = matrix(NA, nrow(loss_triangle), ncol(loss_triangle))
  
  AccidentYear_unique = AccidentYear = as.numeric(rownames(loss_triangle))
  DevelopmentLag_unique = DevelopmentLag = as.numeric(colnames(loss_triangle))
  #DevelopmentLag = DevelopmentLag[-1]
  no_acc = length(AccidentYear); 
  no_dev = length(DevelopmentLag)
  
  AccidentYear = rep(AccidentYear, each = no_dev)
  DevelopmentLag = rep(DevelopmentLag, times = no_acc)
  
  if (virtual_data){
    no_dev_train = no_dev + 1
  } else {
    no_dev_train = no_dev
  }
  
  loss_ratio_triangle = (loss_triangle / premiums)
  #loss_ratio_triangle = loss_ratio_triangle[,-1] - loss_ratio_triangle[,-ncol(loss_ratio_triangle)]
  loss_ratio_triangle = loss_ratio_triangle - cbind(0,loss_ratio_triangle[,-ncol(loss_ratio_triangle)])
  loss_ratio_triangle_log = log(apply((loss_triangle / premiums),2,function(x) pmax(x,eps)))
  loss_ratio_triangle_log = loss_ratio_triangle_log - cbind(0,loss_ratio_triangle_log[,-ncol(loss_ratio_triangle_log)])  
  
  if (virtual_data){
    loss_ratio_triangle = cbind(loss_ratio_triangle, virtual_point)
    loss_ratio_triangle_log = cbind(loss_ratio_triangle_log, virtual_point)
    colnames(loss_ratio_triangle) = 2:(no_dev_train+1)
    colnames(loss_ratio_triangle_log) = 2:(no_dev_train+1)
  }

  
  loss_ratio = as_tibble(loss_ratio_triangle) %>%
    mutate(AccidentYear = as.integer(rownames(loss_ratio_triangle))) %>%
    gather(DevelopmentLag, loss_ratio, -AccidentYear) %>%
    filter(!is.na(loss_ratio)) %>%
    mutate(
      DevelopmentLag = as.integer(DevelopmentLag)
    )
  loss_ratio_log = as_tibble(loss_ratio_triangle_log) %>%
    mutate(AccidentYear = as.integer(rownames(loss_ratio_triangle_log))) %>%
    gather(DevelopmentLag, loss_ratio_log, -AccidentYear) %>%
    filter(!is.na(loss_ratio_log)) %>%
    mutate(
      DevelopmentLag = as.integer(DevelopmentLag)
    )
  loss_ratio = loss_ratio %>%
    left_join(loss_ratio_log, by = c("AccidentYear", "DevelopmentLag"))
  
  if (virtual_data){
    actual_data_flag = loss_ratio$DevelopmentLag <= no_dev_train
  } else {
    actual_data_flag = rep(TRUE,nrow(loss_ratio))
  }

  
  scaler_accident_year = scaler(loss_ratio$AccidentYear)
  scaler_development_lag = scaler(loss_ratio$DevelopmentLag)
  DevelopmentLag_transform = function(x) log(x-min(loss_ratio$DevelopmentLag)+1)
  scaler_loss_ratio = scaler(loss_ratio$loss_ratio[actual_data_flag], l_bound = scaled_lower_bound)
  scaler_loss_ratio_log = scaler(loss_ratio$loss_ratio_log[actual_data_flag], l_bound = scaled_lower_bound)
  #y_sd = sd(loss_ratio$loss_ratio)
  
  if (!take_logs) {
    l_bound = scaler_loss_ratio(l_bound)[1,1]
    mu = c(scaler_loss_ratio(mu))
  } else {
    l_bound = scaler_loss_ratio_log(l_bound)[1,1]
    mu = c(scaler_loss_ratio_log(mu))
  }
  loss_ratio = loss_ratio %>%
    mutate(AccidentYear_scaled = scaler_accident_year(AccidentYear)[,1],
           DevelopmentLag_scaled = scaler_development_lag(DevelopmentLag)[,1], 
           DevelopmentLag_trans = DevelopmentLag_transform(DevelopmentLag),
           y_scaled = scaler_loss_ratio(loss_ratio)[,1],
           y_log_scaled = scaler_loss_ratio_log(loss_ratio_log)[,1]) 
  
  
  # Create inputs for out-of-sample predictions
  loss_ratio_new = tibble(AccidentYear ,DevelopmentLag) %>%
    mutate(AccidentYear_scaled = scaler_accident_year(AccidentYear)[,1],
           DevelopmentLag_scaled = scaler_development_lag(DevelopmentLag)[,1], 
           DevelopmentLag_trans = DevelopmentLag_transform(DevelopmentLag),
           y_scaled = 0,
           y_log_scaled = 0)
  
  formula_x = update(formula_x, . ~ . -1) # remove intercept term
  formula_H = update(formula_H, . ~ . -1) # remove intercept term
  
  stan_list <- list(
    N1 = sum(actual_data_flag),
    N1_virtual = length(actual_data_flag),
    N2 = nrow(loss_ratio_new),
    x1 =  model.matrix(formula_x, loss_ratio[actual_data_flag,]),
    x1_virtual =  model.matrix(formula_x, loss_ratio),
    y1 =  c(loss_ratio[,all.vars(update(formula_x, .~0)) ][[1]])[actual_data_flag],
    y1_virtual = c(loss_ratio[,all.vars(update(formula_x, .~0)) ][[1]]),
    y1 = pmax(l_bound, c(loss_ratio[,all.vars(update(formula_x, .~0)) ][[1]]))[actual_data_flag],
    y1_virtual = pmax(l_bound, c(loss_ratio[,all.vars(update(formula_x, .~0)) ][[1]])),
    dev_lag1 = loss_ratio$DevelopmentLag[actual_data_flag],
    dev_lag1_virtual = loss_ratio$DevelopmentLag,
    x2 = model.matrix(formula_x, loss_ratio_new),
    dev_lag2 = loss_ratio_new$DevelopmentLag,
    mu = mu,
    l_bound=l_bound,
    prior_eta=prior_eta,
    prior_rho = prior_rho,
    prior_theta = prior_theta,
    prior_sigma = prior_sigma,
    prior_beta = prior_beta,
    dimSigma = 11-2 + 1 *(virtual_data == 1) 
  )
  
  if (take_logs & (all.vars(formula_x)[1] != "y_log_scaled")){
    warning("if logs are used, LHS of formula_x must be 'y_log_scaled'")
  } 
  
  stan_list$D = ncol(stan_list$x1)
  
  if(!is.null(formula_H)){
    stan_list$H1 = t(model.matrix(formula_H, loss_ratio[actual_data_flag,]))
    stan_list$H1_virtual = t(model.matrix(formula_H, loss_ratio))
    stan_list$H2 = t(model.matrix(formula_H, loss_ratio_new))
    stan_list$H1 = stan_list$H1 - apply(stan_list$H1,1,min)
    stan_list$H1_virtual = stan_list$H1_virtual - apply(stan_list$H1,1,min)
    stan_list$H2 = stan_list$H2 - apply(stan_list$H1,1,min)
    stan_list$dimH <- nrow(stan_list$H1)
  }
  
  
  if(!is.null(formula_cat)){
    formula_cat = update(formula_cat,  ~ . -1) # remove intercept term
    stan_list$x_cat1 = c(model.matrix(formula_cat, loss_ratio))
    stan_list$x_cat2 = c(model.matrix(formula_cat, loss_ratio_new))
  }
  
  ## generate initial values
  init_list = list()
  for (chain in 1:chains){
    sigma=rep(0, stan_list$dimSigma)
    sigma[1] =exp(max(0.1,rnorm(1,0.05,0.05)))-1
    for (s in 2:stan_list$dimSigma){
      sigma[s] = sigma[s-1] + exp(max(0.1,rnorm(1,0.05,0.05)))-1
    }
    
    init_list[[chain]] =
      list(rho=pmax(0.1, rnorm(n = 2,mean =  c(0.8, 0.8),sd = .2)), 
           sigma  = sigma)
    
  }
  
  
  stanMod <- stan_model(file=paste0(stan_dir, stanFile))
  
  
  if (!use_vb){
    stanRun <- sampling(object = stanMod, data = stan_list, iter=iter, warmup=warmup, chains=chains,init=init_list,
                        control = list(adapt_delta=adapt_delta, max_treedepth =max_treedepth) , verbose= verbose)
    summary_stats = summary(stanRun, pars = c("eta", "rho", "sigma"))$summary
    if (plot_graphs) {
      print(stan_trace(stanRun,pars = c("rho","sigma"))); 
      print(stan_dens(stanRun, pars = "yStar[90]"))
      print(summary_stats) 
    }
  } else{
    stanRun <- vb(object = stanMod, data = stan_list, init=init_list)
  }
  
  sigma = rstan::extract(stanRun, pars= "sigma")[[1]]
  sigma = colMeans(sigma)
  
  if (!take_logs) {
    preds = apply(rstan::extract(stanRun, pars= "yStar")[[1]], 2, FUN = scaler_loss_ratio, unscale = T)
  } else {
    preds = apply(rstan::extract(stanRun, pars= "yStar")[[1]], 2, FUN = scaler_loss_ratio_log, unscale = T)
    #preds = exp(preds)
  }
  preds = cbind(loss_ratio_new[, c("AccidentYear","DevelopmentLag")],t(preds))
  #preds[1:5, 1:5]
  preds = as.data.frame(preds) %>%
    arrange(AccidentYear, DevelopmentLag)
  
  preds = preds %>%
    left_join(loss_ratio %>% select(AccidentYear, DevelopmentLag, loss_ratio, loss_ratio_log), by = c("AccidentYear", "DevelopmentLag"))
  
  if (!take_logs){
    loss_ratio_inc_train = preds$loss_ratio
    lr_cumulative = loss_square/ premiums
  } else {
    loss_ratio_inc_train = preds$loss_ratio_log
    lr_cumulative = log(apply((loss_square / premiums),2,function(x) pmax(x,eps)))
  }
  preds$loss_ratio = NULL
  preds$loss_ratio_log = NULL
  
  # Calculate the percentage ranks for incremental loss ratios
  
  lr_increment = lr_cumulative[,-1] - lr_cumulative[,-ncol(lr_cumulative)]
  lr_cumulative = lr_cumulative[,-1]
  lr_increment = c(t(lr_increment)) #actual incremental loss ratios for full square
  lr_cumulative = c(t(lr_cumulative)) #actual cumulative loss ratios for full square
  #lr_increment - loss_ratio_inc_train
  rank_val_incremental = vector_rank(v = lr_increment, m = as.matrix(preds[,-1:-2]))
  
  # Calculate the cumulative loss ratios
  preds_cum = as.matrix(preds[,-1:-2])
  #dim(preds_cum)
  
  #replaces the predicted loss ratio with the actual for the training period and calculate cumulative claims
  preds_cum[!is.na(loss_ratio_inc_train), ] = loss_ratio_inc_train[!is.na(loss_ratio_inc_train)] 
  #preds[1:5, 1:5]
  #rowMeans(preds_cum) - loss_ratio_inc_train
  
  for (i in 1:length(AccidentYear)){
    if (DevelopmentLag[i] > 2)  {
      preds_cum[i, ] = preds_cum[i, ] + preds_cum[i-1, ]
    } 
  }
  preds_cum[1:5,1:5]
  if (!take_logs){
    preds_cum = preds_cum +  rep(loss_triangle[,1]/premiums, each= no_dev) # add back initial year claims
  } else {
    preds_cum = preds_cum +  rep(log(pmax(eps,loss_triangle[,1]/premiums)), each= no_dev) # add back initial year claims
  }  
  round(rowMeans(preds_cum) - lr_cumulative,4)
  
  # Calculate the percentage ranks for cumulative loss ratios
  rank_val_cumulative = vector_rank(lr_cumulative, preds_cum[,-1:-2])
  #plot(rank_val_cumulative)
  
  # Calculation of simulated ultimate losses
  preds_cum_ult = loss_ratio_new %>%
    select(AccidentYear, DevelopmentLag) %>%
    mutate(premiums = rep(premiums, each=no_dev)) %>%
    cbind(preds_cum * premiums) %>%
    filter(DevelopmentLag == 10)
  sim_ult_losses = colSums(as.matrix(preds_cum_ult[,-1:-3]))
  
  
  loss_ratio_new = loss_ratio_new %>%
    mutate(
      predicted_lr = rowMeans(preds_cum) , #colMeans(preds_cum), 
      predicted_lr_med = apply(preds_cum,1, quantile, 0.5) ,
      lower_lr=apply(preds_cum,1,quantile,probs = 0.05), 
      upper_lr=apply(preds_cum,1,quantile,probs = 0.95),
      actual = lr_cumulative,
      predicted_lr_inc = rowMeans(preds[,-1:-2]) , #colMeans(preds_cum), 
      lower_lr_inc=apply(preds[,-1:-2],1,quantile,probs = 0.05), 
      upper_lr_inc=apply(preds[,-1:-2],1,quantile,probs = 0.95),
      actual_inc = lr_increment
    )
  
  if (!take_logs){
    predicted_loss_ratios = loss_ratio_new$predicted_lr
    loss_square_pred = cbind(loss_triangle[,1],t(matrix(predicted_loss_ratios, no_dev,no_acc))* premiums) 
    loss_square_pred_med = cbind(loss_triangle[,1],t(matrix(loss_ratio_new$predicted_lr_med, no_dev,no_acc))* premiums) 
    loss_square_pred_2 = NA 
  } else {
    predicted_loss_ratios = rowMeans(exp(preds_cum))
    predicted_loss_ratios_med = apply(exp(preds_cum),1,median)
    predicted_loss_ratios_2 = exp(rowMeans(preds_cum))
    loss_square_pred = cbind(loss_triangle[,1],t(matrix(predicted_loss_ratios, no_dev,no_acc))* premiums) 
    loss_square_pred_med = cbind(loss_triangle[,1],t(matrix(predicted_loss_ratios_med, no_dev,no_acc))* premiums) 
    loss_square_pred_2 = cbind(loss_triangle[,1],t(matrix(predicted_loss_ratios_2, no_dev,no_acc))* premiums)
  }
  
  if (plot_graphs){
    gplot = ggplot(loss_ratio_new) + 
      geom_vline(aes(xintercept = 1998-AccidentYear), color = "green", linetype="dashed")+
      geom_line(aes(DevelopmentLag, actual,color="actual"))+ 
      geom_line(aes(DevelopmentLag, predicted_lr,color="pred"))+
      geom_point(aes(DevelopmentLag, actual,color="actual"), size=.5)+ 
      geom_point(aes(DevelopmentLag, predicted_lr,color="pred"), size=.5)+
      geom_ribbon(aes(DevelopmentLag, ymin=lower_lr, ymax=upper_lr), color= "grey", alpha=0.2)+
      facet_wrap(~AccidentYear)+ 
      labs(title = "Cumulative Loss Ratio: Actual vs. Predicted") 
    print(gplot)
    
    gplot = ggplot(loss_ratio_new) + 
      geom_vline(aes(xintercept = 1998-AccidentYear), color = "green", linetype="dashed")+
      geom_hline(aes(yintercept = 0), color = "orange", linetype="dashed") +
      geom_line(aes(DevelopmentLag, actual_inc,color="actual")) + 
      geom_line(aes(DevelopmentLag, predicted_lr_inc,color="pred")) +
      geom_point(aes(DevelopmentLag, actual_inc,color="actual"), size=.5) + 
      geom_point(aes(DevelopmentLag, predicted_lr_inc,color="pred"), size=.5) +
      geom_ribbon(aes(DevelopmentLag, ymin=lower_lr_inc, ymax=upper_lr_inc), color= "grey", alpha=0.2)+
      facet_wrap(~AccidentYear)+ 
      labs(title = "Incremental Loss Ratio: Actual vs. Predicted") 
    print(gplot)
    
    print(stan_plot(stanRun,pars= "sigma"))
    
  }
  if (!calc_ranks){
    return(loss_square = loss_square_pred)
  } else {
    
    errors = error_calcs(premiums = premiums, predicted = loss_square_pred, actual = loss_square)
    
    errors$rmse_mean = rmse_triangle(actual_sq = loss_square/premiums, predicted_sq = loss_square_pred/    premiums,actual_tr = loss_triangle/premiums)
    errors$rmse_med =  rmse_triangle(actual_sq = loss_square/premiums, predicted_sq = loss_square_pred_med/premiums,actual_tr = loss_triangle/premiums)
    
    
    return(list(loss_square=loss_square_pred,
                sim_ult_losses=sim_ult_losses, 
                summary_stats = summary_stats))
  }
}


#' Uses GP with Stan to project Loss Development Factors in order to estimate developed losses
#' @param loss_triangle the Loss triangle and loss square
#' @param premiums Vector of premiums.  Not used for this algorithm
#' @param formula_x formula for square exponential 
#' @param formula_H formula for linear terms in  in mean function
#' @param prior_eta, prior_rho, prior_theta, prior_sigma Priors for various parameters in the Stan model.  If a model does not require one of these priors, then such prior is ignored.
#' @param plot_graphs If TRUE, then sample investigative graphs are plotted.
#' @param stanFile  The stan file used to run the analysis
#' @return If calc_ranks = FALSE, then only the predicted loss_square is returned.  If TRUE, then the percentage ranks, sigma estimates are also returned.

stan_ldf3 = function(loss_triangle, premiums = NULL,
                       formula_x = dev_log ~ AccidentYear_scale + DevelopmentLag_scale -1,
                       formula_H = ~ AccidentYear + log(DevelopmentLag), 
                       eps = 0.001, plot_graphs = FALSE,
                       chains = 6,iter = 600, warmup = iter/4, adapt_delta = 0.95,max_treedepth = 10,
                       prior_eta = 1, prior_rho = c(5.5,5.5), prior_beta = 1, prior_sigma = 1,
                       stanFile = "gpFastPredict04.stan"){
  
  AccidentYear = as.numeric(rownames(loss_triangle))
  DevelopmentLag = as.numeric(colnames(loss_triangle))[-1]
  no_acc = length(AccidentYear); no_dev = length(DevelopmentLag)
  AccidentYear = rep(AccidentYear, each = no_dev)
  DevelopmentLag = rep(DevelopmentLag, times = no_acc)
  
  AccidentYear_mean = mean(AccidentYear); AccidentYear_sd = sd(AccidentYear)
  AccidentYear_scale = (AccidentYear-AccidentYear_mean) / AccidentYear_sd
  DevelopmentLag_mean = mean(DevelopmentLag); DevelopmentLag_sd = sd(DevelopmentLag)
  DevelopmentLag_scale = (DevelopmentLag-DevelopmentLag_mean) / DevelopmentLag_sd
  
  
  dev = dev_factors(loss_triangle) %>%
    as_tibble() %>%
    mutate(AccidentYear = as.integer(rownames(loss_triangle))) %>%
    gather(DevelopmentLag, dev, -AccidentYear) %>%
    filter(!is.na(dev)) %>%
    mutate(dev_log = log(pmax(dev,1)-1+eps),
           DevelopmentLag = as.integer(DevelopmentLag))
  
  AccidentYear_mean = mean(dev$AccidentYear)
  AccidentYear_sd = sd(dev$AccidentYear)
  DevelopmentLag_mean = mean(dev$DevelopmentLag)
  DevelopmentLag_sd = sd(dev$DevelopmentLag)
  
  dev = dev %>%
    mutate(AccidentYear_scale = (AccidentYear - AccidentYear_mean)/AccidentYear_sd,
           DevelopmentLag_scale = (DevelopmentLag - DevelopmentLag_mean)/DevelopmentLag_sd)
  
  mean_dev_log = mean(dev$dev_log)
  sd_dev_log = sd(dev$dev_log)
  dev$dev_log = (dev$dev_log - mean_dev_log)/sd_dev_log
  
  # Create inputs for out-of-sample predictions
  dev_new = tibble(AccidentYear ,DevelopmentLag) %>%
    mutate(AccidentYear_scale = (AccidentYear - AccidentYear_mean)/AccidentYear_sd,
           DevelopmentLag_scale = (DevelopmentLag - DevelopmentLag_mean)/DevelopmentLag_sd)
  
  ########################################################3
  
  # Establish data to be passed to Stan
  N1 = nrow(dev)
  N2 = nrow(dev_new)
  
  dev_new$blank = 0
  
  stan_list <- list(
    N1 = N1,
    N2 = N2,
    x1 = model.matrix(formula_x, dev),
    H1 = model.matrix(formula_H, dev), 
    y1 = c(dev[,all.vars(update(formula_x, .~0)) ][[1]]),
    x2 = model.matrix(update(formula_x, blank ~ .), dev_new),
    H2 = model.matrix(formula_H, dev_new), 
    prior_eta=prior_eta,
    prior_rho=prior_rho,
    prior_sigma=prior_sigma,
    prior_beta=prior_beta
  )
  stan_list$D = ncol(stan_list$x1)
  stan_list$dimH <- ncol(stan_list$H1)
  
  if (ncol(stan_list$H1) != ncol(stan_list$H2)){
    stop("dev and dev_new have different dimensions for mean function ")
  }
  browser()
  
  stanMod <- stan_model(file=paste0(stan_dir, stanFile))
  
  if (!optimizeMod){
    model <- sampling(object = stanMod, data = stan_list, iter=iter, chains=chains, control = list(adapt_delta=adapt_delta, max_treedepth =max_treedepth) )
  } else {
    model <- optimizing(object = stanMod, data = stan_list, verbose=T)
  }
  
  
  unscale = function(x){
    exp(x* sd_dev_log + mean_dev_log) +1- eps  
  }
  
  sims = rstan::extract(model, pars= "yStar")[[1]]
  
  sims= apply(sims, 2, unscale)
  
  dev_preds =colMeans(sims)
  
  preds = tibble(AccidentYear=AccidentYear, DevelopmentLag=DevelopmentLag)
  
  
  ########################################################
  
  # Calc simulated ultimate losses
  sim_ult_losses = rep(NA, iter)
  for (i in 1:iter){
    sim_ldf = matrix(sims[i,], nrow = N, byrow = TRUE)
    sim_loss_square = project_losses(loss_triangle,dev_factors = sim_ldf)
    sim_ult_losses[i] = sum(sim_loss_square[, ncol(sim_loss_square)])
  }
  ldf_sq = matrix(dev_preds, nrow = N, byrow = TRUE)
  loss_square = project_losses(loss_triangle,dev_factors = ldf_sq)
  
  
  # Plot Actual vs. Predicted Factors
  if (plot_graphs){
    print(ggplot(cbind(preds, intercept=1)) + 
            theme_bw() +
            geom_point(aes(DevelopmentLag, dev, color="Actual"))+
            geom_line(aes(DevelopmentLag, dev_preds, color="Predicted")) + 
            facet_wrap(~p0("Acc. Yr ",AccidentYear)) + 
            geom_hline(aes(yintercept = intercept), size = 0.25, linetype = 2) + 
            labs(title = p0("GP Fit to Log Loss Development Factors\nformula = ", format(formula) )))
    
    print(ggplot(cbind(preds, intercept=1)) + 
            theme_bw() +
            geom_point(aes(AccidentYear, dev, color="Actual"))+
            geom_line(aes(AccidentYear, dev_preds, color="Predicted")) + 
            facet_wrap(~p0("Dev. Yr ",sprintf("%02d", DevelopmentLag))) + 
            geom_hline(aes(yintercept = intercept), size = 0.25, linetype = 2) + 
            labs(title = p0("GP Fit to Log Loss Development Factors\nformula = ", format(formula) )))
  }
  
  
  list(loss_square =loss_square, sim_ult_losses = sim_ult_losses) 
}


if (F){
  
}


#' Performs all analytics for a business line for a given algorithm
#' 
#' @param data The desired loss reserve database. This is loaded with the load_losses_subset() function.
#' @param projection_fun Function to run loss development algorithm
#' @param business_line Business line to be analyzed
#' @param include_rank_tests Set to TRUE for analyses that projection_fun returns a distribution of losses in the variable sim_ult_losses
#' @param print_res If TRUE, then interim results are printed
#' @param rank_size Number of cells in square that need to be predicted for rank val tests
#' @results The various performance analytics that are used in the report
business_line_calcs = function(data, projection_fun=stan_incremental_loss_ratio2, business_line="prodliab", include_rank_tests = TRUE,
                               print_res = TRUE, square_size =10, rank_size = 90, ...){

  codes = unlist(unique(data %>% filter(line==business_line) %>% select(GRCODE)))
  N_tests = length(codes)
  Q = square_size # Ultimate column
  N = square_size # number of accident years
  

  results = tibble(code = codes, 
                   claims_ultimate_actual = NA, # 
                   claims_ultimate_predicted = NA, # 
                   mean_loss_ratio_actual  = NA, # 
                   mean_loss_ratio_predicted  = NA, # 
                   mean_loss_ratio_weighted_actual  = NA, # 
                   mean_loss_ratio_weighted_predicted  = NA, #
                   coverage_90 = NA,
                   crps1 = NA,
                   #crps2 = NA,
                   nlpd = NA,
                   rank_ultimate = NA, 
                   #ultimate_claim_percentile = NA, 
                   ks_rejection_level = NA,
                   ks_D = NA,
                   ks_reject_hypothesis  = NA
  ) 
  mat_sa_claims_actual = NULL
  mat_sa_claims_predicted = NULL
  mat_sa_loss_ratio_actual = NULL
  mat_sa_loss_ratio_predicted = NULL
  summary_stats = NULL
  
  #colnames(rmse_step_ahead)=1:9
  # rank_vals_inc = matrix(NA, N_tests, rank_size)
  # rank_vals_cum = matrix(NA, N_tests, rank_size)
  i=1

  for  (code in codes){
    print(paste("Working on... ", i, " of ", N_tests, " ", business_line, code))
    loss = loss_table(data, business_line = business_line, company_code = code)
    premiums = get_premium(data, business_line = business_line, company_code = code)
    

    preds = try(projection_fun(loss_triangle = loss$paid_tr, premiums=premiums,loss_square=loss$paid_sq, ...))

    if (is.error(preds)){
      results[i,] = NA
      step_ahead_claims_actual = rep(NA, N-1)
      step_ahead_claims_predicted = rep(NA, N-1)
      step_ahead_loss_ratio_actual = rep(NA, N-1)
      step_ahead_loss_ratio_predicted = rep(NA, N-1)
      warning(paste("ERROR in Calculation", results[i,]))
    } else {
      
      if(!is.list(preds)) 
        preds = list(loss_square=preds)
      #####  Calculation of performance metrics
      results$claims_ultimate_actual[i] = sum(loss$paid_sq[,Q])
      results$claims_ultimate_predicted[i] = sum(preds$loss_square[,Q])
      results$mean_loss_ratio_actual[i] = mean(loss$paid_sq[,Q]/premiums)
      results$mean_loss_ratio_predicted[i] = mean(preds$loss_square[,Q]/premiums)
      results$mean_loss_ratio_weighted_actual[i] = sum(loss$paid_sq[,Q])/sum(premiums)
      results$mean_loss_ratio_weighted_predicted[i] = sum(preds$loss_square[,Q])/sum(premiums)
      
      if (include_rank_tests){
        #ultimate_losses = sum(preds$loss_square[,ncol(preds$loss_square)])
        ultimate_losses = sum(loss$paid_sq[,Q])
        results$coverage_90[i] = coverage(r =  ultimate_losses, r_distribution = preds$sim_ult_losses, percentage = 0.9)
        results$crps1[i] = crps_sample(y = ultimate_losses, dat = preds$sim_ult_losses)
        #results$crps2[i] = crps_fun(r =  ultimate_losses, r_distribution = preds$sim_ult_losses)
        results$nlpd[i] = nlpd_fun(r =  ultimate_losses, r_distribution  = preds$sim_ult_losses)
        results$rank_ultimate[i] = vector_rank(v = ultimate_losses, m = matrix(preds$sim_ult_losses,nrow = 1))
      }
      
      # Step Ahead Statistics
      step_ahead_claims_actual = rep(0, N-1)
      step_ahead_claims_predicted = rep(0, N-1)
      step_ahead_loss_ratio_actual = rep(0, N-1)
      step_ahead_loss_ratio_predicted = rep(0, N-1)
      
      diag_indicator <- row(preds$loss_square) + col(preds$loss_square) - N - 1

      diag_vec_actual = split(loss$paid_sq,diag_indicator)
      diag_vec_predicted = split(preds$loss_square,diag_indicator)
      prem_vec = split(row(preds$loss_square), diag_indicator)
      prem_vec = map(prem_vec, ~ premiums[.x])

      
      for (n in 1:(N-1)){
        d= as.character(n)
        step_ahead_claims_actual[n] = sum(diag_vec_actual[[d]])
        step_ahead_claims_predicted[n] = sum(diag_vec_predicted[[d]])
        step_ahead_loss_ratio_actual[n] = sum(diag_vec_actual[[d]]) / sum(prem_vec[[d]]) / length(prem_vec[[d]])
        step_ahead_loss_ratio_predicted[n] = sum(diag_vec_predicted[[d]]) / sum(prem_vec[[d]]) / length(prem_vec[[d]])
      }
      
      # rank_vals_inc[i,] = preds$rank_val_incremental
      # rank_vals_cum[i,] = preds$rank_val_cumulative
      
    }
    names(step_ahead_claims_actual) = 1:(N-1)
    names(step_ahead_claims_predicted) = 1:(N-1)
    names(step_ahead_loss_ratio_actual) = 1:(N-1)
    names(step_ahead_loss_ratio_predicted) = 1:(N-1)
    
    # Matrices of step ahead actual and projected metrics
    mat_sa_claims_actual = rbind(mat_sa_claims_actual, step_ahead_claims_actual)
    mat_sa_claims_predicted = rbind(mat_sa_claims_predicted, step_ahead_claims_predicted)
    mat_sa_loss_ratio_actual = rbind(mat_sa_loss_ratio_actual, step_ahead_loss_ratio_actual)
    mat_sa_loss_ratio_predicted = rbind(mat_sa_loss_ratio_predicted, step_ahead_loss_ratio_predicted)
    
    if (!is.error(preds)){
      if (!is.null(preds$summary_stats) ){
        summary_stats = rbind(summary_stats,
                              cbind(tibble(business_line = business_line, code = code,par = rownames(preds$summary_stats)), as_tibble(preds$summary_stats)))
      } else {
        summary_stats = NA
      }
    } else {
      summary_stats = NA
    }
    
    i=i+1
  }
  metrics = list()
  metrics$rmse_claims = rmse(results$claims_ultimate_actual, results$claims_ultimate_predicted, na.rm = T )
  metrics$rmse_loss_ratio = rmse(results$mean_loss_ratio_actual, results$mean_loss_ratio_predicted, na.rm = T )
  metrics$rmse_loss_ratio_weighted = rmse(results$mean_loss_ratio_weighted_actual, results$mean_loss_ratio_weighted_predicted, na.rm = T )
  metrics$rmse_step_ahead_claims = rep(NA, N-1)
  metrics$rmse_step_ahead_loss_ratio = rep(NA, N-1)
  metrics$coverage_90 = mean(results$coverage_90,  na.rm = T)
  metrics$total_crps1 = sum(results$crps1,  na.rm = T)
  metrics$total_nlpd = sum(results$nlpd,  na.rm = T)
  metrics$valid_calc_count = sum(!is.na(results$claims_ultimate_predicted))
  
  for (i in 1:N-1){
    metrics$rmse_step_ahead_claims[i] = rmse(mat_sa_claims_actual[,i], mat_sa_claims_predicted[,i], na.rm = T)
    metrics$rmse_step_ahead_loss_ratio[i] = rmse(mat_sa_loss_ratio_actual[,i], mat_sa_loss_ratio_predicted[,i], na.rm = T)
  }

  if (include_rank_tests){
    rank_test = ks_percentile(p_i = results$rank_ultimate, title=business_line)
    metrics$ks_rejection_level = rank_test$rejection_level
    metrics$ks_D = rank_test$D
    metrics$ks_reject_hypothesis  = rank_test$reject_hypothesis
  } else {
    rank_test = NA
    metrics$ks_rejection_level = NA
    metrics$ks_D = NA
    metrics$ks_reject_hypothesis  = NA
  }

  list(results=results,
       metrics = metrics,
       rank_test = rank_test, 
       summary_stats = summary_stats
       )
}





#' Performs all analytics for all business line for a given algorithm
#' 
#' @param data The desired loss reserve database. This is loaded with the load_losses_subset() function.
#' @param projection_fun Function to run loss development algorithm
#' @param business_lines Business line to be analyzed
#' @param combined_co If TRUE, then the algorithm runs all companies within business line a single loop. Default is false
#' @param print_res If TRUE, then interim results are printed
#' @param rank_size Number of cells in square that need to be predicted for rank val tests
#' @results The various performance analytics that are used in the report
all_business_line_calcs = function(business_lines = c("comauto", "medmal", "othliab","ppauto","prodliab", "wkcomp"), combined_co = FALSE,
                                   data, projection_fun=stan_incremental_loss_ratio2, print_res = TRUE,
                                   include_rank_tests = TRUE,
                                   square_size =10, rank_size = 90, ...){
  results =list()
  

  
  for (bus in business_lines){
    print(paste("Running Business Line:", bus))
    if (!combined_co){
      results[[bus]] =  business_line_calcs(data = data, projection_fun=projection_fun, business_line=bus, print_res = print_res, 
                                          include_rank_tests = include_rank_tests, square_size =square_size, rank_size = rank_size, ...)
    } else {
      results[[bus]] =  projection_fun(data=data, business_line=bus, ...)
    }
  }
  results[[bus]]$business_lines = business_lines
  
  results_all = list()
  results_all$all = cbind(metric = c("RMSE Claims","RMSE LR", "RMSE Weighted LR", "Coverage 90","CRPS1", "NLPD", "KS Rejection Level","KS D", "KS Result", "Valid Count"),   
                   results %>% map_df(~ round(c(.x$metrics$rmse_claims, .x$metrics$rmse_loss_ratio, .x$metrics$rmse_loss_ratio_weighted,
                                                .x$metrics$coverage_90, .x$metrics$total_crps1, .x$metrics$total_nlpd,
                                                .x$metrics$ks_rejection_level, .x$metrics$ks_D, .x$metrics$ks_reject_hypothesis, .x$metrics$valid_calc_count
                                                ),4)))
  results_all$all_step_ahead_claims = cbind(Step_Ahead = 1:length(results[[1]]$metrics$rmse_step_ahead_claims),
                                     results %>% map_df(~ round(c(.x$metrics$rmse_step_ahead_claims),4))) 
  results_all$all_step_ahead_loss_ratio = cbind(Step_Ahead = 1:length(results[[1]]$metrics$rmse_step_ahead_claims),
                                         results %>% map_df(~ round(c(.x$metrics$rmse_step_ahead_loss_ratio),4))) 
  
  results[["all"]] = results_all
  results
}

#' Print key results from function all_business_line_calcs
#' 
#' @param res The output of all_business_line_calcs
print_loss_results = function(res){
  print(res$all$all)
  print(res$all$all_step_ahead_claims)
  print(res$all$all_step_ahead_loss_ratio)
  
  if (! "summary_stats" %in% names(res[[1]])){
    stan_sim = FALSE
  } else {
    if (is.null(dim(res[[1]]$summary_stats))) {
      stan_sim = FALSE
    } else {
      stan_sim = TRUE
    }
    
  }
  
  for (i in 1:length(res)){
    if (names(res)[i] != "all"){
      #print(paste0("Business Line: ", names(res_gp_inc)[i]))
      print(res[[i]]$rank_test$pp_plot)
      print(res[[i]]$rank_test$hist_plot)
      
      if(stan_sim){
        if (is.matrix(res[[i]]$summary_stats)){
          par = rownames(res[[i]]$summary_stats)
          res[[i]]$summary_stats = as_tibble(res[[i]]$summary_stats) %>%
            mutate(par = par,
                   business_line = names(res)[i])
        }
      }
    }
  }
  
  if (stan_sim){
    pars = res
    pars$all = NULL
    
    pars = do.call(rbind, pars %>% map( ~ .x$summary_stats))
    # if (is.null(pars$par)){
    #   pars$par = rownames(pars)
    # }
    
    print(pars)
    
    pars$index = as.numeric(str_extract(pars$par, pattern = "(?<=\\[).+?(?=\\])")) 
    pars$par = ifelse(is.na(str_extract(pars$par, pattern = "^.*(?=(\\[))")) , 
                      pars$par, 
                      str_extract(pars$par, pattern = "^.*(?=(\\[))"))
    
    if(length(res)-1 == nrow(pars %>% filter(par == "eta"))){
      #multi company models
      gg_par_eta = ggplot(pars %>% filter(par == "eta")) +
        theme_minimal()+
        geom_col(aes(business_line, mean,fill=business_line)) + 
        labs(title = "Mean Value of eta by Business Line")
      print(gg_par_eta)
      gg_par_rho_ay = ggplot(pars %>% filter(par == "rho",index ==1)) +
        theme_minimal()+
        geom_col(aes(business_line, mean,fill=business_line)) + 
        labs(title = "Mean Value of rho_AY by Business Line")
      print(gg_par_rho_ay)
      gg_par_rho_dl = ggplot(pars %>% filter(par == "rho",index ==2)) +
        theme_minimal()+
        geom_col(aes(business_line, mean,fill=business_line)) + 
        labs(title = "Mean Value of rho_DL by Business Line")
      print(gg_par_rho_dl)
      gg_par_theta = ggplot(pars %>% filter(par == "theta")) +
        theme_minimal()+
        geom_col(aes(business_line, mean,fill=business_line)) + 
        facet_wrap(~index) +
        labs(title = "Mean Value of theta by Business Line")
      print(gg_par_theta)
      
      
    } else {
      # individual company models
      gg_par_eta = ggplot(pars %>% filter(par == "eta")) +
        theme_minimal()+
        geom_density(aes(mean,fill = business_line), alpha=.4) + 
        guides(fill=guide_legend(title="Business Line"))+
        labs(title = "Distribution of eta by Business Line")
      print(gg_par_eta)
      gg_par_rho_ay = ggplot(pars %>% filter(par == "rho",index ==1)) +
        theme_minimal()+
        geom_density(aes(mean,fill = business_line), alpha=.4) + 
        guides(fill=guide_legend(title="Business Line"))+
        labs(title = "Distribution of rho_AY by Business Line")
      print(gg_par_rho_ay)
      gg_par_rho_dl = ggplot(pars %>% filter(par == "rho",index ==2)) +
        theme_minimal()+
        geom_density(aes(mean,fill = business_line), alpha=.4) + 
        guides(fill=guide_legend(title="Business Line"))+
        labs(title = "Distribution of rho_DL by Business Line")
      print(gg_par_rho_dl)
      gg_par_theta = ggplot(pars %>% filter(par == "theta")) +
        theme_minimal()+
        geom_density(aes(mean,fill = business_line), alpha=.4) +
        facet_wrap(~index) +
        guides(fill=guide_legend(title="Business Line"))+
        labs(title = "Distribution of theta by Business Line")
      print(gg_par_theta)
    }
    
    gg_par_sigma = ggplot(pars %>% 
                            filter(par == "sigma") %>%
                            group_by(index, business_line) %>%
                            summarize(mean = mean(mean)) %>%
                            ungroup() %>%
                            mutate(development_lag = max(index)-index+1)) +
      theme_minimal()+
      geom_line(aes(development_lag,mean , color = business_line)) + 
      guides(fill=guide_legend(title="Business Line"))+
      labs(title = "Distribution of sigma by Business Line", x="Development Lag", y="sigma")
    print(gg_par_sigma)
    gg_rhat = ggplot(pars ) +
      theme_minimal()+
      geom_density(aes(Rhat,fill = business_line), alpha=.4) + 
      guides(fill=guide_legend(title="Business Line"))+
      labs(title = "Distribution of Rhat by Business Line")
    print(gg_rhat)
  }
  
} 



incr_hurdle = function(){
  
  #####################################################
  ## Test Single Case
  
  bus = "comauto"
  codes = unlist(unique(data %>% filter(line==bus) %>% select(GRCODE)))
  company_code = 353
  loss = loss_table(data, business_line = bus, company_code = company_code)
  premiums = get_premium(data, business_line = bus, company_code = company_code)
  
  ##### 
  #Basic  Setup
  adapt_delta = 0.9; iter = 400; warmup=300; max_treedepth=10; chains=6
  a=Sys.time()
  ans_hurdle = stan_incremental_loss_ratio2(loss_triangle = loss$paid_tr,  loss_square = loss$paid_sq, 
                                            premiums = premiums, calc_ranks = T,
                                            prior_sigma = .5, 
                                            mu=9:1 / 9 -.5,
                                            l_bound = 0,
                                            prior_rho = c(9.53, 9.66),
                                            #prior_rho = c(5.5, 5.5), #prior_rho = c(9.4,8.4), #prior_rho = c(26.2,19.7), #prior_rho = c(14.9, 11.2), #prior_rho = c(8.8, 11.8),
                                            plot_graphs = T,
                                            iter = iter, warmup=warmup, adapt_delta = adapt_delta, max_treedepth =max_treedepth, chains=chains, 
                                            formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                            formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                            stanFile = "gp_hurdle_02_phiapprox2.stan"
  )
  b=Sys.time(); b-a
  ans_hurdle
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_hurdle$loss_square/premiums,actual_tr = loss$paid_tr/premiums)
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_hurdle$loss_square_med/premiums,actual_tr = loss$paid_tr/premiums)
  

  ans_hurdle$summary_stats
  ans_hurdle$errors
  plot(ans_hurdle$errors$rmse_step_ahead)
  
  
  # Business line analysis
  ans_bus = business_line_calcs(data = data, projection_fun =  stan_incremental_loss_ratio2, business_line = "medmal", take_logs=F, 
                         prior_sigma = .5, 
                         mu=9:1 / 9 -.5,
                         prior_rho = c(9.4,8.4), 
                         iter = iter, warmup=warmup, chains= chains,
                         calc_ranks = T, plot_graphs = F,
                         formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag -1,
                         formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                         stanFile = "gp_hurdle_02_phiapprox2.stan"
  )
  save(ans_bus, file = paste0(output_dir, "report_medmal_results.RData"))
  
}


