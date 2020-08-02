###  Analysis of Reserves using Multi-Company Loss Triangles in the same model

# Load and prepare data
source('report_fun.R')
#source('comp_gp_01.R')

 
#' Multiple Company Analysis: using an incremental loss ratio model
#' 
#' Uses a vector for mu prior
#' 
#' @param data the full loss development database
#' @param business_line The business line to be analyzed
#' @param calc_ranks If TRUE then additional statistics are calculated
#' @param formula_x formula for square exponential 
#' @param formula_H formula for linear terms in kernel (not in mean function)
#' @param mu The vector for the prior mu values for the GP
#' @param prior_eta, prior_rho, prior_theta, prior_sigma, prior_alpha, prior_beta Priors for various parameters in the Stan model.  If a model does not require one of these priors, then such prior is ignored.
#' @param take_logs If TRUE, then response is logarithm is used as the predictor.  Note, if TRUE, then LHS of formula_x must be y_log_scaled, otherwise use y_scaled
#' @param plot_graphs If TRUE, then sample investigative graphs are plotted.
#' @param virtual_data if TRUE, then a training column for duration 11 is added with incremental loss = 0%
#' @return If calc_ranks = FALSE, then only the predicted loss_square is returned.  If TRUE, then the percentage ranks, sigma estimates are also returned.
stan_ilr_multi_co = function(data, business_line,
                             formula_x = loss_ratio_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled, 
                             formula_H =  loss_ratio_scaled ~ AccidentYear_scaled + DevelopmentLagLog_scaled, 
                             #formula_H =  loss_ratio_scaled ~ AccidentYear_scaled + log(DevelopmentLag_scaled - min(DevelopmentLag_scaled)+1), 
                             formula_cat =  ~ GRCODE,
                             eps = 0.001, plot_graphs = FALSE,
                             chains = 6,iter = 600, warmup =300, adapt_delta = 0.9,max_treedepth = 10,
                             prior_eta = 1, prior_rho = c(5.566019, 5.570317), prior_theta = 1, prior_sigma = 1, 
                             prior_alpha = 1, prior_rho_company = 1, prior_beta = 1, take_logs = FALSE, 
                             mu=10:1/10-.5,
                             
                             scaled_lower_bound = NULL, l_bound = 0, 
                             
                             virtual_data = FALSE,
                             stanFile = "gp_compound_mean_03.stan", 
                             use_vb = FALSE){

  
  # Selecte business line
  block_test = data %>%
    filter(line == business_line)
  
  # Sort block_test and calculate loss development factor
  block_test = block_test %>%
    arrange(GRCODE, AccidentYear, DevelopmentLag) %>%
    group_by(GRCODE, AccidentYear) %>%
    mutate(
      #ldf = CumPaidLoss / lag(CumPaidLoss),
      lag_CumPaidLoss = replace_na(lag(CumPaidLoss),0),
      loss_ratio = (CumPaidLoss - lag_CumPaidLoss) / EarnedPremNet,
      DevelopmentLagLog = log(DevelopmentLag)) 
  
  # %>%
  #   filter(is.finite(ldf),
  #          is.finite(loss_ratio))
  # 
  
  # Make scalers
  scaler_AccidentYear = scaler(x_reference = block_test$AccidentYear)
  scaler_DevelopmentLag = scaler(x_reference = block_test$DevelopmentLag)
  scaler_DevelopmentLagLog = scaler(x_reference = block_test$DevelopmentLagLog)
  
  #scaler_ldf = scaler(x_reference = block_test$ldf[block_test$DevelopmentYear <= 1997])
  
  scaler_loss_ratio = scaler(x_reference = block_test$loss_ratio[block_test$DevelopmentYear <= 1997]) # Scaler is based on upper triangle only
  l_bound = scaler_loss_ratio(l_bound)
  mu = c(scaler_loss_ratio(mu))
  
  
  # Create new training block and test blocks with scaled data
  
  ## Test Data
  block_test = block_test %>% 
    mutate(AccidentYear_scaled = scaler_AccidentYear(AccidentYear),
           DevelopmentLag_scaled = scaler_DevelopmentLag(DevelopmentLag),
           DevelopmentLagLog_scaled = scaler_DevelopmentLagLog(DevelopmentLagLog),
           #ldf_scaled = scaler_ldf(ldf),
           loss_ratio_scaled = scaler_loss_ratio(loss_ratio)
           #ldf_log = log(pmax(ldf,1) -1 +eps),
    ) 
  
  # Training Data
  block_train = block_test %>% 
    filter(DevelopmentYear <= 1997)
  
  
  # mu_vec_train = mu[block_train$DevelopmentLag-1]
  # mu_vec_train = scaler_loss_ratio(mu_vec_train)
  # mu_vec_test = mu[block_test$DevelopmentLag-1]
  # mu_vec_test = scaler_loss_ratio(mu_vec_test)

  # if (virtual_data){
  #   no_dev_train = no_dev + 1
  # } else {
  #   no_dev_train = no_dev
  # }
  

  # if (virtual_data){
  #   loss_ratio_triangle = cbind(loss_ratio_triangle, 0)
  #   loss_ratio_triangle_log = cbind(loss_ratio_triangle_log, 0)
  #   colnames(loss_ratio_triangle) = 2:(no_dev_train+1)
  #   colnames(loss_ratio_triangle_log) = 2:(no_dev_train+1)
  # }
  
  
  
  formula_x = update(formula_x, . ~ . -1) # remove intercept term
  formula_H = update(formula_H,  ~ . -1) # remove intercept term
  
  stan_list <- list(
    N1 = nrow(block_train),
    N2 = nrow(block_test),
    x1 =  model.matrix(formula_x, block_train),
    y1 = pmax(l_bound, c(block_train[,all.vars(update(formula_x, .~0)) ][[1]])),
    dev_lag1 = block_train$DevelopmentLag,
    x2 = model.matrix(formula_x, block_test),
    dev_lag2 = block_test$DevelopmentLag,
    mu = mu,
    l_bound=l_bound[[1]],
    prior_eta=prior_eta,
    prior_rho = prior_rho,
    prior_theta = prior_theta,
    prior_sigma = prior_sigma,
    prior_alpha = prior_alpha,
    prior_rho_company=prior_rho_company,
    prior_beta = prior_beta, #not used
    dimSigma = 11-1 + 1 *(virtual_data == 1)
  )
  
  if (take_logs & (all.vars(formula_x)[1] != "y_log_scaled")){
    warning("if logs are used, LHS of formula_x must be 'y_log_scaled'")
  } 
  
  stan_list$D = ncol(stan_list$x1)
  
  if(!is.null(formula_H)){
    stan_list$H1 = t(model.matrix(formula_H, block_train))
    stan_list$H2 = t(model.matrix(formula_H, block_test))
    #stan_list$H1 = stan_list$H1 - apply(stan_list$H1,1,min)  #forces linear trend terms to be positive - this is done in Stan
    #stan_list$H2 = stan_list$H2 - apply(stan_list$H1,1,min) #forces linear trend terms to be positive - this is done in Stan
    stan_list$dimH <- nrow(stan_list$H1)
  }
  
  
  if(!is.null(formula_cat)){
    formula_cat = update(formula_cat,  ~ . -1) # remove intercept term
    stan_list$x_cat1 = c(model.matrix(formula_cat, block_train))
    stan_list$x_cat2 = c(model.matrix(formula_cat, block_test))
    
    # converts company into consecutive integers starting at 1
    stan_list$x_cat1 =  factor(stan_list$x_cat1) 
    stan_list$x_cat2 =  factor(stan_list$x_cat2, levels = levels(stan_list$x_cat1))
    
    if (any(is.na(stan_list$x_cat2))){
      stop("there are companies in the test set that are not in the training set")
    }
    stan_list$x_cat1 = as.numeric(stan_list$x_cat1)
    stan_list$x_cat2 = as.numeric(stan_list$x_cat2)
    
    stan_list$noCompanies = length(unique(stan_list$x_cat1))
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
      list(rho=pmax(0.1, rnorm(n = ncol(stan_list$x1),mean =  c(0.8, 0.8),sd = .2)), 
           eta=pmax(0.1, rnorm(n = 1,mean =  c(0.8),sd = .2)), 
           alpha=pmax(0.1, rnorm(n = 1,mean =  c(0.2),sd = .2)),
           rho_company=pmax(0.1, rnorm(n = stan_list$noCompanies,mean =  c(0.8, 0.8),sd = .2)), 
           sigma  = sigma)
    
  }
  

  stanMod <- stan_model(file=paste0(stan_dir, stanFile))
  
  
  if (!use_vb){
    stanRun <- sampling(object = stanMod, data = stan_list, iter=iter, warmup=warmup, chains=chains,init=init_list,
                        control = list(adapt_delta=adapt_delta, max_treedepth =max_treedepth) )
    save(stanRun, file = paste0(output_dir, bus, "_multi_co_stan.Rdata"))
    
    if ("alpha" %in% stanRun@sim$pars_oi) {
      summary_stats = summary(stanRun, pars = c("eta", "rho", "theta", "alpha", "sigma"))$summary
      alpha_par = TRUE
    } else {
      summary_stats = summary(stanRun, pars = c("eta", "rho", "theta", "rho_company", "sigma"))$summary
      alpha_par = FALSE
    }
    if (plot_graphs) {
      if (alpha_par) {
        print(stan_trace(stanRun,pars = c("rho","sigma", "alpha")))
      } else {
        print(stan_trace(stanRun,pars = c("rho","sigma", "rho_company")))
      }

      print(stan_dens(stanRun, pars = "yStar[90]"))
      print(summary_stats)
    }
  } else{
    stanRun <- vb(object = stanMod, data = stan_list, init=init_list)
  }
  
  
  sigma = rstan::extract(stanRun, pars= "sigma")[[1]]
  sigma = colMeans(sigma)
  
  ###########################################################################################
  
  preds_scaled = t(rstan::extract(stanRun, pars= "yStar")[[1]])
  str(preds_scaled)
  preds = apply(preds_scaled, 2, FUN = scaler_loss_ratio, unscale = T)
  
  block_test$preds = rowMeans(preds)
  block_test$lower_lr = apply(preds,1, quantile, probs=0.05)
  block_test$upper_lr = apply(preds,1, quantile, probs=0.95)
  
  ggplot(block_test %>% filter(AccidentYear == 1990)) + 
    geom_line(aes(DevelopmentLag,loss_ratio, color="actual"))+
    geom_line(aes(DevelopmentLag,preds, color="predicted"))+
    facet_wrap(~GRCODE, scales= "free")
  ggplot(block_test %>% filter(AccidentYear == 1997)) + 
    geom_line(aes(DevelopmentLag,loss_ratio, color="actual"))+
    geom_line(aes(DevelopmentLag,preds, color="predicted"))+
    facet_wrap(~GRCODE, scales= "free")
  
  company = block_test$GRCODE[1]
  ggplot(block_test %>% filter(GRCODE == company)) + 
    geom_line(aes(DevelopmentLag,loss_ratio, color="actual"))+
    geom_line(aes(DevelopmentLag,preds, color="predicted"))+
    facet_wrap(~AccidentYear, scales= "free")+
    geom_ribbon(aes(DevelopmentLag, ymin=lower_lr, ymax=upper_lr), color= "grey", alpha=0.2)+
    geom_vline(aes(xintercept = 1998-AccidentYear), color = "green", linetype="dashed")+
    labs(title=paste("ILR for Company", company), x = "Development Lag", y="Incremental Loss Ratio")
  
  # Calculate RMSE for each company
  companies = unique(block_test$GRCODE)
  #draws_combined = do.call(rbind, model$draws)
  
  no_sims = ncol(preds)
  

  # results2 = data.frame(line=bus, GRCODE = companies, 
  #                      rmse = rep(NA, length(companies)), 
  #                      rmse_lr = rep(NA, length(companies)))
  
  
  results = tibble(code = companies, 
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
  
  ultimate_losses_sim = rep(NA, no_sims)
  
  mat_sa_claims_actual = NULL
  mat_sa_claims_predicted = NULL
  mat_sa_loss_ratio_actual = NULL
  mat_sa_loss_ratio_predicted = NULL

  
  
  

  loss_ratios_pred = block_test %>%
    dplyr::select(GRCODE, AccidentYear, DevelopmentLag, preds) %>%
    spread(DevelopmentLag, preds)
  N_accident_year = length(unique(block_test$AccidentYear))
  N_development_lag = length(unique(block_test$DevelopmentLag))
  preds_losses=matrix(NA, nrow=nrow(preds), ncol = ncol(preds))
  
  
  
  
  
  i=1
  for (company in companies){
    print(company)
    losses_tr = loss_table(data, business_line = business_line, company_code = company)$paid_tr
    losses_sq = loss_table(data, business_line = business_line, company_code = company)$paid_sq 
    
    premiums = get_premium(data, business_line = business_line, company_code = company)
    
    Q = ncol(losses_sq)
    
    for (sim in 1:no_sims){
      selected_rows = block_test$GRCODE == company
      loss_ratio_mat_sim = matrix(preds[selected_rows,sim], nrow = N_accident_year, ncol = N_development_lag, byrow = T)
      #loss_ratio_mat_sim = cbind("1" = losses_sq[,1]/premiums, loss_ratios_sim)
      
      loss_ratio_mat_sim  = t(apply(loss_ratio_mat_sim,1,cumsum)) # convert incremental to cumulative losses
      
      losses_projected_sim = loss_ratio_mat_sim * premiums
      ultimate_losses_sim[sim] = sum(losses_projected_sim[, ncol(losses_projected_sim)])
      
      preds_losses[selected_rows, sim] = c(t(losses_projected_sim))
    }
    
    loss_ratios_pred_mat = cbind( as.matrix(loss_ratios_pred %>% filter(GRCODE == company))[,-c(1,2)])
    
    loss_ratios_pred_mat  = t(apply(loss_ratios_pred_mat,1,cumsum))
    
    #losses_tr[,-1] / losses_tr[,-10]
    losses_projected = loss_ratios_pred_mat * premiums # mean projected loss ratio

    
    
    ##################
    
    results$claims_ultimate_actual[i] = sum(losses_sq[,Q])
    results$claims_ultimate_predicted[i] = sum(losses_projected[,Q])
    results$mean_loss_ratio_actual[i] = mean(losses_sq[,Q]/premiums)
    results$mean_loss_ratio_predicted[i] = mean(losses_projected[,Q]/premiums)
    results$mean_loss_ratio_weighted_actual[i] = sum(losses_sq[,Q])/sum(premiums)
    results$mean_loss_ratio_weighted_predicted[i] = sum(losses_projected[,Q])/sum(premiums)
    
    
    ultimate_losses = sum(losses_sq[,Q])
    results$coverage_90[i] = coverage(r =  ultimate_losses, r_distribution = ultimate_losses_sim, percentage = 0.9)
    results$crps1[i] = crps_sample(y = ultimate_losses, dat = ultimate_losses_sim)
    #results$crps2[i] = crps_fun(r =  ultimate_losses, r_distribution = ultimate_losses_sim)
    results$nlpd[i] = nlpd_fun(r =  ultimate_losses, r_distribution  = ultimate_losses_sim)
    results$rank_ultimate[i] = vector_rank(v = ultimate_losses, m = matrix(ultimate_losses_sim,nrow = 1))
    
    
    
    #################
    
    
    # results_rmse = rmse_triangle(actual_sq = losses_sq, 
    #                              predicted_sq = losses_projected, 
    #                              actual_tr = losses_tr)
    # results_rmse_lr = rmse_triangle(actual_sq = losses_sq/premiums, 
    #                                 predicted_sq = losses_projected/premiums, 
    #                                 actual_tr = losses_tr)
    # results2[i,c("rmse", "rmse_lr")] = c(results_rmse, results_rmse_lr)
    # 
    # 
    
    
    
    i=i+1
  }
  results
  #results2
  #mean(results2$rmse_lr)
  
  
  max_accident_year = max(block_test$AccidentYear)
  block_test = block_test %>%
    mutate(step_ahead = DevelopmentYear - max_accident_year)

  block_test$quantile = vector_rank(block_test$CumPaidLoss, preds_losses)
  block_test$preds_loss = rowMeans(preds_losses)
  
  # colnames(preds_losses) = paste0("sim_", 1:ncol(preds_losses))
  # preds_losses_df = as.data.frame(preds_losses)
  # preds_losses_df = cbind(as.data.frame(block_test[, c("GRCODE", "GRNAME", "AccidentYear", "DevelopmentYear", "DevelopmentLag","EarnedPremNet", "CumPaidLoss",  "step_ahead")]), preds_losses_df)
  # str(preds_losses_df)
  
  step_ahead = block_test %>%
    filter(step_ahead>0) %>%
    group_by(GRCODE, step_ahead) %>%
    summarize(premiums = sum(EarnedPremNet),
              claims_actual = sum(CumPaidLoss),
              claims_predicted = sum(preds_loss)) %>%
    select(GRCODE, step_ahead, premiums, claims_actual, claims_predicted) %>%
    mutate(loss_ratio_actual = claims_actual/premiums,
           loss_ratio_predicted = claims_predicted/premiums)
  step_ahead
  step_ahead_summary =   step_ahead %>%
    ungroup() %>%
    #select(-GRCODE) %>%
    group_by(step_ahead) %>%
    summarize(rmse_step_ahead_claims = rmse(claims_actual, claims_predicted, na.rm = T),
           rmse_step_ahead_loss_ratio = rmse(loss_ratio_actual, loss_ratio_predicted, na.rm = T))
  step_ahead_summary
  
  
  metrics = list()
  metrics$rmse_claims = rmse(results$claims_ultimate_actual, results$claims_ultimate_predicted, na.rm = T )
  metrics$rmse_loss_ratio = rmse(results$mean_loss_ratio_actual, results$mean_loss_ratio_predicted, na.rm = T )
  metrics$rmse_loss_ratio_weighted = rmse(results$mean_loss_ratio_weighted_actual, results$mean_loss_ratio_weighted_predicted, na.rm = T )
  metrics$rmse_step_ahead_claims = step_ahead_summary$rmse_step_ahead_claims
  metrics$rmse_step_ahead_loss_ratio = step_ahead_summary$rmse_step_ahead_loss_ratio
  metrics$coverage_90 = mean(results$coverage_90,  na.rm = T)
  metrics$total_crps1 = sum(results$crps1,  na.rm = T)
  metrics$total_nlpd = sum(results$nlpd,  na.rm = T)
  metrics$valid_calc_count = sum(!is.na(results$claims_ultimate_predicted))
  

  rank_test = ks_percentile(p_i = results$rank_ultimate, title=business_line)
  metrics$ks_rejection_level = rank_test$rejection_level
  metrics$ks_D = rank_test$D
  metrics$ks_reject_hypothesis  = rank_test$reject_hypothesis

  company = block_test$GRCODE[1]
  ggplot(block_test %>% filter(GRCODE == company)) + 
    geom_line(aes(DevelopmentLag,CumPaidLoss, color="actual"))+
    geom_line(aes(DevelopmentLag,preds_loss, color="predicted"))+
    facet_wrap(~AccidentYear, scales= "free")+
    #geom_ribbon(aes(DevelopmentLag, ymin=lower_lr, ymax=upper_lr), color= "grey", alpha=0.2)+
    geom_vline(aes(xintercept = 1998-AccidentYear), color = "green", linetype="dashed")+
    labs(title=paste("Cumulative Paid Losses for Company", company), x = "Development Lag", y="Incremental Loss Ratio")
  
  return(list(results=results, metrics = metrics, rank_test = rank_test, summary_stats = summary_stats, stanRun=stanRun, block_test = block_test, preds_losses=preds_losses))
  
}  
  
  
if (F){  
  
  ############################################################################################
  if (!take_logs) {
    preds = apply(rstan::extract(stanRun, pars= "yStar")[[1]], 2, FUN = scaler_loss_ratio, unscale = T)
  } else {
    preds = apply(rstan::extract(stanRun, pars= "yStar")[[1]], 2, FUN = scaler_loss_ratio_log, unscale = T)
    #preds = exp(preds)
  }
  preds = cbind(block_test[, c("GRCODE", "AccidentYear","DevelopmentLag")],t(preds))
  #preds[1:5, 1:5]
  preds = as.data.frame(preds) %>%
    arrange(GRCODE, AccidentYear, DevelopmentLag)
  
  preds = preds %>%
    left_join(block_train %>% select(AccidentYear, DevelopmentLag, loss_ratio, loss_ratio_log), by = c("GRCODE", "AccidentYear", "DevelopmentLag"))
  
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
  
  
  
  block_test = block_test %>%
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
    predicted_loss_ratios = block_test$predicted_lr
    loss_square_pred = cbind(loss_triangle[,1],t(matrix(predicted_loss_ratios, no_dev,no_acc))* premiums) 
    loss_square_pred_med = cbind(loss_triangle[,1],t(matrix(block_test$predicted_lr_med, no_dev,no_acc))* premiums) 
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
    gplot = ggplot(block_test) + 
      geom_vline(aes(xintercept = 1998-AccidentYear), color = "green", linetype="dashed")+
      geom_line(aes(DevelopmentLag, actual,color="actual"))+ 
      geom_line(aes(DevelopmentLag, predicted_lr,color="pred"))+
      geom_point(aes(DevelopmentLag, actual,color="actual"), size=.5)+ 
      geom_point(aes(DevelopmentLag, predicted_lr,color="pred"), size=.5)+
      geom_ribbon(aes(DevelopmentLag, ymin=lower_lr, ymax=upper_lr), color= "grey", alpha=0.2)+
      facet_wrap(~AccidentYear)+ 
      labs(title = "Cumulative Loss Ratio: Actual vs. Predicted") 
    print(gplot)
    
    gplot = ggplot(block_test) + 
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
    
    return(list(loss_square=loss_square_pred,loss_square_med=loss_square_pred_med, loss_square_2=loss_square_pred_2, 
                rank_val_incremental=rank_val_incremental,rank_val_cumulative= rank_val_cumulative, sigma = sigma,
                errors = errors, summary_stats = summary_stats))
  }
}

if (F){
  iter=400; warmup=300; chains=5
  
  set.seed(777)
  company_max = 5
  #company_max = Inf
  for (bus in business_lines){
    companies = unique(data$GRCODE[data$line == bus])
    company_count = min(length(companies), company_max)
    all_lines = sample(companies, size = company_count, replace = FALSE) 
    data = data %>%
      filter(line != bus | GRCODE %in% sample(all_lines, size = company_max))
  }
  
  
  bus = "wkcomp"
  stanFile = "gp_hurdle_multi_co_02.stan"
  #stanFile = "gp_hurdle_multi_co_test2.stan"
  stanMod <- stan_model(file=paste0(stan_dir, stanFile))
  (prior_rho = rho_prior_tune(.01, 0.34, 3.1)$par)
  
  ans = stan_ilr_multi_co(data = data, business_line = bus, iter = iter, warmup=warmup, chains=chains, stanFile = stanFile, prior_rho = prior_rho, prior_sigma=0.5)
  ans$results
  ans$metrics
  stan_trace(ans$stanRun)
  
  iter=400; warmup=300; chains=5
  business_lines = c("comauto", "medmal", "othliab","ppauto","prodliab", "wkcomp")
  business_lines = c("wkcomp")
  ans = all_business_line_calcs(projection_fun = stan_ilr_multi_co, business_lines = business_lines, 
                                multi_co = TRUE, data = data,
                                iter = iter, warmup=warmup, chains=chains, stanFile = "gp_hurdle_multi_co_02.stan",
                                prior_eta = 1, prior_rho=prior_rho, prior_theta = 1, prior_sigma = 0.5, 
                                prior_alpha = 1)

  ans
  
ans$wkcomp$summary_stats
}