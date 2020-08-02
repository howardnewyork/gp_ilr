####  Comparison of RMSE and Percentage Ranks Across Different GP Methods

##########################################################################
####  Setup

# Load Functions
source('loss_reserve.R')

#Load Data
data = load_losses_subset()

# Parameters

iter = 200
warmup=100
chains =6
adapt_delta=0.9
max_treedepth =10
business_lines = c("comauto", "medmal", "othliab","ppauto","prodliab", "wkcomp")
#business_lines = c( "medmal","prodliab", "wkcomp")
#business_lines = c( "medmal")

#limit 20 cases per line
if (F){
  table(data[, c("line", "GRCODE")])
  company_max = 20
  for (bus in business_lines){
    all_lines = unique(data$GRCODE[data$line == bus])
    data = data %>%
      filter(line != bus | GRCODE %in% all_lines[1:company_max])
  }
  length(table(data$GRCODE[data$line == "comauto"]))
}


##########################################################################
####  Functions



##########################################################################
####  Incremental Loss Ratio Method with NO trunction

incremental_no_trunc = function(){
  

  ## Test Single Case
  
  bus = "ppauto"
  codes = unlist(unique(data %>% filter(line==bus) %>% select(GRCODE)))
  company_code = 1538
  loss = loss_table(data, business_line = bus, company_code = company_code)
  premiums = get_premium(data, business_line = bus, company_code = company_code)
  
  ans_mean_het = stan_incremental_loss_ratio(loss_triangle = loss$paid_tr, loss_square = loss$paid_sq, premiums = premiums, prior_sigma = 1, plot_graphs = T,
                                             iter = iter, warmup=warmup, calc_ranks = T,
                                             formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                             formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                             stanFile = "gp_compound_mean_het_03.stan"
  )
  
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_mean_het$loss_square/premiums,    actual_tr = loss$paid_tr/premiums)
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_mean_het$loss_square_med/premiums,actual_tr = loss$paid_tr/premiums)
  colMeans(loss$paid_sq[,-1]/premiums - loss$paid_sq[,-10]/premiums <=0)

  ans_warp = stan_incremental_loss_ratio(loss_triangle = loss$paid_tr,  loss_square = loss$paid_sq, premiums = premiums, prior_sigma = 1, plot_graphs = T,
                                             iter = iter, warmup=warmup,calc_ranks = T,
                                             formula_x = y_scaled ~ AccidentYear_scaled + beta_warp(DevelopmentLag, 1,5) -1,
                                             formula_H = y_scaled ~ AccidentYear_scaled + beta_warp(DevelopmentLag, 1,5) -1,
                                             stanFile = "gp_compound_mean_het_03.stan"
  )
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_warp$loss_square/premiums,actual_tr = loss$paid_tr/premiums)
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_warp$loss_square_med/premiums,actual_tr = loss$paid_tr/premiums)
  
  
  ## Analyze All Business Lines
  
  results =list()
  for (bus in business_lines){
    results[[bus]] =  business_line_analysis(projection_fun =  stan_incremental_loss_ratio, business_line = bus, take_logs=F, 
                                             prior_sigma = .5, iter = iter, warmup=warmup, chains= chains,
                                             calc_ranks = T, plot_graphs = F,
                                             formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                             formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                             stanFile = "gp_compound_mean_het_03.stan"
    )
  }
  #save(results, file = paste0(output_dir, "results_incremenal.RData"))
  #load(file = paste0(output_dir, "results_incremenal.RData"))
  # RMSE Vals
  map(results, ~ median(.$rmse_loss_ratio, na.rm = T))
  map(results, ~ mean(.$rmse_loss_ratio, na.rm = T))
  
  # Percentage Rank Plots
  plot_ranks(results)
  
}  


##########################################################################
####  Incremental Loss Ratio Method WITH TRUNCATION

incr_with_trun = function(){
  
  
  ## Test Single Case
  
  bus = "medmal"
  codes = unlist(unique(data %>% filter(line==bus) %>% select(GRCODE)))
  company_code = 669
  loss = loss_table(data, business_line = bus, company_code = company_code)
  premiums = get_premium(data, business_line = bus, company_code = company_code)
  
  # Centered vs. Non-Centered methods
  #Method 1
  adapt_delta = 0.99
  ans_mean_trunc = stan_incremental_loss_ratio(loss_triangle = loss$paid_tr,  loss_square = loss$paid_sq, premiums = premiums, 
                                               prior_sigma = .05, 
                                               prior_rho = c(26.2,19.7), #prior_rho = c(14.9, 11.2),
                                               plot_graphs = T,
                                               iter = iter, warmup=warmup, adapt_delta = adapt_delta, max_treedepth =max_treedepth, 
                                               formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                               formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                               stanFile = "gp_comp_het_trunc_04.stan"
                                               #stanFile = "gp_comp_het_trunc_03.stan"
  )
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_mean_trunc/premiums,actual_tr = loss$paid_tr/premiums)
  
  #Method 2
  ans_mean_trunc = stan_incremental_loss_ratio(loss_triangle = loss$paid_tr,  loss_square = loss$paid_sq, premiums = premiums, prior_sigma = 1, plot_graphs = T,
                                               iter = 800, warmup=700, adapt_delta = 0.99, max_treedepth =10, 
                                               formula_x = y_scaled ~ AccidentYear_scaled + beta_warp(DevelopmentLag, 1,10) -1,
                                               formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                               stanFile = "gp_comp_het_trunc_04.stan"
  )
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_mean_trunc/premiums,actual_tr = loss$paid_tr/premiums)
  
 
  ## Analyze All Business Lines
  
  results =list()
  for (bus in business_lines){
    results[[bus]] =  business_line_analysis(projection_fun =  stan_incremental_loss_ratio, business_line = bus, take_logs=F, 
                                             prior_sigma = .05, prior_rho = c(14.9, 11.2),
                                             iter = iter, warmup=warmup, chains= chains,
                                             calc_ranks = T, plot_graphs = F,
                                             formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag -1,
                                             formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                             stanFile = "gp_compound_mean_het_03.stan"
    )
  }
  #save(results, file = paste0(output_dir, "results_incremenal.RData"))
  #load(file = paste0(output_dir, "results_incremenal.RData"))
  # RMSE Vals
  map(results, ~ median(.$results$rmse_loss_ratio, na.rm = T))
  map(results, ~ mean(.x$results$rmse_loss_ratio, na.rm = T))
  
  # Percentage Rank Plots
  plot_ranks(results)
  
}  
  
##########################################################################
####  Incremental Loss Ratio Method HURDLE MODEL

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
  adapt_delta = 0.9; iter = 600; warmup=500; max_treedepth=10; chains=6
  a=Sys.time()
  ans_hurdle = stan_incremental_loss_ratio(loss_triangle = loss$paid_tr,  loss_square = loss$paid_sq, premiums = premiums, calc_ranks = T,
                                           prior_sigma = .5, mu=-0.75,
                                           l_bound = 0,
                                           #prior_rho = c(5.5, 5.5), 
                                           #prior_rho = c(9.4,8.4), 
                                           prior_rho = c(9.53, 9.66),
                                           #prior_rho = c(26.2,19.7), 
                                           #prior_rho = c(14.9, 11.2),
                                           plot_graphs = T,
                                           iter = iter, warmup=warmup, adapt_delta = adapt_delta, max_treedepth =max_treedepth, 
                                           formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                           formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                           stanFile = "gp_hurdle_01.stan"
                                           #stanFile = "gp_hurdle_zmethod_01.stan"
  )
  b=Sys.time(); b-a
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_hurdle$loss_square/premiums,actual_tr = loss$paid_tr/premiums)
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_hurdle$loss_square_med/premiums,actual_tr = loss$paid_tr/premiums)
  
  ##### 
  # Vector mu
  #adapt_delta = 0.95; iter = 800; warmup=700; max_treedepth=15; chains=10
  adapt_delta = 0.9; iter = 400*2; warmup=300*2; max_treedepth=12; chains=6
  adapt_delta = 0.9; iter = 400; warmup=300; max_treedepth=10; chains=6
  #adapt_delta = 0.7; iter = 800; warmup=700; max_treedepth=15; chains=8
  stanMod <- stan_model(file=paste0(stan_dir, stanFile))
  
  a=Sys.time()
  ans_hurdle = stan_incremental_loss_ratio2(loss_triangle = loss$paid_tr,  loss_square = loss$paid_sq, 
                                            premiums = premiums, calc_ranks = T,
                                           prior_sigma = .5, 
                                           #mu = rep(0,9),
                                           #mu=9:1 / 5 -1.3,
                                           #mu=9:1/8 - .7,
                                           #mu=c(9:5/9-.5,-(0:3 )),
                                           mu=9:1 / 9 -.5,
                                           l_bound = 0,
                                           prior_rho = c(9.53, 9.66),
                                           #prior_rho = c(5.5, 5.5), #prior_rho = c(9.4,8.4), #prior_rho = c(26.2,19.7), #prior_rho = c(14.9, 11.2), #prior_rho = c(8.8, 11.8),
                                           plot_graphs = T,
                                           iter = iter, warmup=warmup, adapt_delta = adapt_delta, max_treedepth =max_treedepth, chains=chains, 
                                           formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                           formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                           # formula_x = y_scaled ~ AccidentYear_scaled + beta_warp(DevelopmentLag, 1,4) -1,
                                           #formula_H = y_scaled ~ AccidentYear_scaled + beta_warp(DevelopmentLag, 1,4) -1,
                                           stanFile = "gp_hurdle_02_phiapprox2.stan"
                                           #stanFile = "gp_hurdle_02.stan"
                                           #stanFile = "gp_hurdle_02_phiapprox.stan"
                                           #stanFile = "gp_hurdle_02_alt_trunc.stan"
                                           #stanFile = "gp_hurdle_zmethod_01.stan"
  )
  b=Sys.time(); b-a
  ans_hurdle$summary_stats
  ans_hurdle$errors
  plot(ans_hurdle$errors$rmse_step_ahead)
  
  
  ##### 
  # Vector mu with virtual data
  adapt_delta = 0.95; iter = 800; warmup=700; max_treedepth=11; chains=6
  #adapt_delta = 0.95; iter = 30; warmup=20; max_treedepth=10; chains=1
  a=Sys.time()
  ans_hurdle_virt = stan_incremental_loss_ratio_virtual(loss_triangle = loss$paid_tr,  loss_square = loss$paid_sq, premiums = premiums, calc_ranks = T,
                                            virtual_data = T, virtual_point = -.3, 
                                            mu=9:0 / 9 -.5,
                                            prior_sigma = .5, 
                                            l_bound = 0,
                                            prior_rho = c(9.53, 9.66),
                                            #prior_rho = c(5.5, 5.5), #prior_rho = c(9.4,8.4), #prior_rho = c(26.2,19.7), #prior_rho = c(14.9, 11.2), #prior_rho = c(8.8, 11.8),
                                            plot_graphs = T,
                                            iter = iter, warmup=warmup, adapt_delta = adapt_delta, max_treedepth =max_treedepth, 
                                            formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                            formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                            stanFile = "gp_hurdle_virtual_02.stan"
  )
  b=Sys.time(); b-a
  ans_hurdle_virt$errors
  ans_hurdle$errors
  ans_hurdle$summary_stats
  ans_hurdle_virt$summary_stats
  ans_hurdle_virt$loss_square / ans_hurdle$loss_square
  plot(ans_hurdle_virt$errors$rmse_step_ahead)
  
  
  # Vector mu with bias factor in hurdle model
  adapt_delta = 0.95; iter = 500; warmup=400; max_treedepth=10; chains=6
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
                                            # formula_x = y_scaled ~ AccidentYear_scaled + beta_warp(DevelopmentLag, 1,4) -1,
                                            #formula_H = y_scaled ~ AccidentYear_scaled + beta_warp(DevelopmentLag, 1,4) -1,
                                            stanFile = "gp_hurdle_02_bias.stan"
                                            #stanFile = "gp_hurdle_zmethod_01.stan"
  )
  b=Sys.time(); b-a
  ans_hurdle$errors
  plot(ans_hurdle$errors$rmse_step_ahead)
  
  
  ##### 
  # Vector mu / NO POSITIVE ORDERED CONSTRAINT ON SIGMA but SUCCESSIVELY MORE AGGRESSIVE SHRINKAGE ON SIGMA
  adapt_delta = 0.95; iter = 300; warmup=200; max_treedepth=12; chains=6
  a=Sys.time()
  ans_hurdle_shrink_sigma = stan_incremental_loss_ratio2(loss_triangle = loss$paid_tr,  loss_square = loss$paid_sq, premiums = premiums, calc_ranks = T,
                                                      prior_sigma = .5, 
                                                      #mu = rep(0,9),
                                                      #mu=9:1 / 5 -1.3,
                                                      #mu=9:1/8 - .7,
                                                      #mu=c(9:5/9-.5,-(0:3 )),
                                                      mu=9:1 / 9 -.5,
                                                      l_bound = 0,
                                                      prior_rho = c(9.53, 9.66),
                                                      #prior_rho = c(5.5, 5.5), #prior_rho = c(9.4,8.4), #prior_rho = c(26.2,19.7), #prior_rho = c(14.9, 11.2), #prior_rho = c(8.8, 11.8),
                                                      plot_graphs = T,
                                                      iter = iter, warmup=warmup, adapt_delta = adapt_delta, max_treedepth =max_treedepth, 
                                                      formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                                      formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                                      # formula_x = y_scaled ~ AccidentYear_scaled + beta_warp(DevelopmentLag, 1,4) -1,
                                                      #formula_H = y_scaled ~ AccidentYear_scaled + beta_warp(DevelopmentLag, 1,4) -1,
                                                      stanFile = "gp_hurdle_02_shrinkage_sigma.stan"
  )
  b=Sys.time(); b-a
  ans_hurdle_shrink_sigma$errors
  plot(ans_hurdle_shrink_sigma$errors$rmse_step_ahead)
  
  
  ##### 
  # Vector mu / NO POSITIVE ORDERED CONSTRAINT ON SIGMA
  adapt_delta = 0.95; iter = 300; warmup=200; max_treedepth=10; chains=6
  a=Sys.time()
  ans_hurdle_unc_sigma = stan_incremental_loss_ratio2(loss_triangle = loss$paid_tr,  loss_square = loss$paid_sq, premiums = premiums, calc_ranks = T,
                                            prior_sigma = .5, 
                                            #mu = rep(0,9),
                                            #mu=9:1 / 5 -1.3,
                                            #mu=9:1/8 - .7,
                                            #mu=c(9:5/9-.5,-(0:3 )),
                                            mu=9:1 / 9 -.5,
                                            l_bound = 0,
                                            prior_rho = c(9.53, 9.66),
                                            #prior_rho = c(5.5, 5.5), #prior_rho = c(9.4,8.4), #prior_rho = c(26.2,19.7), #prior_rho = c(14.9, 11.2), #prior_rho = c(8.8, 11.8),
                                            plot_graphs = T,
                                            iter = iter, warmup=warmup, adapt_delta = adapt_delta, max_treedepth =max_treedepth, 
                                            formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                            formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                            # formula_x = y_scaled ~ AccidentYear_scaled + beta_warp(DevelopmentLag, 1,4) -1,
                                            #formula_H = y_scaled ~ AccidentYear_scaled + beta_warp(DevelopmentLag, 1,4) -1,
                                            stanFile = "gp_hurdle_02_unconstrained_sigma.stan"
  )
  b=Sys.time(); b-a
  ans_hurdle_unc_sigma$errors
  plot(ans_hurdle_unc_sigma$errors$rmse_step_ahead)
  
  
  
  ##### 
  # Vector mu with logistic selection of hurdle
  adapt_delta = 0.95; iter = 500; warmup=400; max_treedepth=10; chains=6
  a=Sys.time()
  ans_hurdle = stan_incremental_loss_ratio2(loss_triangle = loss$paid_tr,  loss_square = loss$paid_sq, premiums = premiums, calc_ranks = T,
                                            prior_sigma = .5, 
                                            #mu=rep(0,9),
                                            mu=9:1 / 5 -1.3,
                                            #mu=c(9:5/9-.5,-(0:3 /1.1)),
                                            l_bound = 0,
                                            #prior_rho = c(5.5, 5.5),      #prior_rho = c(9.4,8.4), 
                                            prior_rho = c(9.53, 9.66),
                                            #prior_rho = c(26.2,19.7),  #prior_rho = c(14.9, 11.2),
                                            plot_graphs = T,
                                            iter = iter, warmup=warmup, adapt_delta = adapt_delta, max_treedepth =max_treedepth, 
                                            formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                            formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                            stanFile = "gp_hurdle_logistic_02.stan"
                                            #stanFile = "gp_hurdle_zmethod_01.stan"
  )
  b=Sys.time(); b-a
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_hurdle$loss_square/premiums,actual_tr = loss$paid_tr/premiums)
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_hurdle$loss_square_med/premiums,actual_tr = loss$paid_tr/premiums)
  colMeans(loss$paid_sq[,-1]/premiums - loss$paid_sq[,-10]/premiums <=0)
  
  
  ##### 
  # No linear component version
  adapt_delta = 0.95; iter = 700; warmup=600; max_treedepth=10; chains=6
  ans_hurdle = stan_incremental_loss_ratio(loss_triangle = loss$paid_tr,  loss_square = loss$paid_sq, premiums = premiums, 
                                           prior_sigma = .25, mu=-0.5,
                                           prior_rho = c(5.5, 5.5), 
                                           #prior_rho = c(9.53, 9.66),
                                           #prior_rho = c(9.4,8.4), 
                                           #prior_rho = c(14.9, 11.2),
                                           #prior_rho = c(26.2,19.7), 
                                           plot_graphs = T,
                                           iter = iter, warmup=warmup, adapt_delta = adapt_delta, max_treedepth =max_treedepth, chains=chains,
                                           formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                           formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                           #stanFile = "gp_hurdle_no_linear_01.stan"
                                           stanFile = "gp_hurdle_01.stan"
  )
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_hurdle/premiums,actual_tr = loss$paid_tr/premiums)

  ## Analyze All Business Lines
  
  adapt_delta = 0.9; iter = 800; warmup=700; max_treedepth=10; chains=6
  results =list()
  for (bus in business_lines){
    results[[bus]] =  business_line_analysis(projection_fun =  stan_incremental_loss_ratio, business_line = bus, take_logs=F, 
                                             prior_sigma = .5, mu=-.5,
                                             #prior_rho = c(26.2,19.7), 
                                             prior_rho = c(9.4,8.4), 
                                             #prior_rho = c(9.53, 9.66),
                                             iter = iter, warmup=warmup, chains= chains,
                                             calc_ranks = T, plot_graphs = F,
                                             formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag -1,
                                             formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                             stanFile = "gp_hurdle_01.stan"
    )
  }
  save(results, file = paste0(output_dir, "results_hurdle.RData"))
  #load(file = paste0(output_dir, "results_hurdle.RData"))
  # RMSE Vals
  summary_rmse = rbind(
    map_df(results, ~ median(.x$results$rmse_loss_ratio, na.rm = T)),
    map_df(results, ~ mean(.x$results$rmse_loss_ratio, na.rm = T)),
    map_df(results, ~ median(.x$results$rmse_loss_ratio_med, na.rm = T)),
    map_df(results, ~ mean(.x$results$rmse_loss_ratio_med, na.rm = T))
  ) %>% as.data.frame()
  rownames(summary_rmse) = c("Median RMSE of Mean LR", "Mean RMSE of Mean LR",  "Median RMSE of Med. LR", "Mean RMSE of Med. LR")
  summary_rmse
  
  # Percentage Rank Plots
  plot_ranks(results)
  
  results2 =list()
  for (bus in business_lines){
    results2[[bus]] =  business_line_analysis(projection_fun =  stan_incremental_loss_ratio, business_line = bus, take_logs=F, 
                                             prior_sigma = .25, mu=-.2,
                                             #prior_rho = c(26.2,19.7), 
                                             prior_rho = c(9.4,8.4), 
                                             
                                             #prior_rho = c(9.53, 9.66),
                                             iter = iter, warmup=warmup, chains= chains,
                                             calc_ranks = T, plot_graphs = F,
                                             formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag -1,
                                             formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                             stanFile = "gp_hurdle_01.stan"
    )
  }
  save(results2, file = paste0(output_dir, "results2_hurdle.RData"))
  #load(file = paste0(output_dir, "results2_hurdle.RData"))
  # RMSE Vals
  map_df(results2, ~ median(.x$results$rmse_loss_ratio, na.rm = T))
  map_df(results2, ~ mean(.x$results$rmse_loss_ratio, na.rm = T))
  map_df(results2, ~ median(.x$results$rmse_loss_ratio_med, na.rm = T))
  map_df(results2, ~ mean(.x$results$rmse_loss_ratio_med, na.rm = T))
  
  
  # Percentage Rank Plots
  plot_ranks(results2)
  
  ################################################################
  ## Test Across Business Lines
  # Vector mu
  adapt_delta = 0.95; iter = 800; warmup=700; max_treedepth=10; chains=6
  results3 =list()
  for (bus in business_lines){
    results3[[bus]] =  business_line_analysis(projection_fun =  stan_incremental_loss_ratio2, business_line = bus, take_logs=F, 
                                             prior_sigma = .5, 
                                             mu=9:1 / 9 -.5,
                                             prior_rho = c(9.4,8.4), 
                                             iter = iter, warmup=warmup, chains= chains,
                                             calc_ranks = T, plot_graphs = F,
                                             formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag -1,
                                             formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                             stanFile = "gp_hurdle_02.stan"
    )
  }
  #save(results3, file = paste0(output_dir, "results3_hurdle.RData"))
  #load(file = paste0(output_dir, "results3_hurdle.RData"))
  # RMSE Vals
  summarize_rmse(results3)
  
  
  # Percentage Rank Plots
  plot_ranks(results3)
  
  
  #####
  # Vector mu with virtual data
  adapt_delta = 0.95; iter = 1000; warmup=900; max_treedepth=11; chains=6
  results4 =list()
  for (bus in business_lines){
    results4[[bus]] =  business_line_analysis(projection_fun =  stan_incremental_loss_ratio_virtual, business_line = bus, take_logs=F, 
                                              prior_sigma = .5, 
                                              virtual_data = T, virtual_point = -.3, 
                                              mu=9:0 / 9 -.5,
                                              prior_rho = c(9.4,8.4), 
                                              iter = iter, warmup=warmup, chains= chains,
                                              calc_ranks = T, plot_graphs = F,
                                              formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag -1,
                                              formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                              stanFile = "gp_hurdle_virtual_02.stan"
    )
  }
  #save(results4, file = paste0(output_dir, "results4_hurdle.RData"))
  #load(file = paste0(output_dir, "results4_hurdle.RData"))
  # RMSE Vals
 
  summarize_rmse(results4)
  # Percentage Rank Plots
  plot_ranks(results4)
  
  allplots = plot_ranks(results4)
  map(allplots, ~.x$plot_step_ahead_cum2)
  map(allplots, ~.x$plot_step_ahead_cum)
  allplots3 = plot_ranks(results3)
  map(allplots3, ~.x$plot_step_ahead_cum2)
  map(allplots3, ~.x$plot_step_ahead_cum)
  
  allplots3$comauto$plot_step_ahead_inc
  allplots$comauto$plot_step_ahead_inc
  
  allplots3$comauto$plot_step_ahead_cum
  allplots3$comauto$plot_step_ahead_cum2
  allplots$comauto$plot_step_ahead_cum
  allplots$comauto$plot_step_ahead_cum2
}  


##########################################################################
####  Incremental LDF Method with HURDLE MODEL

ldf_hurdle = function(){
  
  
  ## Test Single Case
  
  bus = "wkcomp"
  codes = unlist(unique(data %>% filter(line==bus) %>% select(GRCODE)))
  company_code = 353
  loss = loss_table(data, business_line = bus, company_code = company_code)
  premiums = get_premium(data, business_line = bus, company_code = company_code)
  
  ##### 
  # Vector mu
  adapt_delta = 0.95; iter = 500; warmup=400; max_treedepth=10; chains=6
  a=Sys.time()
  ans_hurdle_ldf = stan_ldf2(loss_triangle = loss$paid_tr,  loss_square = loss$paid_sq, premiums = premiums, calc_ranks = T,
                                            prior_sigma = .5, 
                                            mu=9:1 / 9 +.5,
                                            l_bound = 1,
                                            prior_rho = c(9.53, 9.66),
                                            plot_graphs = T,
                                            iter = iter, warmup=warmup, adapt_delta = adapt_delta, max_treedepth =max_treedepth, 
                                            formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                            formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                            stanFile = "gp_hurdle_02.stan"
                                            #stanFile = "gp_hurdle_zmethod_01.stan"
  )
  b=Sys.time(); b-a
  ans_hurdle_ldf$errors
  plot(ans_hurdle_ldf$errors$rmse_step_ahead)

  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_warp/premiums,actual_tr = loss$paid_tr/premiums)
  
  ## Analyze All Business Lines
  
  adapt_delta = 0.9; iter = 800; warmup=700; max_treedepth=10; chains=6
  results =list()
  for (bus in business_lines){
    results[[bus]] =  business_line_analysis(projection_fun =  stan_ldf2, business_line = bus, take_logs=F, 
                                             prior_sigma = .5, 
                                             mu=9:1 / 9 +.5,
                                             l_bound = 1,
                                             #prior_rho = c(26.2,19.7), 
                                             prior_rho = c(9.4,8.4), 
                                             #prior_rho = c(9.53, 9.66),
                                             iter = iter, warmup=warmup, chains= chains,
                                             calc_ranks = T, plot_graphs = F,
                                             formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag -1,
                                             formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                             stanFile = "gp_hurdle_02.stan"
    )
  }
  save(results, file = paste0(output_dir, "results_ldf_hurdle.RData"))
  #load(file = paste0(output_dir, "results_ldf_hurdle.RData"))
  # RMSE Vals
  summarize_rmse(results)
  
  
  # Percentage Rank Plots
  plot_ranks(results)
  
  results2 =list()
  for (bus in business_lines){
    results2[[bus]] =  business_line_analysis(projection_fun =  stan_incremental_loss_ratio, business_line = bus, take_logs=F, 
                                              prior_sigma = .25, mu=-.2,
                                              #prior_rho = c(26.2,19.7), 
                                              prior_rho = c(9.4,8.4), 
                                              
                                              #prior_rho = c(9.53, 9.66),
                                              iter = iter, warmup=warmup, chains= chains,
                                              calc_ranks = T, plot_graphs = F,
                                              formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag -1,
                                              formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                              stanFile = "gp_hurdle_01.stan"
    )
  }
  save(results2, file = paste0(output_dir, "results2_hurdle.RData"))
  #load(file = paste0(output_dir, "results2_hurdle.RData"))
  # RMSE Vals
  map_df(results2, ~ median(.x$results$rmse_loss_ratio, na.rm = T))
  map_df(results2, ~ mean(.x$results$rmse_loss_ratio, na.rm = T))
  map_df(results2, ~ median(.x$results$rmse_loss_ratio_med, na.rm = T))
  map_df(results2, ~ mean(.x$results$rmse_loss_ratio_med, na.rm = T))
  
  
  # Percentage Rank Plots
  plot_ranks(results2)
  
  #####
  # Vector mu
  adapt_delta = 0.95; iter = 500; warmup=400; max_treedepth=10; chains=6
  results3 =list()
  for (bus in business_lines){
    results3[[bus]] =  business_line_analysis(projection_fun =  stan_incremental_loss_ratio2, business_line = bus, take_logs=F, 
                                              prior_sigma = .5, 
                                              mu=9:1 / 9 -.5,
                                              prior_rho = c(9.4,8.4), 
                                              iter = iter, warmup=warmup, chains= chains,
                                              calc_ranks = T, plot_graphs = F,
                                              formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag -1,
                                              formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                              stanFile = "gp_hurdle_02.stan"
    )
  }
  #save(results3, file = paste0(output_dir, "results3_hurdle.RData"))
  #load(file = paste0(output_dir, "results3_hurdle.RData"))
  # RMSE Vals
  summary_rmse = rbind(
    map_df(results3, ~ median(.x$results$rmse_loss_ratio, na.rm = T)),
    map_df(results3, ~ mean(.x$results$rmse_loss_ratio, na.rm = T)),
    map_df(results3, ~ median(.x$results$rmse_loss_ratio_med, na.rm = T)),
    map_df(results3, ~ mean(.x$results$rmse_loss_ratio_med, na.rm = T))
  ) %>% as.data.frame()
  rownames(summary_rmse) = c("Median RMSE of Mean LR", "Mean RMSE of Mean LR",  "Median RMSE of Med. LR", "Mean RMSE of Med. LR")
  (summary_rmse = apply(summary_rmse,2, round, digits = 3))
  
  # Percentage Rank Plots
  plot_ranks(results3)
  
  
}  


##########################################################################
####  LDF Method with NO trunction, WITH Logarithmic Transform

ldf_no_trunc = function(){
  
  
  ## Test Single Case
  
  bus = "ppauto"
  codes = unlist(unique(data %>% filter(line==bus) %>% select(GRCODE)))
  company_code = 15660
  loss = loss_table(data, business_line = bus, company_code = company_code)
  premiums = get_premium(data, business_line = bus, company_code = company_code)
  
  # Using special log transform 
  ans_mean_het = stan_ldf(loss_triangle = loss$paid_tr, loss_square = loss$paid_sq, premiums = premiums, 
                          calc_ranks = T,
                          prior_sigma = 1, plot_graphs = T,
                          take_logs = TRUE,
                           iter = iter, warmup=warmup, 
                           formula_x = y_log_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                           formula_H = y_log_scaled ~ log(DevelopmentLag_scaled - min(DevelopmentLag_scaled)+1) + AccidentYear_scaled-1, 
                           #formula_H = y_log_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                           stanFile = "gp_compound_mean_het_03.stan"
  )
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_mean_het$loss_square/premiums,actual_tr = loss$paid_tr/premiums)
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_mean_het$loss_square_2/premiums,actual_tr = loss$paid_tr/premiums)

  ### Old Method using a mean function
  ans_old = stan_ladder(loss_triangle = loss$paid_tr, 
                     formula_x = dev_log ~ AccidentYear_scale + DevelopmentLag_scale -1,
                     formula_H = ~ log(DevelopmentLag_scale - min(DevelopmentLag_scale)+1) + AccidentYear_scale, 
                     eps = 0.001, adapt_delta = 0.975, iter = 400, prior_theta=1)
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_old/premiums,actual_tr = loss$paid_tr/premiums)
 
  ### Old Method using a compound kernel function:  this should produce similar results to stan_ldf for "loss_square_2"
  ans_old2 = stan_ladder_se(loss_triangle = loss$paid_tr, 
                            formula_H = ~ log(DevelopmentLag_scaled - min(DevelopmentLag_scaled)+1) + AccidentYear_scaled-1, 
                     eps = 0.0005, plot_graphs = FALSE,
                     chains = 6,iter = 400, warmup = 100, adapt_delta = 0.9,max_treedepth = 10,
                     prior_eta = 1, prior_rho = c(5.557,5.564), prior_sigma = 1, 
                     stanFile = "gp_se_03.stan")
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_old2/premiums,actual_tr = loss$paid_tr/premiums)
  
  ans_warp =  stan_ldf(loss_triangle = loss$paid_tr, loss_square = loss$paid_sq, premiums = premiums, 
                                       calc_ranks = T,
                                       prior_sigma = 1, plot_graphs = T,
                                       take_logs = TRUE,
                                       iter = iter, warmup=warmup, 
                                       formula_x = y_log_scaled ~ AccidentYear_scaled + beta_warp(DevelopmentLag, 1,5) -1,
                                       formula_H = y_log_scaled ~ AccidentYear_scaled + beta_warp(DevelopmentLag, 1,5) -1,
                                       stanFile = "gp_compound_mean_het_03.stan"
  )
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_warp$loss_square/premiums,actual_tr = loss$paid_tr/premiums)
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_warp$loss_square_2/premiums,actual_tr = loss$paid_tr/premiums)
  
  
  ## Analyze All Business Lines
  
  results =list()
  for (bus in business_lines){
    results[[bus]] =  business_line_analysis(projection_fun =  stan_ldf, business_line = bus, take_logs=T, 
                                             prior_sigma = .5, iter = iter, warmup=warmup, chains= chains,
                                             calc_ranks = T, plot_graphs = F,
                                             formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                             formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                             stanFile = "gp_compound_mean_het_03.stan"
    )
  }
  #save(results, file = paste0(output_dir, "results_incremenal.RData"))
  #load(file = paste0(output_dir, "results_incremenal.RData"))
  # RMSE Vals
  map(results, ~ median(.$results$rmse_loss_ratio, na.rm = T))
  map(results, ~ mean(.x$results$rmse_loss_ratio, na.rm = T))
  
  # Percentage Rank Plots
  ranks = plot_ranks(results)
  ranks$rank_scores
}  


##########################################################################
####  LDF Method NO trunction, NO Logarithmic Transform

ldf_no_trunc_no_log = function(){
  
  
  ## Test Single Case
  
  bus = "ppauto"
  codes = unlist(unique(data %>% filter(line==bus) %>% select(GRCODE)))
  company_code = 15660
  loss = loss_table(data, business_line = bus, company_code = company_code)
  premiums = get_premium(data, business_line = bus, company_code = company_code)
  
  ans_mean_het = stan_ldf(loss_triangle = loss$paid_tr, loss_square = loss$paid_sq, premiums = premiums, 
                          calc_ranks = T,
                          prior_sigma = 1, plot_graphs = T,
                          take_logs = FALSE,
                          iter = iter, warmup=warmup, 
                          formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                          formula_H = y_scaled ~ log(DevelopmentLag_scaled - min(DevelopmentLag_scaled)+1) + AccidentYear_scaled-1, 
                          #formula_H = y_log_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                          stanFile = "gp_compound_mean_het_03.stan"
  )
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_mean_het$loss_square/premiums,actual_tr = loss$paid_tr/premiums)
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_mean_het$loss_square_2/premiums,actual_tr = loss$paid_tr/premiums)
  
  ans_warp =  stan_ldf(loss_triangle = loss$paid_tr, loss_square = loss$paid_sq, premiums = premiums, 
                       calc_ranks = T,
                       prior_sigma = 1, plot_graphs = T,
                       take_logs = FALSE,
                       iter = iter, warmup=warmup, 
                       formula_x = y_scaled ~ AccidentYear_scaled + beta_warp(DevelopmentLag, 1,5) -1,
                       formula_H = y_scaled ~ AccidentYear_scaled + beta_warp(DevelopmentLag, 1,5) -1,
                       stanFile = "gp_compound_mean_het_03.stan"
  )
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_warp$loss_square/premiums,actual_tr = loss$paid_tr/premiums)
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_warp$loss_square_2/premiums,actual_tr = loss$paid_tr/premiums)
  
  
  ## Analyze All Business Lines
  # With a linear component to kernel
  results =list()
  for (bus in business_lines){
    results[[bus]] =  business_line_analysis(projection_fun =  stan_ldf, business_line = bus, 
                                             take_logs=FALSE, 
                                             prior_sigma = .5, iter = iter, warmup=warmup, chains= chains,
                                             calc_ranks = T, plot_graphs = F,
                                             formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                             formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                             stanFile = "gp_compound_mean_het_03.stan"
    )
  }
  #save(results, file = paste0(output_dir, "results_incremenal.RData"))
  #load(file = paste0(output_dir, "results_incremenal.RData"))
  # RMSE Vals
  map(results, ~ median(.$results$rmse_loss_ratio, na.rm = T))
  map(results, ~ mean(.x$results$rmse_loss_ratio, na.rm = T))
  
  # Percentage Rank Plots
  plot_ranks(results)
  
  
  # No linear component to kernel
  results =list()
  for (bus in business_lines){
    results[[bus]] =  business_line_analysis(projection_fun =  stan_ldf, business_line = bus, 
                                             take_logs=FALSE, 
                                             prior_sigma = .5, iter = iter, warmup=warmup, chains= chains,
                                             calc_ranks = T, plot_graphs = F,
                                             formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                             formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1, # !!! Ignored in Stan
                                             stanFile = "gp_compound_mean_het_no_linear_03.stan"
    )
  }
  #save(results, file = paste0(output_dir, "results_incremenal.RData"))
  #load(file = paste0(output_dir, "results_incremenal.RData"))
  # RMSE Vals
  map(results, ~ median(.$results$rmse_loss_ratio, na.rm = T))
  map(results, ~ mean(.x$results$rmse_loss_ratio, na.rm = T))
  
  # Percentage Rank Plots
  plot_ranks(results)
  
  
}  



##########################################################################
####  LDF Method WITH trunction, NO Logarithmic Transform

ldf_with_trunc_no_log = function(){
  
  
  ## Test Single Case
  
  bus = "ppauto"
  codes = unlist(unique(data %>% filter(line==bus) %>% select(GRCODE)))
  company_code = 15660
  loss = loss_table(data, business_line = bus, company_code = company_code)
  premiums = get_premium(data, business_line = bus, company_code = company_code)
  
  ans_mean_het = stan_ldf(loss_triangle = loss$paid_tr, loss_square = loss$paid_sq, premiums = premiums, 
                          calc_ranks = T,
                          prior_sigma = 1, plot_graphs = T,
                          take_logs = FALSE,adapt_delta=adapt_delta,
                          iter = iter, warmup=warmup, 
                          formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                          formula_H = y_scaled ~ log(DevelopmentLag_scaled - min(DevelopmentLag_scaled)+1) + AccidentYear_scaled-1, 
                          #formula_H = y_log_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                          stanFile = "gp_comp_het_trunc_03.stan"
  )
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_mean_het$loss_square/premiums,actual_tr = loss$paid_tr/premiums)
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_mean_het$loss_square_2/premiums,actual_tr = loss$paid_tr/premiums)
  
  ans_warp =  stan_ldf(loss_triangle = loss$paid_tr, loss_square = loss$paid_sq, premiums = premiums, 
                       calc_ranks = T,
                       prior_sigma = 1, plot_graphs = T,
                       take_logs = FALSE,
                       iter = iter, warmup=warmup, 
                       formula_x = y_scaled ~ AccidentYear_scaled + beta_warp(DevelopmentLag, 1,5) -1,
                       formula_H = y_scaled ~ AccidentYear_scaled + beta_warp(DevelopmentLag, 1,5) -1,
                       stanFile = "gp_comp_het_trunc_03.stan"
  )
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_warp$loss_square/premiums,actual_tr = loss$paid_tr/premiums)
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_warp$loss_square_2/premiums,actual_tr = loss$paid_tr/premiums)
  
  
  ## Analyze All Business Lines
  
  results =list()
  for (bus in business_lines){
    results[[bus]] =  business_line_analysis(projection_fun =  stan_ldf, business_line = bus, 
                                             take_logs=FALSE, 
                                             prior_sigma = .5, iter = iter, warmup=warmup, chains= chains,
                                             calc_ranks = T, plot_graphs = F,
                                             formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                             formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                             stanFile = "gp_comp_het_trunc_03.stan"
    )
  }
  #save(results, file = paste0(output_dir, "results_incremenal.RData"))
  #load(file = paste0(output_dir, "results_incremenal.RData"))
  # RMSE Vals
  map(results, ~ median(.$results$rmse_loss_ratio, na.rm = T))
  map(results, ~ mean(.x$results$rmse_loss_ratio, na.rm = T))
  
  # Percentage Rank Plots
  plot_ranks(results)
  
}  


##########################################################################
####  TWO STEP Approach:  Use Incremental Loss Ratio for Mean Field - then estimate sigma


two_step_residual = function(){
  
  
  ## Test Single Case
  
  bus = "prodliab"
  codes = unlist(unique(data %>% filter(line==bus) %>% select(GRCODE)))
  company_code = 78
  loss = loss_table(data, business_line = bus, company_code = company_code)
  premiums = get_premium(data, business_line = bus, company_code = company_code)
  
  ans_mean_het = residual_mod_incremental_lr(loss_triangle = loss$paid_tr, loss_square = loss$paid_sq, premiums = premiums, 
                                             calc_ranks = T,
                                             prior_sigma = 1, plot_graphs = T,
                                             iter = iter, warmup=warmup,
                                             formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                             formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                             stanFile = "gp_compound_mean_het_03.stan",
                                             stanFileResidual = "residual_01.stan"
  )
  
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_mean_het$loss_square/premiums,actual_tr = loss$paid_tr/premiums)
  
  
  ans_warp = residual_mod_incremental_lr(loss_triangle = loss$paid_tr, loss_square = loss$paid_sq,premiums = premiums, 
                                         calc_ranks = T,prior_sigma = 1, plot_graphs = T,
                                         iter = iter, warmup=warmup,
                                         formula_x = y_scaled ~ AccidentYear_scaled + beta_warp(DevelopmentLag, 1,5) -1,
                                         formula_H = y_scaled ~ AccidentYear_scaled + beta_warp(DevelopmentLag, 1,5) -1,
                                         stanFile = "gp_compound_mean_het_no_linear_03.stan",
                                         stanFileResidual = "residual_01.stan"
  )
  rmse_triangle(actual_sq = loss$paid_sq/premiums, predicted_sq = ans_warp$loss_square/premiums,actual_tr = loss$paid_tr/premiums)
  
  
  ################################
  ## Analyze All Business Lines
  
  ## With Linear Component
  results =list()
  for (bus in business_lines){
    results[[bus]] =  business_line_analysis(projection_fun =  residual_mod_incremental_lr, business_line = bus, take_logs=F, 
                                             prior_sigma = .5, iter = iter, warmup=warmup, chains= chains,
                                             calc_ranks = T, plot_graphs = F,
                                             formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                             formula_H = y_scaled ~ log(DevelopmentLag_scaled - min(DevelopmentLag_scaled)+1) + AccidentYear_scaled-1, 
                                             #formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                             stanFile = "gp_compound_mean_het_03.stan",
                                             stanFileResidual = "residual_01.stan"
    )
  }
  #save(results, file = paste0(output_dir, "results_incremenal.RData"))
  #load(file = paste0(output_dir, "results_incremenal.RData"))
  
  # RMSE Vals
  map(results, ~ median(.$results$rmse_loss_ratio, na.rm = T))
  map(results, ~ mean(.x$results$rmse_loss_ratio, na.rm = T))
  
  # Percentage Rank Plots
  plot_ranks(results)
  
  ## No Linear Component
  results =list()
  for (bus in business_lines){
    results[[bus]] =  business_line_analysis(projection_fun =  residual_mod_incremental_lr, business_line = bus, take_logs=F, 
                                             prior_sigma = .5, iter = iter, warmup=warmup, chains= chains,
                                             calc_ranks = T, plot_graphs = F,
                                             formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                             formula_H = y_scaled ~ log(DevelopmentLag_scaled - min(DevelopmentLag_scaled)+1) + AccidentYear_scaled-1, ## IGNORED in Stan
                                             #formula_H = y_scaled ~ AccidentYear_scaled + log(DevelopmentLag) -1,
                                             stanFile = "gp_compound_mean_het_no_linear_03.stan",
                                             stanFileResidual = "residual_01.stan"
    )
  }
  #save(results, file = paste0(output_dir, "results_incremental_2step.RData"))
  #load(file = paste0(output_dir, "results_incremental_2step.RData"))
  # RMSE Vals
  map(results, ~ median(.$results$rmse_loss_ratio, na.rm = T))
  map(results, ~ mean(.x$results$rmse_loss_ratio, na.rm = T))
  
  # Percentage Rank Plots
  plot_ranks(results)
  
  
}  





