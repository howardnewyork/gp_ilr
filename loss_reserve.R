###########################################################
#################  Loss Reserve Analysis
#################  V 0.3
###########################################################

# Set up libraries

library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(scoringRules)
library(ChainLadder)

#packages for custom GP

### Greta, mvtnorm, coda, bayesplot are required if Greta functions needed
#library(greta)
#library(mvtnorm)
#library(coda)
#library (bayesplot)



library(rstan)
#options(mc.cores = min(10,parallel::detectCores()))
options(mc.cores = min(parallel::detectCores()))
rstan_options(auto_write = TRUE)

library(DiceKriging)
library(invgamma)

library(kolmim) # Gets CDF of the Kolmogorov distribution
library(gridExtra)


# Set up directories
p0=paste0
working_dir = p0(getwd(),"/")
data_dir = p0(working_dir, "data/")
output_dir = p0(working_dir, "output/")
stan_dir = p0(working_dir,"stan/")


## Source file
source('custom_gp.R')
source('utilities.R')

#' Load loss data into memory
#' 
#' @param directory name of data directory where files are stored.  Must end with "/"
#' @param files
load_losses = function(directory = data_dir, files = c("comauto_pos.csv", "medmal_pos.csv", 
                                                 "othliab_pos.csv","ppauto_pos.csv",
                                                 "prodliab_pos.csv", "wkcomp_pos.csv")){

  
  # Get business line names
  business_lines = str_split(files, "_") %>%
    map(1) %>%
    unlist()
  
  read_named = function(file, busines_line){
    data= read_csv(file)
    data$line = busines_line
    cnames = colnames(data)  %>%
      str_split("_") %>%
      map(1) %>%
      unlist()
    colnames(data) = cnames
    data
  }
  
  files = p0(directory,files)
  data = files %>%
    map2(business_lines, read_named)
  data %>% map(colnames)
  head(data)
  data = do.call(rbind, data)
  
  return(data)
}

load_losses_subset = function(min_tot_premium=10000){
  data = load_losses()
  
  # Select cases with 100 cells for product line and company (i.e. data with a full "square")
  
  sum_data = data %>%
    group_by(GRCODE, GRNAME, line) %>%
    mutate(count = sum(EarnedPremNet > 0) + sum(CumPaidLoss >0)) %>%
    filter(count == 200, DevelopmentLag == 1, AccidentYear == 1988) %>%
    dplyr::select(GRCODE, GRNAME, line) 
  
  
  sum_data = sum_data %>%
    mutate(group = paste0(GRCODE,"_",line))
  
  data = data %>%
    mutate(group = paste0(GRCODE,"_",line)) %>%
    filter(group %in% sum_data$group)
  data
}

#' Select a particular business line
#' 
#' @param data Losses database
#' @param business_line Name of business line to be selected.  Options include: "comauto", "medmal", "othliab","ppauto","prodliab", "wkcomp"
select_line = function(data, business_line){
  data %>% filter(line == business_line)
}




#' Create a table of losses, rows = accident year, cols = developmentyear
#' 
#' @param data Losses database
#' @param business_line Name of business line to be selected
#' @param company_code GRCODE value
loss_table = function(data, business_line, company_code){
  incurred_sq = 
    data %>% 
    ungroup() %>%
    select_line(business_line) %>%
    filter(GRCODE == company_code) %>%
    dplyr::select(AccidentYear, DevelopmentLag, IncurLoss) %>%
    spread(DevelopmentLag, IncurLoss) %>%
    as.matrix()
  rownames(incurred_sq) = incurred_sq[,"AccidentYear"]
  incurred_sq = incurred_sq[,-1] 
  
  paid_sq = 
    data %>%
    ungroup() %>%
    select_line(business_line) %>%
    filter(GRCODE == company_code) %>%
    dplyr::select(AccidentYear, DevelopmentLag, CumPaidLoss) %>%
    spread(DevelopmentLag, CumPaidLoss)%>%
    as.matrix()
  rownames(paid_sq) = paid_sq[,"AccidentYear"]
  paid_sq = paid_sq[,-1] 
   
  incurred_tr = incurred_sq
  for (i in 2:nrow(incurred_tr)){
    for (j in ncol(incurred_tr):(ncol(incurred_tr)-i+2))
      incurred_tr[i,j] = NA
  }
  incurred_tr = apply(incurred_tr, 2, function(x) ifelse(x == 0, NA, x))
  
  
  paid_tr = paid_sq
  for (i in 2:nrow(paid_tr)){
    for (j in ncol(paid_tr):(ncol(paid_tr)-i+2))
      paid_tr[i,j] = NA
  }
  paid_tr = apply(paid_tr, 2, function(x) ifelse(x == 0, NA, x))
  
  return(
    list(incurred_sq = incurred_sq, paid_sq= paid_sq, incurred_tr=incurred_tr, paid_tr=paid_tr)
  )
  
}


#' Format loss table in a "Long" format style
#' 
losses_long = function(loss, val_name = "loss"){
  val_name = enquo(val_name)
  loss_transform = function(loss_table){
    loss_table = as.data.frame(loss_table)
    loss_table = loss_table %>%
      mutate(AccidentYear = as.numeric(rownames(loss_table))) %>%
      gather(DevelopmentLag, !!val_name, - AccidentYear) %>%
      mutate(DevelopmentLag = as.numeric(DevelopmentLag))%>% 
      remove_missing()
    
    loss_table
  }
  res = map(loss, loss_transform)
  # res$incurred_tr = res$incurred_tr %>% remove_missing()
  # res$paid_tr = res$paid_tr %>% remove_missing()
  # if (!is.null(res$ldf)){
  #   res$ldf = res$ldf %>% remove_missing()
  # } 
  res
}

#'Calculate the development factors for a loss triangle for use in Chain Ladder Method
#'
#'@param loss_triangle A loss triangle of claims
#'@result The chain ladder development factors
dev_factors_chain = function(loss_triangle){
  
  dev = rep(1, ncol(loss_triangle)-1)
  
  N=nrow(loss_triangle)
  
  for (j in 2:length(dev)){
    dev[j-1] = 
      sum(loss_triangle[1:(N-j+1),j], na.rm = T) /
      sum(loss_triangle[1:(N-j+1),j-1], na.rm = T)
  }
  dev
}

#'Calculate the development factors for a loss triangle for each cell
#'
#'@param loss_triangle A loss triangle of claims
#'@result The development factors for each cell
dev_factors = function(loss_triangle){
  loss_triangle[,-1]/loss_triangle[,-ncol(loss_triangle)] 
}

#' Project a loss triangle to a loss square given a set of loss development factors
#' 
project_losses = function(loss_triangle, dev_factors){
  losses_projected = loss_triangle
  for (j in 2:ncol(losses_projected)){
    losses_projected[,j] = ifelse(is.na(losses_projected[,j]), 
                                  losses_projected[,j-1] * dev_factors[,j-1],
                                  losses_projected[,j]) 
  }
  losses_projected
}

#'Calculate the incremental claims for a loss triangle for each cell
#'
#'@param loss_triangle A loss triangle of claims
#'@result The development factors for each cell
incremental_claims = function(loss_triangle){
  loss_triangle[,-1]-loss_triangle[,-ncol(loss_triangle)] 
}



#' Converts a set of Simulated LDF factors to Simulated losses
#' 
#' @param sims A N x (px(Q-1)) matrix, N= number of simulations, p= number of accident years and Q= number of development lags. Each row must be in the following format: (p_1, q_2), (p_1, q_2),...(p_1,q_Q),....(p_P, q_Q)
#' @result a list containing loss_square = mean expected loss square and sim_ult_losses = ultimate simulated losses and sims_losses = simulated cumulative losses  by accident year and development lag
ldf_sims_to_losses = function(sims, loss_triangle){
  N= nrow(sims)
  P=nrow(loss_triangle)
  QQ = ncol(loss_triangle)
  if (ncol(sims) != P * (QQ-1)) {
    stop("Number of columns of sims must be p x (Q-1), p=#accident years, Q = number of development lags")
  }
  
  sim_ult_losses = rep(NA, N)
  sim_losses = matrix(NA, nrow=N, ncol= P * QQ)
  for (i in 1:N){
    sim_ldf = matrix(sims[i,], nrow = P, byrow = TRUE)
    
    sim_loss_square = project_losses(loss_triangle,dev_factors = sim_ldf)
    sim_losses[i, ] = c(t(sim_loss_square))
    sim_ult_losses[i] = sum(sim_loss_square[, ncol(sim_loss_square)])
  }
  
  loss_square = colMeans(sim_losses)
  loss_square = matrix(loss_square, nrow = P, byrow = TRUE)
  
  list(loss_square=loss_square, sim_ult_losses=sim_ult_losses, sim_losses=sim_losses)
  
}


#' Converts a set of Simulated ILR factors to Simulated losses
#' 
#' @param sims A N x (pxQ) matrix of ILR factors, N= number of simulations, p= number of accident years and Q= number of development lags. Each row must be in the following format: (p_1, q_2), (p_1, q_2),...(p_1,q_Q),....(p_P, q_Q)
#' @result a list containing loss_square = mean expected loss square and sim_ult_losses = ultimate simulated losses and sims_losses = simulated cumulative losses  by accident year and development lag
ilr_sims_to_losses = function(sims, loss_triangle, premiums){
  N= nrow(sims)
  P=nrow(loss_triangle)
  QQ = ncol(loss_triangle)
  if (ncol(sims) != P * (QQ)) {
    stop("Number of columns of sims must be p x (Q), p=#accident years, Q = number of development lags")
  }
  
  sim_ult_losses = rep(NA, N)
  sim_losses = matrix(NA, nrow=N, ncol= P * QQ)
  for (i in 1:N){
    sim_ilr = matrix(sims[i,], nrow = P, byrow = TRUE) * premiums
    
    sim_loss_square = project_losses_incr_claims(loss_triangle = loss_triangle,incr_claims = sim_ilr)
    sim_losses[i, ] = c(t(sim_loss_square))
    sim_ult_losses[i] = sum(sim_loss_square[, ncol(sim_loss_square)])
  }
  
  loss_square = colMeans(sim_losses)
  loss_square = matrix(loss_square, nrow = P, byrow = TRUE)
  
  list(loss_square=loss_square, sim_ult_losses=sim_ult_losses, sim_losses=sim_losses)
  
}




#' Chain Ladder Method for Projecting the Loss Triangle
#' 
#' @param loss_triangle
#' @param premiums A Vector of premiums. Not required or used for this method, but kept for consistency with other methods
chain_ladder = function(loss_triangle, premiums = NULL, ...){
  dev = dev_factors_chain(loss_triangle)
  
  loss_square = loss_triangle
  
  for (j in 2:ncol(loss_square)){
    loss_square[,j] = ifelse(is.na(loss_square[,j]), 
                             loss_square[,j-1] * dev[j-1],
                             loss_square[,j]) 
  }
  loss_square 
}


#' Chain Ladder Method for Projecting the Loss Triangle
#' 
#' @param loss_triangle
#' @param premiums A Vector of premiums. Not required or used for this method, but kept for consistency with other methods
mack_ladder = function(loss_triangle, premiums = NULL, iter = 1000, ...){
  
  results = MackChainLadder(Triangle = as.triangle(loss_triangle))
  loss_square =results$FullTriangle
  
  # See http://users.stat.ufl.edu/~winner/cases/rocknroll_marathon.pptx for a description of method of moments fit
  sample_mean = sum(loss_square[, ncol(loss_square)])
  sample_sd   = results$Total.Mack.S.E
  
  sdlog = sqrt(log((sample_sd^2 + sample_mean^2)/sample_mean^2))
  meanlog = log(sample_mean) - 0.5 * sdlog^2
  sim_ult_losses = rlnorm(n = iter, meanlog = meanlog, sdlog = sdlog)
  
  list(loss_square =loss_square, sim_ult_losses = sim_ult_losses, summary_results = NULL)
}


#' Bootstrap Chain Ladder using ChainLadder package
boot_ladder = function(loss_triangle, premiums = NULL, iter = 1000, loss_square, ...){
  
  results = BootChainLadder(Triangle = as.triangle(loss_triangle), R=iter)
  loss_square_pred =incr2cum(apply(results$IBNR.Triangles, c(1,2), mean))+getLatestCumulative(loss_triangle)
  loss_square_pred_med =incr2cum(apply(results$IBNR.Triangles, c(1,2), median))+getLatestCumulative(loss_triangle)
  #sim_ult_losses <- apply(t(results$IBNR.ByOrigin[,1,] + as.numeric(getLatestCumulative(loss_triangle))),1,sum)
  sim_ult_losses <- results$IBNR.Totals + sum(getLatestCumulative(loss_triangle))
  #sim_ult_losses = rnorm(n = iter, mean = sum(loss_square[, ncol(loss_square)]), sd = results$Total.Mack.S.E)
  
  errors = error_calcs(premiums = premiums, predicted = loss_square_pred, actual = loss_square)
  
  errors$rmse_mean = rmse_triangle(actual_sq = loss_square/premiums, predicted_sq = loss_square_pred/    premiums,actual_tr = loss_triangle/premiums)
  errors$rmse_med =  rmse_triangle(actual_sq = loss_square/premiums, predicted_sq = loss_square_pred_med/premiums,actual_tr = loss_triangle/premiums)
  
  
  list(loss_square =loss_square_pred, sim_ult_losses = sim_ult_losses, loss_square_pred_med = loss_square_pred_med, 
       errors = errors, summary_results = NULL)
}




#' Project a loss triangle to a loss square given a set of incremental claims
#' 
project_losses_incr_claims = function(loss_triangle, incr_claims){
  if (ncol(loss_triangle) != ncol(incr_claims)){
    stop("loss_triangle and incr_claims must have the same number of columns")
  }
  losses_projected = loss_triangle
  for (j in 2:ncol(losses_projected)){
    losses_projected[,j] = ifelse(is.na(losses_projected[,j]), 
                                  losses_projected[,j-1] + incr_claims[,j],
                                  losses_projected[,j]) 
  }
  losses_projected
}

#' Plots the percentage ranks for each business line
plot_ranks = function(results){
  rank_timing = data.frame(AccidentYear = rep(1988:1997,each=9), DevelopmentLag = rep(2:10, times= 10))
  index = rank_timing$AccidentYear+rank_timing$DevelopmentLag - 1998
  business_lines = names(results)
  ans=list()

  for (bus in business_lines){
    rank_inc = data.frame(step_ahead = index, rank  = c(t(results[[bus]]$rank_vals_inc)))
    rank_cum = data.frame(step_ahead = index, rank  = c(t(results[[bus]]$rank_vals_cum)))
    ans[[bus]]$range_score = c(range_50_inc = mean(rank_inc$rank[rank_inc$step_ahead >=1] >=.25 & rank_inc$rank[rank_inc$step_ahead >=1] <= .75),
                                   range_90_inc = mean(rank_inc$rank[rank_inc$step_ahead >=1] >=.05 & rank_inc$rank[rank_inc$step_ahead >=1] <= .95),
                                   range_50_cum = mean(rank_cum$rank[rank_cum$step_ahead >=1] >=.25 & rank_cum$rank[rank_cum$step_ahead >=1] <= .75),
                                   range_90_cum = mean(rank_cum$rank[rank_cum$step_ahead >=1] >=.05 & rank_cum$rank[rank_cum$step_ahead >=1] <= .95))
  
    # rank_fin = 
    #   ans[[bus]] = list()
    ans[[bus]]$range_inc = rank_inc %>%
      group_by(step_ahead) %>%
      summarize(range_50 = mean(rank >=.25 & rank <= .75),
                range_90 = mean(rank >=.05 & rank <= .95),
                range_99 = mean(rank >=.005 & rank <= .995))
    ans[[bus]]$range_cum = rank_cum %>%
      group_by(step_ahead) %>%
      summarize(range_50 = mean(rank >=.25 & rank <= .75),
                range_90 = mean(rank >=.05 & rank <= .95),
                range_99 = mean(rank >=.005 & rank <= .995))
    
    
    ans[[bus]]$plot_inc = (ggplot(rank_inc[rank_inc$step_ahead>0,]) + 
                             geom_density(aes(rank,fill = as.factor(step_ahead)), alpha=.2) + 
                             labs(fill = "Steps Ahead", title = paste("Percentage Ranks for Incremental Loss Ratios\nBusiness Line:", bus)))
    
    ans[[bus]]$plot_cum = (ggplot(rank_cum[rank_inc$step_ahead>0,]) + 
                             geom_density(aes(rank,fill = as.factor(step_ahead)), alpha=.2) + 
                             labs(fill = "Steps Ahead", title = paste("Percentage Ranks for Cumulative Loss Ratios\nBusiness Line:", bus)))
    
    ans[[bus]]$plot_inc_train = (ggplot(rank_inc[rank_inc$step_ahead<=0,]) + 
                                   geom_density(aes(rank,fill = as.factor(step_ahead)), alpha=.2) + 
                                   labs(fill = "Steps Ahead", title = paste("Training Percentage Ranks for Incremental Loss Ratios\nBusiness Line:", bus)))
    ans[[bus]]$plot_fin = (ggplot(rank_cum[rank_timing$DevelopmentLag == 10,]) + 
                             geom_density(aes(rank), fill = "blue", alpha=.2) + 
                             labs(fill = "Steps Ahead", title = paste("Percentage Ranks for Cumulative Loss Ratios\nBusiness Line:", bus)))
    ans[[bus]]$plot_step_ahead_inc = ggplot(ans[[bus]]$range_inc) +geom_line(aes(step_ahead, range_50, color="50% Range"))+geom_line(aes(step_ahead, range_90, color="90% Range")) + 
      geom_hline(aes(yintercept = 0.9, color="90% Range"),linetype=2)+
      geom_hline(aes(yintercept = 0.5, color="50% Range"),linetype=2) + 
      geom_vline(aes(xintercept = 1), color="green", linetype=2) + 
      labs(title = paste0(bus, ": Incremental Loss Ratios\nProportion of Samples Falling in Inner 50% an 90% Ranges"), y= "% in Range", x="Step Ahead")+
      scale_x_continuous(breaks=min(ans[[bus]]$range_inc$step_ahead): max(ans[[bus]]$range_inc$step_ahead))
    
    ans[[bus]]$plot_step_ahead_cum = ggplot(ans[[bus]]$range_cum %>% filter(step_ahead>=1)) +
      geom_line(aes(step_ahead, range_50, color="50% Range"))+
      geom_line(aes(step_ahead, range_90, color="90% Range")) + 
      geom_hline(aes(yintercept = 0.9, color="90% Range"),linetype=2)+
      geom_hline(aes(yintercept = 0.5, color="50% Range"),linetype=2) + 
      geom_vline(aes(xintercept = 1), color="green", linetype=2) + 
      labs(title = paste0(bus, ": Cumulative Loss Ratios\nProportion of Samples Falling in Inner 50% an 90% Ranges"), y= "% in Range", x="Step Ahead")+
      scale_x_continuous(breaks=min(ans[[bus]]$range_cum$step_ahead): max(ans[[bus]]$range_cum$step_ahead))
    
    ans[[bus]]$plot_step_ahead_inc2 = ggplot(ans[[bus]]$range_inc) +
      geom_line(aes(step_ahead, range_99, color="99% Range"))+
      geom_line(aes(step_ahead, range_90, color="90% Range")) + 
      geom_hline(aes(yintercept = 0.9, color="90% Range"),linetype=2)+
      geom_hline(aes(yintercept = 0.99, color="99% Range"),linetype=2)+
      geom_vline(aes(xintercept = 1), color="green", linetype=2) + 
      labs(title = paste0(bus, ": Incremental Loss Ratios\nProportion of Samples Falling in Inner 50% an 90% Ranges"), y= "% in Range", x="Step Ahead")+
      scale_x_continuous(breaks=min(ans[[bus]]$range_inc$step_ahead): max(ans[[bus]]$range_inc$step_ahead))
    
    ans[[bus]]$plot_step_ahead_cum2 = ggplot(ans[[bus]]$range_cum %>% filter(step_ahead>=1)) +
      geom_line(aes(step_ahead, range_99, color="99% Range"))+
      geom_line(aes(step_ahead, range_90, color="90% Range")) + 
      geom_hline(aes(yintercept = 0.9, color="90% Range"),linetype=2)+
      geom_hline(aes(yintercept = 0.99, color="99% Range"),linetype=2) + 
      geom_vline(aes(xintercept = 1), color="green", linetype=2) + 
      labs(title = paste0(bus, ": Cumulative Loss Ratios\nProportion of Samples Falling in Inner 90% an 99% Ranges"), y= "% in Range", x="Step Ahead")+
      scale_x_continuous(breaks=min(ans[[bus]]$range_cum$step_ahead): max(ans[[bus]]$range_cum$step_ahead))
    
    step_ahead = data.frame(step_ahead = 1:ncol(results[[bus]]$rmse_step_ahead),
                            mean = colMeans(results[[bus]]$rmse_step_ahead, na.rm=T),
                            median = apply(results[[bus]]$rmse_step_ahead,2,median, na.rm=T),
                            quantile_10 = apply(results[[bus]]$rmse_step_ahead,2,quantile, probs=0.1),
                            quantile_90 = apply(results[[bus]]$rmse_step_ahead,2,quantile, probs=0.9))

    ans[[bus]]$plot_step_ahead_rmse = ggplot(step_ahead) +
      geom_line(aes(step_ahead, mean, color="Mean"))+
      geom_line(aes(step_ahead, median, color="Median")) + 
      geom_ribbon(aes(step_ahead, ymin=quantile_10,ymax=quantile_90), fill= "grey", alpha=0.2)+
      labs(title = paste0(bus, ": Step Ahead RMSE"), y= "RMSE on Loss Ratio", x="Step Ahead")+
      scale_x_continuous(breaks=min(ans[[bus]]$range_cum$step_ahead): max(ans[[bus]]$range_cum$step_ahead))
    
    
  }

  #Scores
  ans$rank_scores = map_df(ans, ~.x$range_score)
  ans$rank_scores$score = names(ans[[1]]$range_score)
  
  
  ans
}


#' Overdispersed Poisson Method for Projecting the Loss Triangle
#' 
#' @param loss_triangle
#' @param premiums A Vector of premiums. 
od_poisson = function(loss_triangle, premiums = NULL, loss_square = NULL, overdispersion = TRUE,...){
  incr_claims = incremental_claims(loss_triangle) %>% as.data.frame()
  incr_claims = cbind("1" = loss_triangle[,1], incr_claims, premiums = premiums)
  
  incr_claims = incr_claims %>%
    mutate(AccidentYear = factor(as.numeric(rownames(incr_claims)))) %>%
    gather(DevelopmentLag, y, - AccidentYear, -premiums) %>%
    mutate(DevelopmentLag = as.factor(DevelopmentLag)) 
  
  new_data = incr_claims[is.na(incr_claims$y),]
  
  
  incr_claims = incr_claims %>%
    remove_missing(na.rm = TRUE) %>%
    mutate(y = pmax(y, 0))

  
  if (!overdispersion) {
    # + offset(log(premiums)
    odp_model = stats::glm(formula = y ~ AccidentYear + DevelopmentLag + offset(log(premiums)), 
                           family = stats::poisson(), data = incr_claims)
  } else {
    odp_model = stats::glm(formula = y ~ AccidentYear + DevelopmentLag + offset(log(premiums)), 
                           family = quasipoisson, data = incr_claims) 
  } 
  
  new_data$y = predict(odp_model, newdata = new_data, type = "response")
  preds = rbind(incr_claims, new_data) %>%
    mutate(DevelopmentLag = as.numeric(as.character(DevelopmentLag)),
           AccidentYear = as.numeric(as.character(AccidentYear))) %>%
    spread(DevelopmentLag, y) %>%
    arrange(AccidentYear) %>%
    select(-premiums, - AccidentYear)

  losses_projected = project_losses_incr_claims(loss_triangle, incr_claims = preds)
  losses_projected
  
}




#' Overdispersed Poisson Method using ChainLadder Package for Projecting the Loss Triangle
#' 
#' @param loss_triangle
#' @param premiums A Vector of premiums. 
od_poisson_glm = function(loss_triangle, premiums = NULL, iter = 1000, ...){
  
  
  # Adjusts triangle to ensure no negative losses
  loss_triangle_adj = loss_triangle
  for (i in 2:ncol(loss_triangle_adj)){
    loss_triangle_adj[,i] = pmax(loss_triangle_adj[,i], loss_triangle_adj[,i-1])
  }
  results = glmReserve(as.triangle(loss_triangle_adj),  mse.method = "formula", nsim = iter,var.power = 1, cum = TRUE)
  loss_square =results$FullTriangle
  sim_ult_losses = rnorm(n = iter, mean = sum(loss_square[, ncol(loss_square)]), sd = results$summary$S.E[length(results$summary$S.E)])
  
  # latest = sum(diag(loss_triangle[,ncol(loss_triangle):1])) # calculatest the latest cumulative paid losses
  # pr <- as.data.frame(results$sims.reserve.pred) #predicted reserves
  # sim_ult_losses = rowSums(pr) + latest # total cumulative ultimate losses

  list(loss_square =loss_square, sim_ult_losses = sim_ult_losses)
}







#' Calculates the Root Mean Square Error between two vectors or matrices of the same size
#' 
#' @param x1, x2 Vectors or matrices x1, x2
#' @para, na.rm Remove na values
#' @return The RMSE 
rmse = function(x1, x2,na.rm = F){
  sqrt(mean(c((x1-x2)^2), na.rm = na.rm))
}


#' Calculates the Root Mean Square Error of a Predicted set of Rates from a Triangle
#' 
#' @param actual_sq The actual losses in a square format
#' @param predicted_sq The predicted losses in a square format
#' @param actual_tr The actual triangle used to calculate the predicted values
#' 
#' @return The RMSE of the prediction.  The actual_tr is required.  The RMSE is only calculated for those cells for which there is not a value in actual_tr
rmse_triangle = function(actual_sq, predicted_sq, actual_tr){
  sqrt(mean(c((actual_sq - predicted_sq)^2)[is.na(c(actual_tr))], na.rm = T))
}

#' Loss Projection using Gaussian Process to Project CUMULATIVE Loss Development Factors
#' 
#' @param loss_triangle A loss triangle of data, it must include rownames (AccidentYear) and column names (DevelopmentYear)
#' @param premiums A Vector of premiums. Not required or used for this method, but kept for consistency with other methods
gp_ladder_cum = function(loss_triangle, formula = ~ AccidentYear + log(DevelopmentYear), premiums = NULL){
  eps = 0.05
  N=nrow(loss_triangle)
  if (is.null(rownames(loss_triangle))){
    stop("The rownames of the loss triangle must be the accident year")
  }
  if (is.null(colnames(loss_triangle))){
    stop("The colnames of the loss triangle must be the development year")
  }
  
  AccidentYear = as.numeric(rownames(loss_triangle))
  DevelopmentYear = as.numeric(colnames(loss_triangle)[-1])
  no_acc = length(AccidentYear); no_dev = length(DevelopmentYear)
  AccidentYear = rep(AccidentYear, each = no_dev)
  DevelopmentYear = rep(DevelopmentYear, times = no_acc)
  
  cum_dev = dev_factors(loss_triangle)
  cum_dev = t(apply(cum_dev,1,cumprod))
  
  dev = cum_dev %>%
    as_tibble() %>%
    mutate(AccidentYear = as.integer(rownames(loss_triangle))) %>%
    gather(DevelopmentYear, dev, -AccidentYear) %>%
    filter(!is.na(dev)) %>%
    mutate(dev_log = log(dev),
           DevelopmentYear = as.integer(DevelopmentYear))
  
  mean_dev_log = mean(dev$dev_log)
  sd_dev_log = sd(dev$dev_log)
  dev$dev_log = (dev$dev_log - mean_dev_log)/sd_dev_log
  
  
  ######################################################
  #### Fit GP
  
  model_nug <- km(formula = formula, 
                  design = dev[, c("AccidentYear", "DevelopmentYear")], 
                  response = dev[, c("dev_log")],
                  nugget.estim=TRUE,
                  covtype="gauss",
                  optim.method="gen",
                  # the "control" parameters below handle speed versus risk of
                  # converging to local minima.  See "rgenoud" package for details
                  control=list(max.generations=100,pop.size=100,wait.generations=8,
                               solution.tolerance=1e-5))
  model_nug
  
  nug <- model_nug@covariance@nugget
  model <- km(formula = formula, 
              design = model_nug@X, response = model_nug@y,
              noise.var = rep(nug,model_nug@n),
              coef.trend = model_nug@trend.coef,  # re-use obtained hyperparameters
              coef.cov = model_nug@covariance@range.val,
              coef.var = model_nug@covariance@sd2, 
              covtype = model_nug@covariance@name)
  print(model)
  
  # Do predictions on upper and lower triangles
  preds = tibble(AccidentYear ,
                 DevelopmentYear)
  dev_log_preds <- predict(model, newdata=preds,
                           cov.compute=TRUE,
                           se.compute=TRUE,type="UK")
  
  # Unscale and take exponent
  preds$dev_preds = exp(dev_log_preds$mean * sd_dev_log + mean_dev_log)
  
  # Join with original data
  preds = preds %>%
    left_join(dev)
  
  preds
  
  # Plot Actual vs. Predicted Factors
  print(ggplot(cbind(preds, intercept=1)) + 
    geom_point(aes(DevelopmentYear, dev, color="Actual"))+
    geom_line(aes(DevelopmentYear, dev_preds, color="Predicted")) + 
    facet_wrap(~p0("Acc. Yr ",AccidentYear)) + 
    geom_hline(aes(yintercept = intercept), size = 0.25, linetype = 2) + 
    labs(title = p0("GP Fit to Log Loss Development Factors\nformula = ", format(formula) )))
  
  
  ldf_sq =   preds %>%
    select(AccidentYear, DevelopmentYear, dev_preds) %>%
    spread(DevelopmentYear, dev_preds) %>%
    as.matrix()
  rownames(ldf_sq) = ldf_sq[, "AccidentYear"]
  ldf_sq = ldf_sq[,-1]
  ldf_sq
  
  loss_square = loss_triangle
  for (j in 2:ncol(loss_square)){
    if (j==2){
      loss_square[,j] = ifelse(is.na(loss_triangle[,j]),loss_square[,j-1] * ldf_sq[,j-1], loss_triangle[,j])
    } else {
      loss_square[,j] = ifelse(is.na(loss_triangle[,j]),loss_square[,j-1] * ldf_sq[,j-1]/ldf_sq[,j-2], loss_triangle[,j])
    }
  } 
   
  loss_square 
}


#' Loss Projection using Gaussian Process to Project ANNUAL Loss Development Factors
#' 
#' @param loss_triangle A loss triangle of data, it must include rownames (AccidentYear) and column names (DevelopmentYear)
#' @param premiums A Vector of premiums. Not required or used for this method, but kept for consistency with other methods
#' @param loss_square Not required or used for this method, but kept for consistency with other methods
gp_ladder = function(loss_triangle, formula = ~ AccidentYear + DevelopmentLag,
                     eps = 0.0005, plot_graphs = TRUE, premiums=NULL, loss_square = NULL, iter=1000){
  N=nrow(loss_triangle)
  if (is.null(rownames(loss_triangle))){
    stop("The rownames of the loss triangle must be the accident year")
  }
  if (is.null(colnames(loss_triangle))){
    stop("The colnames of the loss triangle must be the development year")
  }
  
  AccidentYear = as.numeric(rownames(loss_triangle))
  DevelopmentLag = as.numeric(colnames(loss_triangle)[-1])
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
           # AccidentYear_scale = (AccidentYear - AccidentYear_mean)/AccidentYear_sd,
           # DevelopmentLag_scale = (DevelopmentLag - DevelopmentLag_mean)/DevelopmentLag_sd)
  
  mean_dev_log = mean(dev$dev_log)
  sd_dev_log = max(eps^2,sd(dev$dev_log))
  dev$dev_log = (dev$dev_log - mean_dev_log)/sd_dev_log
  
  
  ######################################################
  #### Fit GP
  
  model_nug <- km(formula = formula, 
                  design = dev[, c("AccidentYear", "DevelopmentLag")], 
                  response = dev[, c("dev_log")],
                  nugget.estim=TRUE,
                  covtype="gauss",
                  optim.method="gen",
                  # the "control" parameters below handle speed versus risk of
                  # converging to local minima.  See "rgenoud" package for details
                  control=list(max.generations=100,pop.size=100,wait.generations=8,
                               solution.tolerance=1e-5, print.level=0))
  model_nug
  
  nug <- model_nug@covariance@nugget
  model <- km(formula = formula, 
              design = model_nug@X, response = model_nug@y,
              noise.var = rep(nug,model_nug@n),
              coef.trend = model_nug@trend.coef,  # re-use obtained hyperparameters
              coef.cov = model_nug@covariance@range.val,
              coef.var = model_nug@covariance@sd2, 
              covtype = model_nug@covariance@name)
  print(model)
  
  # Do predictions on upper and lower triangles
  preds = tibble(AccidentYear=AccidentYear, DevelopmentLag=DevelopmentLag)
  #AccidentYear_scale=AccidentYear_scale, DevelopmentLag_scale=DevelopmentLag_scale)
  dev_log_preds <- predict(model, newdata=preds,
                           cov.compute=TRUE,
                           se.compute=TRUE,type="UK")
  
  sims = simulate(object = model, newdata=preds, nsim = iter,  nugget.sim=nug/100, cond = TRUE)
  
  
  # Unscale and take exponent
  sims = exp(sims * sd_dev_log + mean_dev_log)+1-eps
  
  
  # Calc simulated  losses
  results = ldf_sims_to_losses(sims = sims, loss_triangle = loss_triangle)
  
  #preds$dev_preds = exp(dev_log_preds$mean * sd_dev_log + mean_dev_log)+1-eps
  preds$dev_preds = c(t(dev_factors(results$loss_square)))
  
      
  # Join with original data
  
  
  # Plot Actual vs. Predicted Factors
  if (plot_graphs){
    preds = preds %>%
      left_join(dev, by = c("AccidentYear", "DevelopmentLag"))
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
  
  
  # ldf_sq =   preds %>%
  #   select(AccidentYear, DevelopmentLag, dev_preds) %>%
  #   spread(DevelopmentLag, dev_preds) %>%
  #   as.matrix()
  # rownames(ldf_sq) = ldf_sq[, "AccidentYear"]
  # ldf_sq = ldf_sq[,-1]
  # ldf_sq
  # 
  #loss_square = project_losses(loss_triangle = loss_triangle,dev_factors = ldf_sq)

  list(loss_square =results$loss_square, sim_ult_losses = results$sim_ult_losses) 
}



#' Loss Projection using Gaussian Process to Project ANNUAL Incremental Loss Factors, using DiceKriging Package
#' 
#' @param loss_triangle A loss triangle of data, it must include rownames (AccidentYear) and column names (DevelopmentYear)
#' @param premiums A Vector of premiums. Not required or used for this method, but kept for consistency with other methods
#' @param loss_square Not required or used for this method, but kept for consistency with other methods
gp_ilr = function(loss_triangle, formula = ~ AccidentYear + DevelopmentLag,
                     eps = 0.0005, plot_graphs = TRUE, premiums=NULL, loss_square = NULL, iter=1000){
  
  N=nrow(loss_triangle)
  if (is.null(rownames(loss_triangle))){
    stop("The rownames of the loss triangle must be the accident year")
  }
  if (is.null(colnames(loss_triangle))){
    stop("The colnames of the loss triangle must be the development year")
  }
  
  AccidentYear = as.numeric(rownames(loss_triangle))
  DevelopmentLag = as.numeric(colnames(loss_triangle))
  no_acc = length(AccidentYear); no_dev = length(DevelopmentLag)
  AccidentYear = rep(AccidentYear, each = no_dev)
  DevelopmentLag = rep(DevelopmentLag, times = no_acc)
  
  scaler_accident_year  = scaler(AccidentYear)
  AccidentYear_scaled = c(scaler_accident_year(AccidentYear))
  scaler_development_lag  = scaler(DevelopmentLag)
  DevelopmentLag_scaled = c(scaler_accident_year(DevelopmentLag))
  
  
  ilr = incremental_claims(cbind(0,loss_triangle/premiums)) %>%
    as_tibble() %>%
    mutate(AccidentYear = as.integer(rownames(loss_triangle))) %>%
    gather(DevelopmentLag, ilr, -AccidentYear) %>%
    #filter(!is.na(ilr)) %>%
    mutate(DevelopmentLag = as.integer(DevelopmentLag),
           AccidentYear_scaled = c(scaler_accident_year(AccidentYear)),
           DevelopmentLag_scaled = c(scaler_development_lag(DevelopmentLag))) %>%
    arrange(AccidentYear, DevelopmentLag)
  
  scaler_ilr = scaler(ilr$ilr[!is.na(ilr$ilr)])
  ilr$ilr_scaled = c(scaler_ilr(ilr$ilr))
  
  
  preds = ilr 
  ilr = ilr %>%
    filter(!is.na(ilr))
  design_ilr = ilr[, c("AccidentYear", "DevelopmentLag")]
  design_preds = preds[, c("AccidentYear", "DevelopmentLag")]
  
  ######################################################
  #### Fit GP
  
  model_nug <- km(formula = formula, 
                  design = design_ilr, 
                  response = ilr[, c("ilr_scaled")],
                  nugget.estim=TRUE,
                  covtype="gauss",
                  optim.method="gen",
                  # the "control" parameters below handle speed versus risk of
                  # converging to local minima.  See "rgenoud" package for details
                  control=list(max.generations=100,pop.size=100,wait.generations=8,
                               solution.tolerance=1e-5, print.level=0))
  model_nug
  
  nug <- model_nug@covariance@nugget
  model <- km(formula = formula, 
              design = model_nug@X, response = model_nug@y,
              noise.var = rep(nug,model_nug@n),
              coef.trend = model_nug@trend.coef,  # re-use obtained hyperparameters
              coef.cov = model_nug@covariance@range.val,
              coef.var = model_nug@covariance@sd2, 
              covtype = model_nug@covariance@name)
  print(model)
  
  # Do predictions on upper and lower triangles
  ilr_preds <- predict(model, newdata=design_preds,
                           cov.compute=TRUE,
                           se.compute=TRUE,type="UK")
  
  # Simulate results
  sims = simulate(object = model, newdata=design_preds, nsim = iter,  nugget.sim=nug/100, cond = TRUE)
  
  # unscale
  sims = apply(sims,2, FUN = function(x) c(scaler_ilr(x, unscale = T)))
  
  # Calc simulated  losses
  results = ilr_sims_to_losses(sims = sims, loss_triangle = loss_triangle, premiums = premiums)
  
  
  #preds$dev_preds = exp(dev_log_preds$mean * sd_dev_log + mean_dev_log)+1-eps
  preds$ilr_preds = c(t(incremental_claims(cbind(0,results$loss_square/premiums))))
  
  #preds$ilr_preds = c(scaler_ilr(ilr_preds$mean, unscale = T))


  
  # Plot Actual vs. Predicted Factors
  if (plot_graphs){
    print(ggplot(cbind(preds, intercept=0)) + 
            theme_bw() +
            geom_point(aes(DevelopmentLag, ilr, color="Actual"))+
            geom_line(aes(DevelopmentLag, ilr_preds, color="Predicted")) + 
            facet_wrap(~p0("Acc. Yr ",AccidentYear)) + 
            geom_hline(aes(yintercept = intercept), size = 0.25, linetype = 2) + 
            labs(title = p0("GP Fit to Log Loss Development Factors\nformula = ", format(formula) )))
    
    print(ggplot(cbind(preds, intercept=0)) + 
            theme_bw() +
            geom_point(aes(AccidentYear, ilr, color="Actual"))+
            geom_line(aes(AccidentYear, ilr_preds, color="Predicted")) + 
            facet_wrap(~p0("Dev. Yr ",sprintf("%02d", DevelopmentLag))) + 
            geom_hline(aes(yintercept = intercept), size = 0.25, linetype = 2) + 
            labs(title = p0("GP Fit to Log Loss Development Factors\nformula = ", format(formula) )))
  }
  
  
  
  list(loss_square =results$loss_square, sim_ult_losses = results$sim_ult_losses) 
}




#' Plots actual to projected square loss developments
#' 
#' @param actual_sq, predicted_sq Actual and predicted square losses with appropriate rownames (AccidentYear) and column (DevelopmentLag) names
plot_ae = function(actual_sq, predicted_sq, predicted_sq2 = NULL, title = "NULL"){
  
  
  if (is.null(title))
    title= "Actual to Predicted Loss Development"
  
  AccidentYear_actual = rownames(actual_sq)
  AccidentYear_predicted = rownames(predicted_sq)
  
  actual_sq = actual_sq %>%
    as_tibble() %>%
    mutate(AccidentYear = AccidentYear_actual) %>%
    gather(DevelopmentLag, Actual, -AccidentYear) %>%
    mutate(AccidentYear = as.numeric(AccidentYear),
           DevelopmentLag = as.numeric(DevelopmentLag))
  
  actual_sq
  
  predicted_sq = predicted_sq %>%
    as_tibble() %>%
    mutate(AccidentYear = AccidentYear_predicted) %>%
    gather(DevelopmentLag, Predicted, -AccidentYear) %>%
    mutate(AccidentYear = as.numeric(AccidentYear),
           DevelopmentLag = as.numeric(DevelopmentLag))
  
  predicted_sq
  
  ae = left_join(actual_sq, predicted_sq)

  if (!is.null(predicted_sq2)){
    predicted_sq2 = predicted_sq2 %>%
      as_tibble() %>%
      mutate(AccidentYear = AccidentYear_predicted) %>%
      gather(DevelopmentLag, Predicted2, -AccidentYear) %>%
      mutate(AccidentYear = as.numeric(AccidentYear),
             DevelopmentLag = as.numeric(DevelopmentLag))
    
    predicted_sq2
    ae = left_join(ae, predicted_sq2)
    
  }
  
  gg  =ggplot(ae) +
    theme_bw() +
    geom_line(aes(DevelopmentLag, Actual, color = "Actual"))+
    geom_line(aes(DevelopmentLag, Predicted, color = "Predicted")) + 
    facet_wrap(~AccidentYear, scales = "free" ) + 
    labs(title = title, x= "Development Lag", y= "Developed Losses")
  
  if (!is.null(predicted_sq2)){
    gg = gg + 
      geom_line(aes(DevelopmentLag, Predicted2, color = "Predicted2")) 
  } 
  
  gg
}


#' Predicts a Loss Square from a Loss Triangle using Predicted Development Factors
#' 
#' @param loss_triangle Triangle of losses
#' @param dev_factors A data frame of loss development factors
#' @param AccidentYear_name,DevelopmentYear_name,dev_name The column names in dev_factors designating the various required fields
#' @return A square of losses
make_loss_square = function(loss_triangle, dev_factors,
                            AccidentYear_name = "AccidentYear", 
                            DevelopmentLag_name = "DevelopmentLag",
                            dev_name = "pred_dev"){ 
  
  
  ldf_sq =   dev_factors %>%
    select_(AccidentYear_name, DevelopmentLag_name, dev_name) %>%
    spread_(DevelopmentLag_name, dev_name) %>%
    as.matrix()
  rownames(ldf_sq) = ldf_sq[, AccidentYear_name]
  ldf_sq = ldf_sq[,-1]
  ldf_sq
  
  loss_square = loss_triangle
  
  for (j in 2:ncol(loss_square)){
    loss_square[,j] = ifelse(is.na(loss_triangle[,j]),loss_square[,j-1] * ldf_sq[,j-1], loss_triangle[,j])
  }
  loss_square
}

#' Selects all the companies with full loss squares for a given business line
#' 
#' @param data loss data input
#' @param business_line Name of business line
#' @result Subset of the data that are produce full squares for a particular business line.
#' @examples 
#' data = load_losses()
#' (ans = get_full_squares(data))
get_full_squares = function(data, business_line="prodliab"){
  codes = unique(data$GRCODE[data$line==business_line])
  codes = codes[sapply(codes, function(x) sum(!is.na(loss_table(data = data, business_line = business_line, company_code = x)$paid_tr[,1]))) == 10]
  
  data = data %>% 
    filter(line == business_line, GRCODE %in% codes) 
  
  cell_count = data %>% 
    filter(line == business_line, GRCODE %in% codes) %>%
    group_by(GRCODE) %>%
    summarize( count = n())
  if (any(cell_count$count != 100)){
    stop("There are incompete squares")
  }
  
  data 
}

#' Get the premiums by accident year for selected business lines and company codes
#' 
#' @param data The input loss data
#' @param business_line The selected business line (If NULL, the all lines are selected)
#' @param company_code The selected company code (If NULL, the all company codes are selected)
#' @return the premiums for the selected data
get_premium = function(data, business_line = NULL, company_code=NULL) {
  
  if (!is.null(business_line)){
    data = data %>% filter(line == business_line)
  }
  if (!is.null(company_code)){
    data = data %>% filter(GRCODE == company_code)
  }
  
  prems = data %>% 
    group_by(AccidentYear) %>%
    summarize(premium = min(EarnedPremNet), premium_max = max(EarnedPremNet))
  
  if (any(prems$premium != prems$premium_max)){
    warning("A loss table has premiums that vary by development year")
  }
  return(prems$premium)
  
}


#' Calculates the approximate Continuous Ranked Probability Score
#' 
#' @param r Actual measured value
#' @param r_distribution Simulated distribution of r
#' @result CRPS Value
crps_fun = function(r, r_distribution){
  r_mean = mean(r_distribution)
  r_sd = sd(r_distribution)
  z=(r-r_mean)/r_sd
  ans = -r_sd * (1/sqrt(pi) -2 * dnorm(z) - z * (2 * pnorm(q = z)-1) )
  ans
}

#' Calculates the approximate Negative Log Probability Density
#' 
#' @param r Actual measured value
#' @param r_distribution Simulated distribution of r
#' @result NLPD Value
nlpd_fun = function(r, r_distribution){
  r_mean = mean(r_distribution)
  r_sd = sd(r_distribution)
  z=(r-r_mean)/r_sd
  ans = z^2 + log(r_sd^2)
  ans
}


#' Calculates the percentage of actuals falling with a percentage interval of the data
#' 
#' @param r Actual measured value
#' @param r_distribution Simulated distribution of r
#' @param percentage The cover interval percentage
#' @result NLPD Value
coverage = function(r , r_distribution , percentage = 0.9){
  
  sum((r >= quantile(r_distribution, probs = (1-percentage)/2)) & (r <= quantile(r_distribution, probs = percentage + (1-percentage)/2)))
  
}


#' Calculate the RMSE for all Companies that have full squares for a particular business line
#' 
#' @param projection_fun The projection function that returns a square for a given loss triangle
#' @param data the data to be used for training
#' @param business_line The selected business line
#' @param ... any additional parameters to be given to projection_fun
#' @return a tibble containing the RMSE for each selected company for the projected loss square, and the RMSE for the projected cumulative loss ratio.  Cumulative loss ratio = projected square / annual premiums  ("EarnedPremNet")
rmse_business_line = function(projection_fun=chain_ladder, data, business_line="prodliab", print_res = TRUE,...){
  data = get_full_squares(data, business_line)
  codes = unique(data$GRCODE)
  N_tests = length(codes)
  results = tibble(code = codes, rmse = rep(NA, N_tests), rmse_loss_ratio  = NA)
  i=1
  for  (code in codes){
    print(paste("Working on... ", i, " of ", N_tests, " ", business_line, code))
    loss  =  loss_table(data = data, business_line = business_line, company_code = code)
    premiums = get_premium(data = data, company_code = code)
    
    pred_sq = try(projection_fun(loss$paid_tr, premiums=premiums,...))
    if (is.error(pred_sq)){
      results$rmse[i] = NA
      results$rmse_loss_ratio[i] = NA
      print(paste("ERROR in Calculation", results[i,]))
    } else {
      
      results$rmse[i] = rmse_triangle(actual_sq = loss$paid_sq, predicted_sq = pred_sq, actual_tr = loss$paid_tr)
      results$rmse_loss_ratio[i] = rmse_triangle(actual_sq = loss$paid_sq/ premiums, predicted_sq = pred_sq/ premiums, actual_tr = loss$paid_tr)
      if (print_res) print(results[i,])
    }
    i=i+1
  }
  results
}




#########################################################################
##  Comparison of Business Lines including Percentage Rank Analysis
#########################################################################


#' Comparison of business lines with Percentage Rank Analysis
#' 
#' @param projection_fun Function to produce RMSE and Percentage Ranks
#' @param business_line Business line to be analyzed
#' @param print_res If TRUE, then interim results are printed
#' @param rank_size Number of cells in square that need to be predicted for rank val tests
#' @results The RMSE calculations, the percentage ranks per cell and the percentage ranks for cumulative claims
business_line_analysis = function(projection_fun=chain_ladder, business_line="prodliab", print_res = TRUE, rank_size = 90, ...){
  #data = load_losses_subset()
  codes = unlist(unique(data %>% filter(line==business_line) %>% select(GRCODE)))
  N_tests = length(codes)
  results = tibble(code = codes, rmse = rep(NA, N_tests), rmse_loss_ratio  = NA, rmse_last_col = NA)
  rmse_step_ahead = matrix(NA,length(codes), 9)
  colnames(rmse_step_ahead)=1:9
  rank_vals_inc = matrix(NA, N_tests, rank_size)
  rank_vals_cum = matrix(NA, N_tests, rank_size)
  i=1
  for  (code in codes){
    print(paste("Working on... ", i, " of ", N_tests, " ", business_line, code))
    loss = loss_table(data, business_line = business_line, company_code = code)
    premiums = get_premium(data, business_line = business_line, company_code = code)
    preds = try(projection_fun(loss_triangle = loss$paid_tr, premiums=premiums,loss_square=loss$paid_sq, ...))
    if (is.error(preds)){
      results$rmse[i] = NA
      results$rmse_loss_ratio[i] = NA
      rmse_step_ahead[i,]=NA
  
      print(paste("ERROR in Calculation", results[i,]))
    } else {
      
      results$rmse[i] = rmse_triangle(actual_sq = loss$paid_sq, predicted_sq = preds$loss_square, actual_tr = loss$paid_tr)
      results$rmse_loss_ratio[i] = preds$errors$rmse_mean #rmse_triangle(actual_sq = loss$paid_sq/ premiums, predicted_sq = preds$loss_square/ premiums, actual_tr = loss$paid_tr)
      results$rmse_loss_ratio_med[i] = preds$errors$rmse_med #rmse_triangle(actual_sq = loss$paid_sq/ premiums, predicted_sq = preds$loss_square_med/ premiums, actual_tr = loss$paid_tr)
      results$rmse_last_col[i] = preds$errors$rmse_last_col
      
      rmse_step_ahead[i,] = preds$errors$rmse_step_ahead
      
      rank_vals_inc[i,] = preds$rank_val_incremental
      rank_vals_cum[i,] = preds$rank_val_cumulative
      
      if (print_res) print(results[i,])
    }
    i=i+1
  }
  list(results=results,rmse_step_ahead = rmse_step_ahead, rank_vals_inc=rank_vals_inc, rank_vals_cum=rank_vals_cum)
}



#' Calculation of various error measures
#' 
#' @param premiums Vector of premiums
#' @param actual Actual loss square
#' @param predicted predicted loss square
#' @return RMSE for (i) final column, (ii) for each step ahead

error_calcs = function(premiums, actual, predicted){
  
  actual = actual / premiums
  predicted = predicted / premiums
  
  N = nrow(actual)
  if (ncol(actual) != N | N != nrow(predicted) | N != ncol(predicted) ){
    stop("actual and predicted must be squares of the same size")
  }
  
  errors_sq = (actual - predicted)^2
  
  # Loss Ratio RMSE on whole trianlge
  
  # RMSE of final column
  rmse_last_col = sqrt(mean(errors_sq[, N]))
  
  # RMSE for each step ahead
  rmse_step_ahead = rep(0, N-1)
  
  diag_indicator <- row(errors_sq) + col(errors_sq) - N - 1
  
  # diag_vecs_actual = split(actual, diag_indicator)
  # diag_vecs_predicted = split(predicted, diag_indicator)
  diag_vecs = split(errors_sq,diag_indicator)
  for (n in 1:(N-1)){
    d= as.character(n)
    rmse_step_ahead[n] = sqrt(mean(diag_vecs[[d]]))
  }
  names(rmse_step_ahead) = 1:(N-1)
  return(list(rmse_last_col = rmse_last_col, rmse_step_ahead = rmse_step_ahead))
  
}



#' Summarize the RMSE on Loss Ratio
#' 
#' @param results the list of outputs of the function business_line_analysis
summarize_rmse = function(results){
  summary_rmse = rbind(
    map_df(results, ~ median(.x$results$rmse_loss_ratio, na.rm = T)),
    map_df(results, ~ mean(.x$results$rmse_loss_ratio, na.rm = T)),
    map_df(results, ~ median(.x$results$rmse_loss_ratio_med, na.rm = T)),
    map_df(results, ~ mean(.x$results$rmse_loss_ratio_med, na.rm = T)),
    map_df(results, ~ median(.x$results$rmse_last_col, na.rm = T)),
    map_df(results, ~ mean(.x$results$rmse_last_col, na.rm = T))
  ) %>% as.data.frame()
  rownames(summary_rmse) = c("Median RMSE of Mean LR", "Mean RMSE of Mean LR",  
                             "Median RMSE of Med. LR", "Mean RMSE of Med. LR",  
                             "Median RMSE of of Last Col LR", "Mean RMSE of Last Col LR")
  summary_rmse = apply(summary_rmse,2, round, digits = 3)
  summary_rmse
} 


#########################################################################
##  Examples of the various utility and modeling functions
#########################################################################

general_gp_functions_wrapper = function(){
  
  #  Load Data
  data = load_losses()
  
  # List company codes for a particular line of business
  
  # Select a single line of business from the data
  select_line(data, "comauto")
  
  # List company codes applicable for a specified line of business
  (company_medmal = unique(select_line(data, "comauto")$GRCODE))
  (company_medmal = unique(select_line(data, "prodliab")$GRCODE))
  
  
  # Develop Square and Triangular loss development tables
  loss_table(data, "comauto", company_code = 266)
  loss_table(data, "otherliab", company_code = 266)
  loss_table(data, "otherliab", company_code = 337)

  # Calculate square and triangle tables
  loss  =  loss_table(data, "prodliab", company_code = 78) #78 86   248     388   620   667   715 
  
  # e.g. Paid loss triangle
  loss$paid_tr
  
  #Chain ladder development factors:  one factor per development year
  dev_factors_chain(loss$paid_tr)

  #Cell level development factors
  dev_factors(loss$paid_tr)
  
  # Chain Ladder Projection
  pred_sq = chain_ladder(loss$paid_tr)
  pred_sq
  rmse_triangle(actual_sq = loss$paid_sq, predicted_sq = pred_sq, actual_tr = loss$paid_tr)
  plot_ae(loss$paid_sq, pred_sq )
  
  # Chain Ladder for companies in a given business line
  res_chain_ladder = rmse_business_line(projection_fun = chain_ladder, data=data, business_line = "prodliab")
  mean(res_chain_ladder$rmse); mean(res_chain_ladder$rmse_loss_ratio)
  median(res_chain_ladder$rmse); median(res_chain_ladder$rmse_loss_ratio)
  
  
 
  # GP Ladder Projection
  loss  =  loss_table(data, "comauto", company_code = 38997) #78 86   248     388   620   667   715 
  prem  =  get_premium(data, "comauto", company_code = 38997) #78 86   248     388   620   667   715 
  pred_sq = gp_ladder(loss_triangle = loss$paid_tr, 
                      formula = ~ log(DevelopmentLag) + (AccidentYear), eps=0.0001)
  pred_sq
  rmse_triangle(actual_sq = loss$paid_sq, predicted_sq = pred_sq, actual_tr = loss$paid_tr)
  plot_ae(loss$paid_sq, pred_sq )
  
  
  # GP  Ladder for companies in a given business line
  res = rmse_business_line(projection_fun = gp_ladder, data=data, business_line = "prodliab", 
                           formula = ~ log(DevelopmentLag) + (AccidentYear), eps=0.1, plot_graphs= FALSE)
  mean(res$rmse); mean(res$rmse_loss_ratio)
  median(res$rmse); median(res$rmse_loss_ratio)
  res_combined = cbind(res_chain_ladder, res)
  res_combined$chain_minus_gp = res_combined[,3] - res_combined[,6]
  plot(res_combined$chain_minus_gp[abs(res_combined$chain_minus_gp)<1])
  

  # Stan GP  Ladder for companies in a given business line
  
  res_stan_ladder = rmse_business_line(projection_fun = stan_ladder, data=data)
  mean(res_stan_ladder$rmse); mean(res_stan_ladder$rmse_loss_ratio)
  median(res_stan_ladder$rmse); median(res_stan_ladder$rmse_loss_ratio)
  
  
  # GP Ladder Fitted to CUMULATIVE LDFs
  # This does not work very well
  pred_sq = gp_ladder_cum(loss_triangle = loss$paid_tr, formula = ~ log(DevelopmentLag) + AccidentYear)
  pred_sq
  rmse_triangle(actual_sq = loss$paid_sq, predicted_sq = pred_sq, actual_tr = loss$paid_tr)
  plot_ae(loss$paid_sq, pred_sq )
     
}

#' Check for an error
#' 
is.error <- function(x) inherits(x, "try-error")

#' Compare the results of different methods accross all business lines
#' 
#' Business lines include chain ladder, GP using Dice Kriging, GP using Stan HMC
#' 
#' @param business_lines A string vector containing the code names for each of the business lines to be investigated.  Choices are: c("comauto", "medmal", "othliab","ppauto","prodliab", "wkcomp")
#' @return A summary of the error rates for each algorithm
compare_methods = function(business_lines = c("comauto", "medmal", "othliab","ppauto","prodliab", "wkcomp"),
                           methods = c("chain_ladder","gp_ladder","stan_ladder", "stan_loss_ratio_incr", "stan_ladder_se", 
                                       "stan_ladder_compound", "stan_compound_log", "stan_comp_het", "stan_inc_het", "stan_inc_vb")){
 
  
  #  Load Data
  data = load_losses()
  
  # Set business lines to analyze

  #business_lines = c("prodliab")
  all_res = list()
  final_res = NULL
  for (bus in business_lines){
    all_res[[bus]] = list()
    # Method 1.1:  Chain Ladder 
    if ("chain_ladder" %in% methods){
      all_res[[bus]][["chain_ladder"]] = rmse_business_line(projection_fun = chain_ladder, data=data, business_line = bus)
    }
    
    # Method 1.2:  Chain Ladder 
    if ("od_poisson" %in% methods){
      all_res[[bus]][["od_poisson"]] = rmse_business_line(projection_fun = od_poisson, data=data, business_line = bus)
    }
    
    # Method 2:  GP Ladder 
    if ("gp_ladder" %in% methods){
      all_res[[bus]][["gp_ladder"]] = rmse_business_line(projection_fun = gp_ladder, data=data, business_line = bus, 
                             formula = ~ log(DevelopmentLag) + (AccidentYear), eps=0.1, plot_graphs= FALSE)
    }
    
    # Method 3:  Stan Ladder 
    if ("stan_ladder" %in% methods){
      all_res[[bus]][["stan_ladder"]] = rmse_business_line(projection_fun = stan_ladder, data=data, business_line = bus,
                                         formula_x = dev_log ~ AccidentYear_scale + DevelopmentYear_scale -1,
                                         formula_H = ~ log(DevelopmentYear_scale - min(DevelopmentYear_scale)+1) + AccidentYear_scale, 
                                         eps = 0.001, adapt_delta = 0.975, iter = 400, prior_theta=1)
    }
    # Method 4:  Stan Ladder Incremental Loss Ratio 
    if ("stan_loss_ratio_incr" %in% methods){
      all_res[[bus]][["stan_loss_ratio_incr"]]= rmse_business_line(projection_fun = stan_loss_ratio_incr, data=data, business_line = bus,
                                             formula_x = loss_ratio_scale ~ AccidentYear_scale + DevelopmentYear_scale -1,
                                             formula_H = ~ I(1/ (DevelopmentYear_scale - min(DevelopmentYear_scale)+1)) + AccidentYear_scale, 
                                             eps = 0.001, adapt_delta = 0.95, iter = 500, 
                                             prior_theta=1, prior_sigma = 0.5, 
                                             max_treedepth = 13, take_logs = FALSE,
                                             stanFile = "gpFastPredict04.stan")
    }
    
    # Method 5:  Stan Ladder using SE kernel, Zero Incremental Loss Ratio 
    if ("stan_ladder_se" %in% methods){
      
      all_res[[bus]][["stan_ladder_se"]]= rmse_business_line(projection_fun = stan_ladder_se,data=data, 
                                                             business_line = bus,eps = 0.0005, plot_graphs = FALSE,
                                                             chains = 6,iter = 400, warmup = 100, adapt_delta = 0.9,max_treedepth = 10,
                                                             prior_eta = 1, prior_rho = c(5.557,5.564), prior_sigma = 1, 
                                                             stanFile = "gp_se_03.stan")
    }
    
    if ("stan_compound" %in% methods){
      
      all_res[[bus]][["stan_compound"]]= rmse_business_line(projection_fun = stan_ladder_se,data=data, 
                                                            formula_H = ~ AccidentYear_scaled + DevelopmentLag_scale-1,
                                                            business_line = bus,eps = 0.0005, plot_graphs = FALSE,
                                                            chains = 6,iter = 400, warmup = 100, adapt_delta = 0.9,max_treedepth = 10,
                                                            prior_eta = 1, prior_rho = c(5.557,5.564), prior_theta = 1, prior_sigma = 1, 
                                                            stanFile = "gp_compound_03.stan")
    }
    
    if ("stan_compound_log" %in% methods){
      
      all_res[[bus]][["stan_compound_log"]]= rmse_business_line(projection_fun = stan_ladder_se,data=data, 
                                                            formula_H = ~ AccidentYear_scaled + DevelopmentLag_trans-1,
                                                            business_line = bus,eps = 0.0005, plot_graphs = FALSE,
                                                            chains = 6,iter = 400, warmup = 100, adapt_delta = 0.9,max_treedepth = 10,
                                                            prior_eta = 1, prior_rho = c(5.557,5.564), prior_theta = 1, prior_sigma = 1, 
                                                            stanFile = "gp_compound_03.stan")
    }
    if ("stan_comp_het" %in% methods){
      
      all_res[[bus]][["stan_comp_het"]]= rmse_business_line(projection_fun = stan_ladder_se,data=data, 
                                                            formula_H = ~ AccidentYear_scaled + DevelopmentLag_trans-1,
                                                            include_dev_lag = TRUE,
                                                            business_line = bus,eps = 0.0005, plot_graphs = FALSE,
                                                            chains = 6,iter = 400, warmup = 100, adapt_delta = 0.9,max_treedepth = 10,
                                                            prior_eta = 1, prior_rho = c(5.557,5.564), prior_theta = 1, prior_sigma = 1, prior_lambda = .05,
                                                            stanFile = "gp_comp_het.stan")
    }
    if ("stan_inc_het" %in% methods){
      
      all_res[[bus]][["stan_comp_het"]]= rmse_business_line(projection_fun = odp_gp,data=data, 
                                                            formula_H = ~ AccidentYear_scaled + DevelopmentLag_trans-1,
                                                            include_dev_lag = FALSE,
                                                            business_line = bus,eps = 0.0005, plot_graphs = FALSE,
                                                            chains = 1,iter = 100, warmup = 50, adapt_delta = 0.9,max_treedepth = 10,
                                                            prior_eta = 1, prior_rho = c(5.557,5.564), prior_theta = 1, prior_sigma = 1, prior_lambda = .05,
                                                            stanFile = "gp_comp_het.stan")
    }
    if ("stan_inc_vb" %in% methods){
      
      all_res[[bus]][["stan_inc_vb"]]= rmse_business_line(projection_fun = odp_gp_fac,data=data,use_vb = TRUE, 
                                                            include_dev_lag = FALSE,
                                                            business_line = bus,eps = 0.0005, plot_graphs = FALSE,
                                                            chains = 6,iter = 400, warmup = 200, adapt_delta = 0.9,max_treedepth = 10,
                                                            stanFile = "odp_gp_03.stan")
    }
    
    codes = all_res[[bus]][[1]]$code
    bus_res = map(all_res[[bus]], ~ .x[,-1])
    # bus_res = map2(bus_res, methods, function(x,y) {
    #   colnames(x)=paste0(y,"_", colnames(x))
    #   x
    # })
    bus_res = cbind(code = codes, line = bus, do.call(cbind, bus_res))
    

    if (is.null(final_res)){
      final_res = bus_res
    } else {
      final_res = rbind(final_res, bus_res)
    }
    save(final_res, file = paste0(output_dir, "interim_compare_all4.RData"))
    
    
  }
  final_res 
}





if (F){
  data = load_losses()
  bus = "prodliab"

  res_stan_ladder = rmse_business_line(projection_fun = stan_ladder, data=data, business_line = bus,
                                       formula_x = dev_log ~ AccidentYear_scale + DevelopmentYear_scale -1,
                                       formula_H = ~ log(DevelopmentYear) + AccidentYear, 
                                       eps = 0.001, adapt_delta = 0.99, iter = 400, prior_theta=1,
                                       max_treedepth = 13,
                                       stanFile = "gpFastPredict05.stan")
  
  #############################
  
  data = load_losses()
  bus = "wkcomp"
  dy_mean = mean(log(data$DevelopmentYear))
  dy_sd = sd(log(data$DevelopmentYear))
  
  res_stan_ladder = rmse_business_line(projection_fun = stan_ladder, data=data %>% filter(GRCODE == 86), business_line = bus,
                                       formula_x = dev_log ~ AccidentYear_scale + DevelopmentYear_scale -1,
                                       formula_H = ~ I((log(DevelopmentYear)-dy_mean)/dy_sd) + AccidentYear_scale, 
                                       eps = 0.001, adapt_delta = 0.99, iter = 400, prior_theta=1,
                                       max_treedepth = 13,
                                       stanFile = "gpFastPredict04.stan")
  
}

###   Run comparison across business lines
compare_cl_gp_stan = function(){
  
  # yyy

  a=Sys.time()
  ans = compare_methods(business_lines = c("comauto", "medmal", "othliab","ppauto","prodliab", "wkcomp"))  
  b=Sys.time(); b-a
  save(ans, file = paste0(output_dir, "compare_all3.RData"))
  #load(file = paste0(output_dir, "compare_all.RData"))
  #load(file = paste0(output_dir, "compare_all3.RData"))
  
  ans$code=as.factor(ans$code)
  ggplot(ans) + 
    geom_point(aes(code, cl_rmse_LR, color="Chain Ladder"))+ 
    #geom_point(aes(code, gp_rmse_LR, color="Dice Kriging"))+ 
    geom_point(aes(code, stan_rmse_LR, color="Stan MCMC")) + 
    labs(title= "Comparison of RMSE as a % of Premium", x="Company code", y="Error") + ylim(0, 1.5) + 
    facet_wrap(~line)
  ggplot(ans) + 
    geom_point(aes(code, cl_rmse_LR-gp_rmse_LR, color="Chain Ladder"))+ 
    #geom_point(aes(code, , color="Dice Kriging"))+ 
    #geom_point(aes(code, stan_rmse_LR, color="Stan MCMC")) + 
    labs(title= "Comparison of RMSE as a % of Premium", x="Company code", y="Error") + 
    #ylim(0, 1.5) + 
    facet_wrap(~line)
  ggplot(ans) + 
    geom_density(aes(cl_rmse_LR-stan_rmse_LR))+ 
    labs(title= "Density of Excess of Chain Ladder over Stan Error\nRMSE as a % of Premium", x="Excess of CL Error") + 
    xlim(-.25, .25) + 
    facet_wrap(~line)
  ggplot(ans) + 
    geom_density(aes(cl_rmse_LR-gp_rmse_LR))+ 
    labs(title= "Density of Excess of Chain Ladder over Dice Kriging Error\nRMSE as a % of Premium", x="Excess of CL Error") + 
    xlim(-.25, .25) + 
    facet_wrap(~line)
  
  ans %>%
    select(cl_rmse_LR, gp_rmse_LR, stan_rmse_LR, line) %>%
    group_by(line) %>%
    summarize(cl_med = median(cl_rmse_LR, na.rm = T), 
              gp_med = median(gp_rmse_LR, na.rm = T), 
              stan_med = median(stan_rmse_LR, na.rm = T))
  ans %>%
    select(cl_rmse_LR, gp_rmse_LR, stan_rmse_LR, line) %>%
    group_by(line) %>%
    summarize(cl_mean = mean(cl_rmse_LR, na.rm = T), 
              gp_mean = mean(gp_rmse_LR, na.rm = T), 
              stan_mean = mean(stan_rmse_LR, na.rm = T))
  ans %>%
    select(cl_rmse_LR, gp_rmse_LR, stan_rmse_LR, line) %>%
    group_by(line) %>%
    summarize(cl_sd = sd(cl_rmse_LR, na.rm = T), 
              gp_sd = sd(gp_rmse_LR, na.rm = T), 
              stan_sd = sd(stan_rmse_LR, na.rm = T))
}


##########  Greta Example
## Manually step through this function
greta_example = function(){
  
  ## Source file
  source('custom_gp.R')
  
  #  Load Data
  data = load_losses()
  
  # Calculate square and triangle tables
  loss  =  loss_table(data, "medmal", company_code = 669) #669, 683, 7854, 32514
  loss_triangle = loss$paid_tr
  
  # Create training arrays
  eps = 0.1
  AccidentYear = as.numeric(rownames(loss_triangle))
  DevelopmentYear = as.numeric(colnames(loss_triangle)[-1])
  no_acc = length(AccidentYear); no_dev = length(DevelopmentYear)
  AccidentYear = rep(AccidentYear, each = no_dev)
  DevelopmentYear = rep(DevelopmentYear, times = no_acc)
  
  accident_scaler = scaler(AccidentYear)
  development_scaler = scaler(DevelopmentYear)
  
  dev = dev_factors(loss_triangle) %>%
    as_tibble() %>%
    mutate(AccidentYear = as.integer(rownames(loss_triangle))) %>%
    gather(DevelopmentYear, dev, -AccidentYear) %>%
    filter(!is.na(dev)) %>%
    mutate(dev_log = log(dev-1+eps),
           DevelopmentYear = as.integer(DevelopmentYear),
           AccidentYear_scaled = accident_scaler(AccidentYear),
           DevelopmentYear_scaled = accident_scaler(DevelopmentYear))
  
  dev_scaler = scaler(dev$dev_log)
  dev$dev_log = c(dev_scaler(dev$dev_log))

  x2 = tibble(AccidentYear_scaled=c(accident_scaler(AccidentYear)) ,
              DevelopmentYear_scaled=c(development_scaler(DevelopmentYear)))
  
  model = gp_greta(x1 = dev[,c("AccidentYear_scaled", "DevelopmentYear_scaled")],n_samples = 200,
                   y =  dev[,c("dev_log")],
                   x2=x2)

  model = model %>% 
    left_join(dev) %>%
    mutate(pred_dev = dev_scaler(predicted, unscale = TRUE),
           pred_lower = dev_scaler(lower, unscale = TRUE),
           pred_upper = dev_scaler(upper, unscale = TRUE)
    )
  
  model
  
  ggplot(model) + 
    geom_point(aes(DevelopmentYear, dev, color = "actual")) +
    geom_line(aes(DevelopmentYear, pred_dev, color = "predicted")) +
    facet_wrap(~AccidentYear,scales = "free")
  
  
  
  loss_square = make_loss_square(loss_triangle = loss_triangle, dev_factors = model)
  
  rmse_triangle( actual_sq = loss$paid_sq, predicted_sq = loss_square, actual_tr = loss_triangle)
  
  plot_ae(actual_sq = loss$paid_sq, predicted_sq = loss_square)
  
  ### Compare to Chain Ladder
  # Chain Ladder Projection
  pred_sq = chain_ladder(loss$paid_tr)
  pred_sq
  rmse_triangle(actual_sq = loss$paid_sq, predicted_sq = pred_sq, actual_tr = loss$paid_tr)
  plot_ae(loss$paid_sq, pred_sq )
  
}



# Uses GP with Stan to project Loss Development Factors in order to estimate developed losses

stan_ladder = function(loss_triangle, 
                       formula_x = dev_log ~ AccidentYear_scale + DevelopmentLag_scale -1,
                       formula_H = ~ AccidentYear + log(DevelopmentLag), 
                       eps = 0.001, plot_graphs = FALSE,
                       chains = 6,iter = 600, warmup = iter/4, adapt_delta = 0.95,max_treedepth = 10,
                       prior_eta_sq = 1, prior_theta = 1, prior_sigma_sq = 1,
                       stanFile = "gpFastPredict04.stan", premiums = NULL){

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
  
  
  model =   gp_stan(data = dev, data_new = dev_new, optimizeMod = FALSE,chains = chains,iter = iter, warmup=warmup, 
                    formula_x = formula_x,
                    formula_H = formula_H,
                    prior_eta_sq = prior_eta_sq, prior_theta = prior_theta, prior_sigma_sq = prior_sigma_sq,
                    adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                    stanFile = stanFile) 
  
  unscale = function(x){
    exp(x* sd_dev_log + mean_dev_log) +1- eps  
  }
  
  preds = rstan::extract(model, pars= "yStar")[[1]]
  
  preds_summary = tibble(AccidentYear = AccidentYear, DevelopmentLag = DevelopmentLag, 
                             predicted = colMeans(preds), lower=quantile(preds, probs = 0.05), upper=quantile(preds, probs = 0.95))

  preds_summary = preds_summary %>% 
    left_join(dev) %>%
    mutate(pred_dev = unscale(predicted),
           pred_lower = unscale(lower),
           pred_upper = unscale(upper)
    )
  
  if (plot_graphs){
    print(ggplot(preds_summary) + 
      geom_point(aes(DevelopmentLag, dev, color = "actual")) +
      geom_line(aes(DevelopmentLag, pred_dev, color = "predicted")) +
      facet_wrap(~AccidentYear,scales = "free"))
    
    #plot_ae(actual_sq = loss$paid_sq, predicted_sq = loss_square)
    
  }

  loss_square = make_loss_square(loss_triangle = loss_triangle, dev_factors = preds_summary)
  
  loss_square
  
}


# Uses GP with Stan to project Loss Development Factors in order to estimate developed losses
# Utlizes Square Exponential Kernel, Zero Mean Prior

stan_ladder_se = function(loss_triangle, 
                       #formula_x = dev_log ~ AccidentYear_scale + DevelopmentLag_scale -1,
                       formula_H = NULL, include_dev_lag = FALSE, 
                       eps = 0.0005, plot_graphs = FALSE,
                       chains = 6,iter = 400, warmup = 100, adapt_delta = 0.9,max_treedepth = 10,
                       prior_eta = 1, prior_rho = c(5.557,5.564), prior_theta = 1, prior_sigma = 1, prior_lambda=1,
                       optimizeMod = FALSE,
                       stanFile = "gp_se_03.stan", premiums = NULL){
  
  # calculate LDFs
  ldf = dev_factors(loss_triangle)
  
  # Convert Tables to long format
  loss_long = losses_long(list(loss_triangle))[[1]]
  ldf_long = losses_long(list(ldf), val_name = "ldf")[[1]]
  
  #Join LDF and Loss Tables
  loss_long = left_join(loss_long, ldf_long) %>%
    remove_missing()
  
  # Scale Data
  ## scaler is a utlity function for scaling and unscaling data
  scaler_accident_year = scaler(loss_long$AccidentYear)
  scaler_development_lag = scaler(loss_long$DevelopmentLag)
  scaler_ldf = scaler(loss_long$ldf)
  DevelopmentLag_transform = function(x) log(x-min(loss_long$DevelopmentLag)+1) 
  scaler_df = function(df, eps = 0.0005){
    df$AccidentYear_scaled = scaler_accident_year(df$AccidentYear)[,1]
    df$DevelopmentLag_scaled = scaler_development_lag(df$DevelopmentLag)[,1]
    df$DevelopmentLag_trans = DevelopmentLag_transform(df$DevelopmentLag)
    df$ldf_scaled = scaler_ldf(df$ldf)[,1]
    df$ldf_trans = log(pmax(df$ldf,1)-1+eps)
    df
  }
  
  loss_long = scaler_df(loss_long, eps=eps)
  
  
  #############################
  AccidentYear = as.numeric(rownames(loss_triangle))
  DevelopmentLag = as.numeric(colnames(loss_triangle))[-1]
  no_acc = length(AccidentYear); no_dev = length(DevelopmentLag)
  AccidentYear = rep(AccidentYear, each = no_dev)
  DevelopmentLag = rep(DevelopmentLag, times = no_acc)
  
  data_new = data.frame(AccidentYear_scaled = scaler_accident_year(AccidentYear)[,1],
                        DevelopmentLag_scaled = scaler_development_lag(DevelopmentLag)[,1],
                        DevelopmentLag=DevelopmentLag) %>%
    mutate(DevelopmentLag_trans = DevelopmentLag_transform(DevelopmentLag))
  
  
  a= Sys.time()
  model = gp_stan_se (data=loss_long, data_new = data_new, 
                      formula_x = ldf_trans ~  AccidentYear_scaled + DevelopmentLag_scaled -1,
                      formula_H =  formula_H, 
                      include_dev_lag = include_dev_lag,
                      #formula_cat =  ~ GRCODE, 
                      prior_eta = prior_eta, prior_rho = prior_rho, prior_sigma = prior_sigma, 
                      prior_lambda = prior_lambda, prior_theta = prior_theta,
                      adapt_delta=adapt_delta, max_treedepth = max_treedepth,
                      optimizeMod = optimizeMod, iter =iter, chains = chains, warmup =warmup,
                      stanFile = stanFile)
  b= Sys.time(); b-a
  
  
  preds = rstan::extract(model$draws, pars= "yStar")[[1]]
  
  preds_summary = tibble(AccidentYear = AccidentYear, DevelopmentLag = DevelopmentLag, 
                         predicted = colMeans(preds)) %>%
    mutate(pred_dev = exp(predicted) +1-eps)
  
  loss_square = make_loss_square(loss_triangle = loss_triangle, dev_factors = preds_summary)
  
  loss_square
  
}




# Uses GP with Stan to project cumulative loss in order to estimate developed losses
# Models the incremental loss ratio

stan_loss_ratio = function(loss_triangle, premiums=NULL,
                       formula_x = loss_ratio_scale ~ AccidentYear_scale + DevelopmentYear_scale -1,
                       formula_H = ~ AccidentYear + log(DevelopmentYear), 
                       eps = 0.001, plot_graphs = FALSE,
                       chains = 6,iter = 600, adapt_delta = 0.95,max_treedepth = 10,
                       prior_eta_sq = 1, prior_theta = 1, prior_sigma_sq = 1,
                       stanFile = "gpFastPredict04.stan"){
  
  AccidentYear = as.numeric(rownames(loss_triangle))
  DevelopmentYear = as.numeric(colnames(loss_triangle))[-1]
  no_acc = length(AccidentYear); no_dev = length(DevelopmentYear)
  AccidentYear = rep(AccidentYear, each = no_dev)
  DevelopmentYear = rep(DevelopmentYear, times = no_acc)
  
  loss_ratio = as_tibble(loss_triangle / premiums) %>%
    mutate(AccidentYear = as.integer(rownames(loss_triangle))) %>%
    gather(DevelopmentYear, loss_ratio, -AccidentYear) %>%
    filter(!is.na(loss_ratio)) %>%
    mutate(loss_ratio_log = log(pmax(eps,loss_ratio) ),
           DevelopmentYear = as.integer(DevelopmentYear))
  
  AccidentYear_mean = mean(loss_ratio$AccidentYear)
  AccidentYear_sd = sd(loss_ratio$AccidentYear)
  DevelopmentYear_mean = mean(loss_ratio$DevelopmentYear)
  DevelopmentYear_sd = sd(loss_ratio$DevelopmentYear)
  
  loss_ratio = loss_ratio %>%
    mutate(AccidentYear_scale = (AccidentYear - AccidentYear_mean)/AccidentYear_sd,
           DevelopmentYear_scale = (DevelopmentYear - DevelopmentYear_mean)/DevelopmentYear_sd)
  
  mean_loss_ratio = mean(loss_ratio$loss_ratio)
  sd_loss_ratio = sd(loss_ratio$loss_ratio)
  loss_ratio$loss_ratio_scale = (loss_ratio$loss_ratio - mean_loss_ratio)/sd_loss_ratio
  
  # Create inputs for out-of-sample predictions
  loss_ratio_new = tibble(AccidentYear ,DevelopmentYear) %>%
    mutate(AccidentYear_scale = (AccidentYear - AccidentYear_mean)/AccidentYear_sd,
           DevelopmentYear_scale = (DevelopmentYear - DevelopmentYear_mean)/DevelopmentYear_sd)
  
  
  model =   gp_stan(data = loss_ratio, data_new = loss_ratio_new, optimizeMod = FALSE,chains = chains,iter = iter,
                    formula_x = formula_x,
                    formula_H = formula_H,
                    prior_eta_sq = prior_eta_sq, prior_theta = prior_theta, prior_sigma_sq = prior_sigma_sq,
                    adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                    stanFile = stanFile) 
  pairs(model, pars = c("eta_sq", "theta"))
  unscale = function(x){
    x* sd_loss_ratio + mean_loss_ratio
  }
  
  preds = unscale(rstan::extract(model, pars= "yStar")[[1]])
  
  preds_summary = tibble(AccidentYear = AccidentYear, DevelopmentYear = DevelopmentYear, 
                         premiums = rep(premiums, each = no_dev),
                         predicted = colMeans(preds), 
                         lower=apply(preds,2,quantile,probs = 0.05), upper=apply(preds,2,quantile,probs = 0.95))
  
  preds_summary = preds_summary %>% 
    left_join(loss_ratio) %>%
    mutate(pred_loss = (predicted) * premiums,
           pred_lower = (lower)* premiums,
           pred_upper = (upper)* premiums,
           actual = loss_ratio *premiums
    )
  
  if (plot_graphs){
    print(ggplot(preds_summary) + 
            geom_point(aes(DevelopmentYear, actual , color = "actual")) +
            geom_line(aes(DevelopmentYear, pred_loss, color = "predicted")) +
            facet_wrap(~AccidentYear,scales = "free"))
    print(ggplot(preds_summary) + 
            geom_point(aes(DevelopmentYear, actual , color = as.factor(AccidentYear))) +
            geom_line(aes(DevelopmentYear, pred_loss, color = as.factor(AccidentYear))) + 
            labs(title = "GP Fit to Cumulative Loss Ratio/nCumulative Losses Shown")
    )
    print(ggplot(preds_summary) + 
            geom_point(aes(DevelopmentYear, loss_ratio , color = as.factor(AccidentYear))) +
            geom_line(aes(DevelopmentYear, predicted, color = as.factor(AccidentYear))) + 
            labs(title = "GP Fit to Cumulative Loss Ratio/Cumlative Loss Ratios Shown")
    )
    
    
    #plot_ae(actual_sq = loss$paid_sq, predicted_sq = loss_square)
    
  }
  
  loss_square = preds_summary %>%
    select(AccidentYear, DevelopmentYear, pred_loss) %>%
    spread(DevelopmentYear, pred_loss)
  
  loss_square
  
  }


# Uses GP with Stan to project cumulative loss in order to estimate developed losses

stan_loss_ratio_incr = function(loss_triangle, premiums=NULL,
                           formula_x = loss_ratio_scale ~ AccidentYear_scale + DevelopmentYear_scale -1,
                           formula_H = ~ AccidentYear + log(DevelopmentYear), 
                           eps = 0.001, plot_graphs = FALSE,
                           chains = 6,iter = 600, warmup=200, adapt_delta = 0.95,max_treedepth = 10,
                           prior_eta_sq = 1, prior_theta = 1, prior_sigma_sq = 1, take_logs = FALSE, 
                           stanFile = "gpFastPredict04.stan"){

  AccidentYear_unique = AccidentYear = as.numeric(rownames(loss_triangle))
  DevelopmentYear_unique = DevelopmentYear = as.numeric(colnames(loss_triangle))
  DevelopmentYear = DevelopmentYear[-1]
  no_acc = length(AccidentYear); 
  no_dev = length(DevelopmentYear)
  AccidentYear = rep(AccidentYear, each = no_dev)
  DevelopmentYear = rep(DevelopmentYear, times = no_acc)
  
  loss_ratio_triangle = (loss_triangle / premiums)
  loss_ratio_triangle = loss_ratio_triangle[,-1] - loss_ratio_triangle[, -ncol(loss_ratio_triangle)]
  loss_ratio_triangle_log = log(apply(loss_triangle / premiums,2,function(x) pmax(x,eps)))
  loss_ratio_triangle_log = loss_ratio_triangle_log[,-1] - loss_ratio_triangle_log[, -ncol(loss_ratio_triangle_log)]
  
  
  loss_ratio = as_tibble(loss_ratio_triangle) %>%
    mutate(AccidentYear = as.integer(rownames(loss_ratio_triangle))) %>%
    gather(DevelopmentYear, loss_ratio, -AccidentYear) %>%
    filter(!is.na(loss_ratio)) %>%
    mutate(
      DevelopmentYear = as.integer(DevelopmentYear)
    )
  loss_ratio_log = as_tibble(loss_ratio_triangle_log) %>%
    mutate(AccidentYear = as.integer(rownames(loss_ratio_triangle_log))) %>%
    gather(DevelopmentYear, loss_ratio_log, -AccidentYear) %>%
    filter(!is.na(loss_ratio_log)) %>%
    mutate(
      DevelopmentYear = as.integer(DevelopmentYear)
    )
  loss_ratio = loss_ratio %>%
    left_join(loss_ratio_log)
  
  #ggplot(loss_ratio_log[loss_ratio_log$DevelopmentYear>1,]) + geom_line(aes(DevelopmentYear, loss_ratio_log, color=as.factor(AccidentYear)))
  #ggplot(loss_ratio) + geom_line(aes(as.integer(DevelopmentYear), loss_ratio, color=as.factor(AccidentYear)))
  
  AccidentYear_mean = mean(loss_ratio$AccidentYear)
  AccidentYear_sd = sd(loss_ratio$AccidentYear)
  DevelopmentYear_mean = mean(loss_ratio$DevelopmentYear)
  DevelopmentYear_sd = sd(loss_ratio$DevelopmentYear)
  
  loss_ratio = loss_ratio %>%
    mutate(AccidentYear_scale = (AccidentYear - AccidentYear_mean)/AccidentYear_sd,
           DevelopmentYear_scale = (DevelopmentYear - DevelopmentYear_mean)/DevelopmentYear_sd)
  
  mean_loss_ratio = mean(loss_ratio$loss_ratio)
  sd_loss_ratio = sd(loss_ratio$loss_ratio)
  loss_ratio$loss_ratio_scale = (loss_ratio$loss_ratio - mean_loss_ratio)/sd_loss_ratio
  mean_loss_ratio_log = mean(loss_ratio$loss_ratio_log)
  sd_loss_ratio_log = sd(loss_ratio$loss_ratio_log)
  loss_ratio$loss_ratio_log_scale = (loss_ratio$loss_ratio_log - mean_loss_ratio_log)/sd_loss_ratio_log
  
  # Create inputs for out-of-sample predictions
  loss_ratio_new = tibble(AccidentYear ,DevelopmentYear) %>%
    mutate(AccidentYear_scale = (AccidentYear - AccidentYear_mean)/AccidentYear_sd,
           DevelopmentYear_scale = (DevelopmentYear - DevelopmentYear_mean)/DevelopmentYear_sd)
  loss_ratio_dev_year_min = min(loss_ratio_new$DevelopmentYear_scale)
  loss_ratio_new$DevelopmentYear_scale_log = loss_ratio_new$DevelopmentYear_scale - loss_ratio_dev_year_min + 1
  loss_ratio$DevelopmentYear_scale_log = loss_ratio$DevelopmentYear_scale - loss_ratio_dev_year_min + 1
  
  
  model =   gp_stan(data = loss_ratio, data_new = loss_ratio_new, optimizeMod = FALSE,chains = chains,iter = iter,
                    formula_x = formula_x,
                    formula_H = formula_H,
                    prior_eta_sq = prior_eta_sq, prior_theta = prior_theta, prior_sigma_sq = prior_sigma_sq,
                    adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                    stanFile = stanFile) 
  #pairs(model, pars = c("eta_sq", "theta"))
  unscale = function(x){
    if (take_logs){
      x=x* sd_loss_ratio_log + mean_loss_ratio_log
    } else {
      x=x* sd_loss_ratio + mean_loss_ratio
    }
    x
  }
  
  preds = unscale(rstan::extract(model, pars= "yStar")[[1]])
  
  if (take_logs){
    loss_ratio_first = log(pmax(eps, loss_triangle[,1]/ premiums))
  } else {
    loss_ratio_first = loss_triangle[,1]/ premiums
  }  
  
  if (F){ #check
    aa= matrix(colMeans(preds), ncol = ncol(loss_ratio_triangle_log), byrow = T)
    aa
    aa/ loss_ratio_triangle_log
    aa= cbind(loss_triangle[,1]/ premiums, exp(aa))
    t(apply(aa,1,cumprod))
  }
  
  no_iter = nrow(preds)

  preds_cum = matrix(NA,nrow=no_iter, ncol = no_acc * (no_dev+1))
  
  for (acc_year in (1:no_acc-1)){
    for (dev_year in 1:no_dev){
      j = acc_year * no_dev + dev_year
      j_star = acc_year * (no_dev+1) + dev_year
      if (dev_year == 1){
        preds_cum[, j_star] = loss_ratio_first[acc_year+1]
      }
      preds_cum[,j_star+1] = preds[,j] + preds_cum[,j_star]
    }
  }
  
  if (take_logs){
    preds_cum = exp(preds_cum)
  }
  
  AccidentYear = rep(AccidentYear, each = no_dev)
  DevelopmentYear = rep(DevelopmentYear, times = no_acc)
  
  ## Following uses median prediction, not mean.
  preds_summary = tibble(AccidentYear = rep(AccidentYear_unique, each = no_dev+1), 
                         DevelopmentYear = rep(DevelopmentYear_unique, times = no_acc), 
                         premiums = rep(premiums, each = no_dev+1),
                         predicted_lr = apply(preds_cum,2,quantile,probs = 0.5), #colMeans(preds_cum), 
                         lower_lr=apply(preds_cum,2,quantile,probs = 0.05), upper_lr=apply(preds_cum,2,quantile,probs = 0.95))
  


  loss_triangle_long = as_tibble(loss_triangle/premiums) %>%
    mutate(AccidentYear = AccidentYear_unique) %>%
    gather(DevelopmentYear, actual_lr, -AccidentYear)  %>%
    mutate(DevelopmentYear = as.numeric(DevelopmentYear))

  preds_summary = preds_summary %>%
    left_join(loss_triangle_long) %>%
    mutate(predicted = predicted_lr * premiums,
           lower = lower_lr * premiums,
           upper = upper_lr * premiums,
           actual = actual_lr *premiums)
    
  if (plot_graphs){
    print(ggplot(preds_summary) + 
            geom_point(aes(DevelopmentYear, actual , color = "actual")) +
            geom_line(aes(DevelopmentYear, predicted, color = "predicted")) +
            facet_wrap(~AccidentYear,scales = "free"))
    print(ggplot(preds_summary) + 
            geom_point(aes(DevelopmentYear, actual , color = as.factor(AccidentYear))) +
            geom_line(aes(DevelopmentYear, predicted, color = as.factor(AccidentYear))) + 
            labs(title = "GP Fit to Cumulative Loss Ratio/nCumulative Losses Shown")
    )
    print(ggplot(preds_summary) + 
            geom_point(aes(DevelopmentYear, actual_lr , color = as.factor(AccidentYear))) +
            geom_line(aes(DevelopmentYear, predicted_lr, color = as.factor(AccidentYear))) + 
            labs(title = "GP Fit to Cumulative Loss Ratio\nCumlative Loss Ratios Shown")
    )
    
    
    #plot_ae(actual_sq = loss$paid_sq, predicted_sq = loss_square)
    
  }
  
  loss_square = preds_summary %>%
    select(AccidentYear, DevelopmentYear, predicted) %>%
    spread(DevelopmentYear, predicted) %>%
    select(-AccidentYear) %>%
    as.matrix()
  
  loss_square
  
}



#' Fit of  a Gaussian Process with a Trend Function and a single categorical variable.
#' 
#' Prediction can be done on out-of-sample datasets
#' 
#' @param data A data frame containing the input training set.  
#' @param data_new The out-of-sample prediction set.  If NULL, data_new is set to data.  No target values are required.
#' @param formula_x The GP regression formula
#' @param forumla_H The right hand side formula for the predictors of the mean function
#' @param forumla_cat The right hand side formula for the single categorical predictor.  The designated column must contain real numbers.
#' @param prior_eta_sq Prior for eta_sq parameter
#' @param prior_theta Prior for theta parameter
#' @param prior_sigma_sq Prior for sigma_sq parameter
#' @param optimizeMod Default is \code{TRUE}.  If \TRUE}, then an optimized set of parameters are provided.  If \code{FALSE}, then MCMC is conducted
#' @param iter Number of iterations if MCMC is run
#' @param chains Number of chains if MCMC is run
#' @return The MCMC chains or the optimized valus of the parameters of the Gaussian Process model
gp_stan_cat2 <- function(data, data_new = NULL, formula_x, formula_H = NULL, formula_cat = NULL, 
                         prior_eta = 1, prior_rho = c(5.5,5.5), prior_theta = 1, prior_alpha = 1, prior_sigma = 1, 
                         adapt_delta=0.9, max_treedepth = 10,verbose=TRUE,
                         optimizeMod = FALSE, iter =100, chains = 1, warmup = 75,
                         stanFile = "gp_comp_cat.stan"){
  
  if (is.null(data_new)){
    data_new = data
  }
  
  data=as_tibble(data)
  data_new= as_tibble(data_new)
  
  # Establish data to be passed to Stan
  N1 = nrow(data)
  N2 = nrow(data_new)
  
  data_new$blank = 0
  #data$blank = 0
  
  formula_x = update(formula_x, . ~ . -1) # remove intercept term
  
  stan_list <- list(
    N1 = N1,
    N2 = N2,
    x1 = model.matrix(formula_x, data),
    H1 = t(model.matrix(formula_H, data)), 
    y1 = c(data[,all.vars(update(formula_x, .~0)) ][[1]]),
    x2 = model.matrix(update(formula_x, blank ~ .), data_new),
    H2 = t(model.matrix(formula_H, data_new)), 
    prior_eta=prior_eta,
    prior_theta=prior_theta,
    prior_rho=prior_rho,
    prior_alpha=prior_alpha,
    prior_sigma=prior_sigma
  )
  
  if(!is.null(formula_cat)){
    formula_cat = update(formula_cat,  ~ . -1) # remove intercept term
    stan_list$C1 = (model.matrix(formula_cat, data))
    stan_list$C2 = (model.matrix(formula_cat, data_new))
    stan_list$dimC <- ncol(stan_list$C1)
  }
  
  stan_list$D = ncol(stan_list$x1)
  stan_list$dimH <- nrow(stan_list$H1)

  stanMod <- stan_model(file=paste0(stan_dir, stanFile))
  
  if (!optimizeMod){
    stanRun <- sampling(object = stanMod, data = stan_list, iter=iter, warmup=warmup, chains=chains, 
                        verbose=verbose, control = list(adapt_delta=adapt_delta, max_treedepth =max_treedepth) )
  } else {
    stanRun <- optimizing(object = stanMod, data = stan_list, verbose=verbose)
  }
  #stanRun
  #plot(stanRun, pars = "yStar")
  if(F){
    summary(stanRun, pars = c("eta", "rho", "theta", "alpha"))$summary
    
    preds = colMeans(rstan::extract(stanRun, pars = "yStar")[[1]])
    str(preds)
    data_new$preds =  preds
    head(data_new)
    
    ##### Loss Reserve Check
    ggplot(data_new) + 
      geom_line(aes(DevelopmentYear, preds, color = as.factor(AccidentYear))) + 
      facet_wrap(~ as.factor(GRCODE))
    ggplot(data_new) + 
      geom_line(aes(DevelopmentLag, preds, color = as.factor(AccidentYear))) + 
      facet_wrap(~ as.factor(GRCODE))
    
  }
  
  ans = list(draws = stanRun, mu_pred = rstan::extract(stanRun, pars = "yStar")[[1]], tau_pred = NA)
  ans
}


#' Fit of  a Gaussian Process with a Trend Function and a single categorical variable.
#' 
#' Prediction can be done on out-of-sample datasets
#' 
#' @param data A data frame containing the input training set.  
#' @param data_new The out-of-sample prediction set.  If NULL, data_new is set to data.  No target values are required.
#' @param formula_x The GP regression formula
#' @param forumla_H The right hand side formula for the predictors of the mean function
#' @param forumla_cat The right hand side formula for the single categorical predictor.  The designated column must contain real numbers.
#' @param prior_eta_sq Prior for eta_sq parameter
#' @param prior_theta Prior for theta parameter
#' @param prior_sigma_sq Prior for sigma_sq parameter
#' @param optimizeMod Default is \code{TRUE}.  If \TRUE}, then an optimized set of parameters are provided.  If \code{FALSE}, then MCMC is conducted
#' @param iter Number of iterations if MCMC is run
#' @param chains Number of chains if MCMC is run
#' @return The MCMC chains or the optimized valus of the parameters of the Gaussian Process model
gp_stan_cat_kron <- function(data, data_new = NULL, 
                             first_dim_vars = c("AccidentYear_scaled", "DevelopmentLag_scaled"),
                             basis_vars = c("AccidentYear_scaled", "DevelopmentLag_scaled"),
                             cat_var = "GRCODE", 
                             dev_lag_var = NULL,
                             response_var = "ldf_trans",
                             prior_eta = 1, prior_rho = c(5.5,5.5), prior_theta = 1, prior_alpha = 1, 
                             prior_sigma = 1,prior_lambda = 1, 
                             adapt_delta=0.9, max_treedepth = 10,
                             optimizeMod = FALSE, iter =100, chains = 1, warmup = 75,
                             stanFile = "gp_comp_cat_kron.stan", 
                             return_fStar = FALSE){
  
  if (is.null(data_new)){
    data_new = data
  }
  data=as_tibble(data)
  data_new_long = data_new= as_tibble(data_new)
  
  
  # Establish data to be passed to Stan
  union_vars = union(union(first_dim_vars, basis_vars), dev_lag_var)
  
  #data = data[, c(first_dim_vars, cat_var,response_var)] 
  data = data %>%
    select_(.dots = c(union_vars, cat_var,response_var))  %>%
    spread_(cat_var,response_var)
  
  data_new = data_new %>%
    select_(.dots = c(union_vars, cat_var,response_var))  %>%
    spread_(cat_var,response_var)
  
  head(data)
  head(data_new)
  
  cat_classes_train = factor(colnames(data %>% select(- one_of(union_vars))))
  cat_classes_test = factor(colnames(data_new %>% select(- one_of(union_vars))), levels = levels(cat_classes_train))
  
  if (!all (cat_classes_test %in% cat_classes_train)){
    stop("Test set has categories that are not in the training set")
  }
  
  
  # prepare the data
  stan_list <- list(data, 
                    N_train_1 = nrow(data),
                    N_train_2 = ncol(data) - length(union_vars),
                    N_test_1 = nrow(data_new), 
                    N_test_2 = ncol(data_new) - length(union_vars),
                    
                    x_train_1 = data[, first_dim_vars],
                    H_train_1 = t(data[, basis_vars]), 
                    cat_train_2 = one_hot(cat_classes_train),
                    
                    y = t(as.matrix(select(data, -one_of(union_vars)))),
                    
                    x_test_1 = data_new[, first_dim_vars],
                    H_test_1 = t(data_new[, basis_vars]),
                    cat_test_2 = one_hot(cat_classes_test),
                    
                    prior_eta=prior_eta,
                    prior_theta=prior_theta,
                    prior_rho=prior_rho,
                    prior_alpha=prior_alpha,
                    prior_sigma=prior_sigma,
                    prior_lambda=prior_lambda
  ) 
  
  if (!is.null(dev_lag_var)){
    stan_list$dev_lag_train_1 = unlist(data[, dev_lag_var])
    stan_list$dev_lag_test_1 = unlist(data_new[, dev_lag_var])
    stan_list$dev_lag_length = max(stan_list$dev_lag_train_1)
    
  }
  stan_list$D = ncol(stan_list$x_train_1)
  stan_list$dimH <- nrow(stan_list$H_train_1)
  stan_list$dimC <- ncol(stan_list$cat_train_2)
  
  stanMod <- stan_model(file=paste0(stan_dir, stanFile))
  
  if (!optimizeMod){
    stanRun <- sampling(object = stanMod, data = stan_list, iter=iter, warmup=warmup, chains=chains, control = list(adapt_delta=adapt_delta, max_treedepth =max_treedepth) )
  } else {
    stanRun <- optimizing(object = stanMod, data = stan_list, verbose=T)
  }
  summary(stanRun)$summary
  plot(stanRun, pars = "yStar")
  
  if (return_fStar){
    preds = t(rstan::extract(stanRun, pars = "fStar")[[1]])
  } else {
    preds = t(rstan::extract(stanRun, pars = "yStar")[[1]])
  }
  str(preds)
  
  reorder_output = function(vec, dim_1, dim_2){
    vec = c(t(matrix(vec, dim_2, dim_1)))
    vec
  }
  
  preds = apply(preds, 2, reorder_output, dim_1 = stan_list$N_train_1, dim_2 = stan_list$N_train_2)
  colnames(preds) = paste0("sim_",1:ncol(preds))
  data_long = (data_new %>%
                 gather(GRCODE, ldf_trans, -one_of(union_vars)))
  
  data_long = cbind(data_long, pred_trans = rowMeans(preds), preds)
  
  str(data_long)
  
  if(F){
    stan_trace(stanRun)
    summary(stanRun, pars = c("eta", "rho", "theta", "alpha"))$summary
    
    ##### Loss Reserve Check
    
    unique(data_new_long$GRCODE)[1:9]
    ggplot(data_new_long %>% filter(GRCODE %in% unique(data_new_long$GRCODE)[1:9])) + 
      geom_line(aes(DevelopmentLag_scaled, ldf_trans, color = as.factor(AccidentYear_scaled))) + 
      facet_wrap(~ as.factor(GRCODE), scales = "free")
    
    
  }
  
  ans = list(draws = stanRun, pred = data_long, tau_pred = NA)
  ans
}



#' Overdispersed Poisson Method for Projecting the Loss Triangle
#' using additive SE kernels for AY and DL factor variables.
#' 
#' @param loss_triangle
#' @param premiums A Vector of premiums. 
odp_gp_fac = function(loss_triangle, premiums = NULL, 
                      formula_H =  AccidentYear_scaled + DevelopmentLag_scaled, 
                      include_dev_lag = FALSE, include_premiums = TRUE,
                      prior_eta = 1, prior_rho = c(5.5,5.5), prior_sigma = 1, 
                      prior_lambda = .2, prior_theta = 1, prior_mean = 0, 
                      adapt_delta=.9, max_treedepth = 10,
                      optimizeMod = FALSE, iter =100, chains = 6, warmup =50,
                      stanFile = "odp_gp_03.stan",use_vb = FALSE,
                      ...){
  
  
  data = incremental_claims(loss_triangle) %>% 
    as.data.frame() %>%
    map_df( ~ pmax(.x,0))
  data = cbind(data, premiums = premiums)
  
  data = data %>%
    mutate(AccidentYear = as.numeric(rownames(data))) %>%
    gather(DevelopmentLag, y, - AccidentYear, -premiums) %>%
    mutate(DevelopmentLag = as.numeric(DevelopmentLag)) 
  
  data_new = data[is.na(data$y),]
  
  
  data = data %>%
    remove_missing() %>%
    mutate(y = pmax(y, 0))
  
  scaler_accident_year = scaler(data$AccidentYear)
  scaler_development_lag = scaler(data$DevelopmentLag)
  DevelopmentLag_transform = function(x) log(x-min(data$DevelopmentLag)+1) 
  
  
  data = data %>%
    mutate(AccidentYear_scaled = scaler_accident_year(AccidentYear)[,1],
           DevelopmentLag_scaled = scaler_development_lag(DevelopmentLag)[,1], 
           DevelopmentLag_trans = DevelopmentLag_transform(DevelopmentLag))
  data_new = data_new %>%
    mutate(AccidentYear_scaled = scaler_accident_year(AccidentYear)[,1],
           DevelopmentLag_scaled = scaler_development_lag(DevelopmentLag)[,1], 
           DevelopmentLag_trans = DevelopmentLag_transform(DevelopmentLag)) 
  #################################
  
  data=as_tibble(data)
  data_new= as_tibble(data_new)
  
  # Establish data to be passed to Stan
  
  x_ay1 = one_hot(factor(data$AccidentYear))
  x_dl1 = one_hot(factor(data$DevelopmentLag))
  x_ay2 = one_hot(factor(data_new$AccidentYear))
  x_dl2 = one_hot(factor(data_new$DevelopmentLag))
  
  formula_x = y ~ AccidentYear_scaled + DevelopmentYear_lag -1
  formula_H = ~ AccidentYear_scaled + DevelopmentLag_scaled
  formula_cat = NULL
  
  prior_eta=1
  prior_rho_ay= c(5.5, 5.5)
  prior_rho_dl= c(5.5, 5.5)
  
  data_new$blank = 0
  formula_x = update(formula_x, . ~ . -1) # remove intercept term
  
  stan_list <- list(
    N1 = nrow(data),
    N2 = nrow(data_new),
    D_ay = ncol(x_ay1),
    D_dl = ncol(x_dl1),
    #x1 = model.matrix(formula_x, data),
    x_ay1 = x_ay1,
    x_dl1 = x_dl1,
    y1 = c(data[,all.vars(update(formula_x, .~0)) ][[1]]),
    #x2 = model.matrix(update(formula_x, blank ~ .), data_new),
    x_ay2 = x_ay2,
    x_dl2 = x_dl2,
    prior_eta=prior_eta,
    prior_rho_ay=prior_rho_ay,
    prior_rho_dl = prior_rho_dl,
    prior_theta = prior_theta,
    prior_mean = prior_mean
    #prior_sigma=prior_sigma,
    #prior_lambda = prior_lambda
  )
  
  if(!is.null(formula_H)){
    stan_list$H1 = t(model.matrix(formula_H, data))
    stan_list$H2 = t(model.matrix(formula_H, data_new))
    stan_list$dimH <- nrow(stan_list$H1)
  }
  
  if (include_dev_lag){
    stan_list$dev_lag1 = data$DevelopmentLag
    stan_list$dev_lag2 = data_new$DevelopmentLag
  }
  if (include_premiums){
    stan_list$premiums1 = data$premiums
    stan_list$premiums2 = data_new$premiums
  }
  
  if(!is.null(formula_cat)){
    formula_cat = update(formula_cat,  ~ . -1) # remove intercept term
    stan_list$x_cat1 = c(model.matrix(formula_cat, data))
    stan_list$x_cat2 = c(model.matrix(formula_cat, data_new))
  }
  
  
  stanMod <- stan_model(file=paste0(stan_dir, stanFile))
  
  
  if (!use_vb){
    stanRun <- sampling(object = stanMod, data = stan_list, iter=iter, warmup=warmup, chains=chains,
                        control = list(adapt_delta=adapt_delta, max_treedepth =max_treedepth) )
  } else{
    stanRun <- vb(object = stanMod, data = stan_list )
  }
  
  
  if(F){
    summary(stanRun, pars = c("eta_ay","eta_dl","rho_ay","rho_dl"))$summary
    summary(stanRunVB, pars = c("eta_ay","eta_dl","rho_ay","rho_dl"))$summary
    stan_trace(stanRun)
    data_new$y = colMeans(rstan::extract(stanRunVB, pars= "yStar")[[1]])
    
  }
  
  data_new$y = colMeans(rstan::extract(stanRun, pars= "yStar")[[1]])
  data_new$blank = NULL
  preds = rbind(data, data_new) %>%
    #mutate(DevelopmentLag = as.numeric(as.character(DevelopmentLag))) %>%
    select(AccidentYear, DevelopmentLag, y) %>%
    spread(DevelopmentLag, y) %>%
    select( - AccidentYear) %>%
    as.matrix()
  
  losses_projected = project_losses_incr_claims(loss_triangle, incr_claims = preds)
  losses_projected
  
}




#' Calculate the predicted loss ratio and supporting analytics using an incremental loss ratio model
#' 
#' @param loss_triangle the Loss triangle and loss square
#' @param premiums Vector of premiums
#' @param calc_ranks If TRUE then additional statistics are calculated
#' @param formula_x formula for square exponential 
#' @param formula_H formula for linear terms in kernel (not in mean function)
#' @param take_logs If TRUE, then response is logarithm is used as the predictor.  Note, if TRUE, then LHS of formula_x must be y_log_scaled, otherwise use y_scaled
#' @param plot_graphs If TRUE, then sample investigative graphs are plotted.
#' @return If calc_ranks = FALSE, then only the predicted loss_square is returned.  If TRUE, then the percentage ranks, sigma estimates are also returned.
stan_incremental_loss_ratio = function(loss_triangle, premiums=NULL, loss_square = NULL, calc_ranks = FALSE, 
                                       formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                       formula_H = y_scaled ~  AccidentYear_scaled + DevelopmentLag_scaled -1 , formula_cat = NULL,
                                       eps = 0.001, plot_graphs = FALSE,
                                       chains = 6,iter = 600, warmup =300, adapt_delta = 0.9,max_treedepth = 10,
                                       prior_eta = 1, prior_rho = c(5.566019, 5.570317), prior_theta = 1, prior_sigma = 1, take_logs = FALSE, 
                                       scaled_lower_bound = NULL, l_bound = 0, mu=0,
                                       stanFile = "gp_compound_mean_03.stan", use_vb = FALSE){
  
  if(is.null(loss_square)) 
    loss_square = matrix(NA, nrow(loss_triangle), ncol(loss_triangle))
  
  AccidentYear_unique = AccidentYear = as.numeric(rownames(loss_triangle))
  DevelopmentLag_unique = DevelopmentLag = as.numeric(colnames(loss_triangle))
  DevelopmentLag = DevelopmentLag[-1]
  no_acc = length(AccidentYear); 
  no_dev = length(DevelopmentLag)
  AccidentYear = rep(AccidentYear, each = no_dev)
  DevelopmentLag = rep(DevelopmentLag, times = no_acc)
  
  loss_ratio_triangle = (loss_triangle / premiums)
  loss_ratio_triangle = loss_ratio_triangle[,-1] - loss_ratio_triangle[,-ncol(loss_ratio_triangle)]
  loss_ratio_triangle_log = log(apply((loss_triangle / premiums),2,function(x) pmax(x,eps)))
  loss_ratio_triangle_log = loss_ratio_triangle_log[,-1] - loss_ratio_triangle_log[,-ncol(loss_ratio_triangle_log)]
  
  
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
  #y_sd = sd(loss_ratio$loss_ratio)
  if (!take_logs) {
    l_bound = scaler_loss_ratio(l_bound)[1,1]
    mu = scaler_loss_ratio(mu)[1,1]
  } else {
    l_bound = scaler_loss_ratio_log(l_bound)[1,1]
    mu = scaler_loss_ratio_log(mu)[1,1]
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
    N1 = nrow(loss_ratio),
    N2 = nrow(loss_ratio_new),
    x1 =  model.matrix(formula_x, loss_ratio),
    y1 = pmax(l_bound, c(loss_ratio[,all.vars(update(formula_x, .~0)) ][[1]])),
    dev_lag1 = loss_ratio$DevelopmentLag,
    x2 = model.matrix(formula_x, loss_ratio_new),
    dev_lag2 = loss_ratio_new$DevelopmentLag,
    mu = mu,
    l_bound=l_bound,
    prior_eta=prior_eta,
    prior_rho = prior_rho,
    prior_theta = prior_theta,
    prior_sigma = prior_sigma,
    dimSigma = 11-2 
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
    stanRun <- sampling(object = stanMod, data = stan_list, iter=iter, warmup=warmup, chains=chains,init=init_list,
                        control = list(adapt_delta=adapt_delta, max_treedepth =max_treedepth) )
    if (plot_graphs) {
      print(stan_trace(stanRun,pars = c("rho","sigma"))); 
      #print(stan_trace(stanRun,pars = c("beta1", "beta2"))); 
      print(stan_dens(stanRun, pars = "yStar[90]"))
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
    return(list(loss_square=loss_square_pred,loss_square_med=loss_square_pred_med, loss_square_2=loss_square_pred_2, 
                rank_val_incremental=rank_val_incremental,rank_val_cumulative= rank_val_cumulative, sigma = sigma))
  }
}




#' Calculate the predicted cumulative claims using the LDF factor model
#' 
#' @param loss_triangle the Loss triangle and loss square
#' @param premiums Vector of premiums
#' @param calc_ranks If TRUE then additional statistics are calculated
#' @param formula_x formula for square exponential 
#' @param formula_H formula for linear terms in kernel (not in mean function)
#' @param take_logs If TRUE, then the exponential transform of the LDF is used
#' @param plot_graphs If TRUE, then sadmple investigative graphs are plotted.
#' @return If calc_ranks = FALSE, then only the predicted loss_square is returned.  If TRUE, then the percentage ranks, sigma estimates are also returned.

stan_ldf = function(loss_triangle, premiums=NULL, loss_square = NULL, calc_ranks = FALSE, 
                    formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                    formula_H = y_scaled ~  AccidentYear_scaled + DevelopmentLag_scaled -1 , formula_cat = NULL,
                    eps = 0.001, eps_log = 0.001,plot_graphs = FALSE,
                    chains = 6,iter = 600, warmup =300, adapt_delta = 0.9,max_treedepth = 10,
                    prior_eta = 1, prior_rho = c(5.566019, 5.570317), prior_theta = 1, prior_sigma = 1, 
                    take_logs = TRUE, 
                    scaled_lower_bound = NULL,
                    stanFile = "gp_compound_mean_03.stan", use_vb = FALSE){
  
  
  
  if(is.null(loss_square)) 
    loss_square = matrix(NA, nrow(loss_triangle), ncol(loss_triangle))
  
  AccidentYear_unique = AccidentYear = as.numeric(rownames(loss_triangle))
  DevelopmentLag_unique = DevelopmentLag = as.numeric(colnames(loss_triangle))
  DevelopmentLag = DevelopmentLag[-1]
  no_acc = length(AccidentYear); 
  no_dev = length(DevelopmentLag)
  AccidentYear = rep(AccidentYear, each = no_dev)
  DevelopmentLag = rep(DevelopmentLag, times = no_acc)
  
  ldf_triangle = loss_triangle[,-1] / loss_triangle[,-ncol(loss_triangle)]
  
  ldf_transform_fun = function(dev) {
    log(pmax(dev,1)-1+eps_log)
  }
  
  ldf_inv_transform_fun = function(dev) {
    exp(dev)+1-eps_log
  }
  
  
  
  ldf = as_tibble(ldf_triangle) %>%
    mutate(AccidentYear = as.integer(rownames(ldf_triangle))) %>%
    gather(DevelopmentLag, ldf, -AccidentYear) %>%
    filter(!is.na(ldf)) %>%
    mutate(
      DevelopmentLag = as.integer(DevelopmentLag),
      ldf_log = ldf_transform_fun(ldf)
    )
  
  scaler_accident_year = scaler(ldf$AccidentYear)
  scaler_development_lag = scaler(ldf$DevelopmentLag)
  DevelopmentLag_transform = function(x) log(x-min(ldf$DevelopmentLag)+1)
  scaler_ldf = scaler(ldf$ldf, l_bound = scaled_lower_bound)
  scaler_ldf_log = scaler(ldf$ldf_log, l_bound = scaled_lower_bound)
  #y_sd = sd(ldf$ldf)
  if (!take_logs) {
    l_bound = scaler_ldf(1)[1,1]
  } else {
    l_bound = scaler_ldf_log(ldf_transform_fun(1))[1,1]
  }
  ldf = ldf %>%
    mutate(AccidentYear_scaled = scaler_accident_year(AccidentYear)[,1],
           DevelopmentLag_scaled = scaler_development_lag(DevelopmentLag)[,1], 
           DevelopmentLag_trans = DevelopmentLag_transform(DevelopmentLag),
           y_scaled = scaler_ldf(ldf)[,1],
           y_log_scaled = scaler_ldf_log(ldf_log)[,1]) 
  
  
  # Create inputs for out-of-sample predictions
  ldf_new = tibble(AccidentYear ,DevelopmentLag) %>%
    mutate(AccidentYear_scaled = scaler_accident_year(AccidentYear)[,1],
           DevelopmentLag_scaled = scaler_development_lag(DevelopmentLag)[,1], 
           DevelopmentLag_trans = DevelopmentLag_transform(DevelopmentLag),
           y_scaled = 0,
           y_log_scaled = 0)
  
  formula_x = update(formula_x, . ~ . -1) # remove intercept term
  formula_H = update(formula_H, . ~ . -1) # remove intercept term
  
  stan_list <- list(
    N1 = nrow(ldf),
    N2 = nrow(ldf_new),
    x1 =  model.matrix(formula_x, ldf),
    y1 = pmax(l_bound, c(ldf[,all.vars(update(formula_x, .~0)) ][[1]])),
    dev_lag1 = ldf$DevelopmentLag-1, #index adjusted to start at 1
    x2 = model.matrix(formula_x, ldf_new),
    dev_lag2 = ldf_new$DevelopmentLag-1, #index adjusted to start at 1
    mu = l_bound,
    l_bound=l_bound,
    prior_eta=prior_eta,
    prior_rho = prior_rho,
    prior_theta = prior_theta,
    prior_sigma = prior_sigma,
    dimSigma = 11-2 
  )
  
  if (take_logs & (all.vars(formula_x)[1] != "y_log_scaled")){
    warning("if logs are used, LHS of formula_x must be 'y_log_scaled'")
  } 
  
  stan_list$D = ncol(stan_list$x1)
  
  if(!is.null(formula_H)){
    stan_list$H1 = t(model.matrix(formula_H, ldf))
    stan_list$H2 = t(model.matrix(formula_H, ldf_new))
    stan_list$H1 = stan_list$H1 - apply(stan_list$H1,1,min)
    stan_list$H2 = stan_list$H2 - apply(stan_list$H1,1,min)
    stan_list$dimH <- nrow(stan_list$H1)
  }
  
  
  if(!is.null(formula_cat)){
    formula_cat = update(formula_cat,  ~ . -1) # remove intercept term
    stan_list$x_cat1 = c(model.matrix(formula_cat, ldf))
    stan_list$x_cat2 = c(model.matrix(formula_cat, ldf_new))
  }
  stanMod <- stan_model(file=paste0(stan_dir, stanFile))
  
  
  if (!use_vb){
    stanRun <- sampling(object = stanMod, data = stan_list, iter=iter, warmup=warmup, chains=chains,
                        control = list(adapt_delta=adapt_delta, max_treedepth =max_treedepth) )
    
    summary_stats = summary(stanRun, pars = c("eta", "rho", "sigma"))$summary
    #stan_trace(stanRun)
    #plot(stanRun, pars="yStar")
  } else{
    stanRun <- vb(object = stanMod, data = stan_list )
  }

  
  ### old method :  delete?
  if (F){
    sigma = rstan::extract(stanRun, pars= "sigma")[[1]]
    sigma = colMeans(sigma)
    
    if (!take_logs) {
      preds = apply(rstan::extract(stanRun, pars= "yStar")[[1]], 2, FUN = scaler_ldf, unscale = T)
    } else {
      preds = apply(rstan::extract(stanRun, pars= "yStar")[[1]], 2, FUN = scaler_ldf_log, unscale = T)
    }
    
    preds = cbind(ldf_new[, c("AccidentYear","DevelopmentLag")],t(preds))
    preds = as.data.frame(preds) %>%
      arrange(AccidentYear, DevelopmentLag)
    
    preds = preds %>%
      left_join(ldf %>% select(AccidentYear, DevelopmentLag, ldf, ldf_log), by = c("AccidentYear", "DevelopmentLag"))
    #preds[1:15, 1:5]
    
    ldf_train = preds$ldf #makes sure that ldf and preds are in the same order
    preds$ldf = NULL
    preds$ldf_log = NULL
    
    # Calculate the percentage ranks for incremental loss ratios
    
    ldf_increment = loss_square[,-1] / loss_square[,-ncol(loss_square)]
    ldf_increment = c(t(ldf_increment)) #actual incremental ldf factors for full square
    if (take_logs){
      ldf_increment = ldf_transform_fun(ldf_increment)
    }
    # Calculate rates on LDF factors
    rank_val_incremental = vector_rank(v = ldf_increment, m = as.matrix(preds[,-1:-2]))
    
    
    actual_claims = c(t(loss_square[,-1])) #actual cumulative paid losses for full square
    
    # Calculate the cumulative loss ratios
    preds_cum = as.matrix(preds[,-1:-2])
    
    #dim(preds_cum)
    
    # replaces the predicted loss ratio with the actual for the training period and calculate cumulative claims
    preds_cum[!is.na(ldf_train), ] = ldf_increment[!is.na(ldf_train)] #!!!!!  CHECK!!!!
    #preds[1:5, 1:5]
    #preds_cum[1:20, 1:5]
    #rowMeans(preds_cum) - ldf_train
    
    # Reverse log transform for cumulative cases
    if (take_logs){
      preds_cum = ldf_inv_transform_fun(preds_cum) 
    }
    
    for (i in 1:length(AccidentYear)){
      if (DevelopmentLag[i] > 2)  {
        preds_cum[i, ] = preds_cum[i, ] * preds_cum[i-1, ]
      } 
    }
    #preds_cum[1:20,1:5]
    preds_cum = preds_cum *  rep(loss_triangle[,1], each= no_dev) # multiply back initial year claims
    
    #round(rowMeans(preds_cum) - actual_claims,2)
    
    # Calculate the percentage ranks for cumulative loss ratios
    rank_val_cumulative = vector_rank(actual_claims, preds_cum[,-1:-2])
    #plot(rank_val_cumulative)
    
    
    
    ldf_new = ldf_new %>%
      mutate(
        predicted_claims = rowMeans(preds_cum) , #colMeans(preds_cum), 
        lower_claims = apply(preds_cum,1,quantile,probs = 0.05), 
        upper_claims = apply(preds_cum,1,quantile,probs = 0.95),
        actual = actual_claims,
        predicted_ldf_inc = rowMeans(preds[,-1:-2]) , #colMeans(preds_cum), 
        predicted_ldf_median=apply(preds[,-1:-2],1,quantile,probs = 0.5), 
        lower_ldf_inc=apply(preds[,-1:-2],1,quantile,probs = 0.05), 
        upper_ldf_inc=apply(preds[,-1:-2],1,quantile,probs = 0.95),
        actual_ldf = ldf_increment
      )
    
    # Calculate Loss Square
    loss_square_pred = cbind(loss_triangle[,1],t(matrix(ldf_new$predicted_claims, no_dev, no_acc)))
    
    # Calculate mean/median log version
    pred_ldf2 = ldf_new$predicted_ldf_inc # mean LDF
    pred_ldf_med = ldf_new$predicted_ldf_median # mean LDF
    if (take_logs){
      pred_ldf2 = ldf_inv_transform_fun(pred_ldf2)
      pred_ldf_med = ldf_inv_transform_fun(pred_ldf_med)
    }
    
    pred_ldf2[!is.na(ldf_train)] = ldf_train[!is.na(ldf_train)] # Replace training sims with actual ldf
    loss_square_pred_2 = cbind(loss_triangle[,1],t(matrix(pred_ldf2, no_dev, no_acc)))
    loss_square_pred_2 = t(apply(loss_square_pred_2,1,cumprod))
    
    pred_ldf_med[!is.na(ldf_train)] = ldf_train[!is.na(ldf_train)] # Replace training sims with actual ldf
    loss_square_pred_med = cbind(loss_triangle[,1],t(matrix(pred_ldf_med, no_dev, no_acc)))
    loss_square_pred_med = t(apply(loss_square_pred_med,1,cumprod))
    loss_square_pred_2 / loss_square_pred_med
    
    
    # if (!take_logs){
    #   predicted_loss_ratios = ldf_new$predicted_ldf_inc
    #   loss_square_pred = cbind(loss_triangle[,1],t(matrix(predicted_loss_ratios, no_dev,no_acc))* premiums)
    #   loss_square_pred_med = NA
    #   loss_square_pred_2 = NA
    # } else {
    #   predicted_loss_ratios = rowMeans(exp(preds_cum))
    #   predicted_loss_ratios_med = apply(exp(preds_cum),1,median)
    #   predicted_loss_ratios_2 = exp(rowMeans(preds_cum))
    #   loss_square_pred = cbind(loss_triangle[,1],t(matrix(predicted_loss_ratios, no_dev,no_acc))* premiums)
    #   loss_square_pred_med = cbind(loss_triangle[,1],t(matrix(predicted_loss_ratios_med, no_dev,no_acc))* premiums)
    #   loss_square_pred_2 = cbind(loss_triangle[,1],t(matrix(predicted_loss_ratios_2, no_dev,no_acc))* premiums)
    # }
    
    if (plot_graphs){
      gplot = ggplot(ldf_new) + 
        geom_vline(aes(xintercept = 1998-AccidentYear), color = "green", linetype="dashed")+
        geom_line(aes(DevelopmentLag, actual,color="actual"))+ 
        geom_line(aes(DevelopmentLag, predicted_claims,color="pred"))+
        geom_point(aes(DevelopmentLag, actual,color="actual"), size=.5)+ 
        geom_point(aes(DevelopmentLag, predicted_claims,color="pred"), size=.5)+
        geom_ribbon(aes(DevelopmentLag, ymin=lower_claims, ymax=upper_claims), color= "grey", alpha=0.2)+
        facet_wrap(~AccidentYear, scales="free")+ 
        labs(title = "Cumulative Paid Losses: Actual vs. Predicted\nusing LDF Model") 
      print(gplot)
      
      gplot = ggplot(ldf_new) + 
        geom_vline(aes(xintercept = 1998-AccidentYear), color = "green", linetype="dashed")+
        geom_line(aes(DevelopmentLag, actual_ldf,color="actual"))+ 
        geom_line(aes(DevelopmentLag, predicted_ldf_inc,color="pred"))+
        geom_point(aes(DevelopmentLag, actual_ldf,color="actual"), size=.5)+ 
        geom_point(aes(DevelopmentLag, predicted_ldf_inc,color="pred"), size=.5)+
        geom_ribbon(aes(DevelopmentLag, ymin=lower_ldf_inc, ymax=upper_ldf_inc), color= "grey", alpha=0.2)+
        facet_wrap(~AccidentYear, scales="free")+ 
        labs(title = "Loss Development Factors: Actual vs. Predicted\nusing LDF Model") 
      print(gplot)
      
      print(stan_plot(stanRun,pars= "sigma"))
      
    }
    
    if (!calc_ranks){
      return(loss_square = loss_square_pred)
    } else {
      return(list(loss_square=loss_square_pred,loss_square=loss_square_pred, loss_square_2=loss_square_pred_2,
                  loss_square_med=loss_square_pred_med, 
                  rank_val_incremental=rank_val_incremental,rank_val_cumulative= rank_val_cumulative, sigma = sigma, 
                  pred_LDF = apply(preds,2, ldf_inv_transform_fun)))
    }
  }
  
  #################
  # Extract and unscale predictions
  sims = rstan::extract(stanRun, pars= "yStar")[[1]]
  if (!take_logs) {
    sims =apply(sims, 2, FUN = scaler_ldf, unscale = T)
  } else {
    sims =apply(sims, 2, FUN = scaler_ldf_log, unscale=T)
    sims =apply(sims, 2, FUN = ldf_inv_transform_fun)
  }
  
  # Calc simulated  losses
  results = ldf_sims_to_losses(sims = sims, loss_triangle = loss_triangle)
  
  #preds$dev_preds = exp(dev_log_preds$mean * sd_dev_log + mean_dev_log)+1-eps
  
  
  # Plot Actual vs. Predicted Factors
  if (plot_graphs){
    ldf_new$dev_preds = c(t(dev_factors(results$loss_square)))
    ldf_new$predicted = c(t((results$loss_square[,-1])))
    
    ldf_new$actual = c(t(loss_triangle[,-1]))
    
    dl_index = rep(DevelopmentLag_unique, times = length(AccidentYear_unique))
    sims_lag2_plus = results$sim_losses[, dl_index !=1]
    
    ldf_new$lower_claims = apply(sims_lag2_plus, 2, FUN = quantile, probs=0.05)
    ldf_new$upper_claims = apply(sims_lag2_plus, 2, FUN = quantile, probs=0.95)
    
    print(ggplot(cbind(ldf_new, intercept=1)) + 
            theme_bw() +
            geom_point(aes(DevelopmentLag, actual, color="Actual"))+
            geom_line(aes(DevelopmentLag, predicted, color="Predicted")) + 
            facet_wrap(~p0("Acc. Yr ",AccidentYear), scales="free") + 
            geom_ribbon(aes(DevelopmentLag, ymin=lower_claims, ymax=upper_claims), color= "grey", alpha=0.2)+
            #geom_hline(aes(yintercept = intercept), size = 0.25, linetype = 2) + 
            labs(title = p0("GP Fit to Log Loss Development Factors\nformula = ", format(formula) )))
    
    # print(ggplot(cbind(ldf_new, intercept=1)) + 
    #         theme_bw() +
    #         geom_point(aes(DevelopmentLag, dev, color="Actual"))+
    #         geom_line(aes(DevelopmentLag, dev_preds, color="Predicted")) + 
    #         facet_wrap(~p0("Acc. Yr ",AccidentYear)) + 
    #         geom_hline(aes(yintercept = intercept), size = 0.25, linetype = 2) + 
    #         labs(title = p0("GP Fit to Log Loss Development Factors\nformula = ", format(formula) )))
    # 
    # print(ggplot(cbind(preds, intercept=1)) + 
    #         theme_bw() +
    #         geom_point(aes(AccidentYear, dev, color="Actual"))+
    #         geom_line(aes(AccidentYear, dev_preds, color="Predicted")) + 
    #         facet_wrap(~p0("Dev. Yr ",sprintf("%02d", DevelopmentLag))) + 
    #         geom_hline(aes(yintercept = intercept), size = 0.25, linetype = 2) + 
    #         labs(title = p0("GP Fit to Log Loss Development Factors\nformula = ", format(formula) )))
  }
  
  list(loss_square =results$loss_square, sim_ult_losses = results$sim_ult_losses, sim_losses=results$sim_losses, summary_stats=summary_stats) 
  
  
  
} 



#' Calculate the predicted cumulative claims using the LDF wiht HURDLE factor model
#' 
#' @param loss_triangle the Loss triangle and loss square
#' @param premiums Vector of premiums
#' @param calc_ranks If TRUE then additional statistics are calculated
#' @param formula_x formula for square exponential 
#' @param formula_H formula for linear terms in kernel (not in mean function)
#' @param take_logs If TRUE, then the exponential transform of the LDF is used
#' @param plot_graphs If TRUE, then sample investigative graphs are plotted.
#' @return If calc_ranks = FALSE, then only the predicted loss_square is returned.  If TRUE, then the percentage ranks, sigma estimates are also returned.
stan_ldf2 = function(loss_triangle, premiums=NULL, loss_square = NULL, calc_ranks = FALSE, 
                     formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                     formula_H = y_scaled ~  AccidentYear_scaled + DevelopmentLag_scaled -1 , formula_cat = NULL,
                     eps = 0.001, eps_log = 0.001,plot_graphs = FALSE,
                     chains = 6,iter = 600, warmup =300, adapt_delta = 0.9,max_treedepth = 10,
                     prior_eta = 1, prior_rho = c(5.566019, 5.570317), prior_theta = 1, prior_sigma = 1, take_logs = FALSE, 
                     l_bound =1, mu=9:1/9+.5,
                     scaled_lower_bound = NULL,
                     stanFile = "gp_hurdle_02.stan", use_vb = FALSE){
  
  
  if(is.null(loss_square)) 
    loss_square = matrix(NA, nrow(loss_triangle), ncol(loss_triangle))
  
  AccidentYear_unique = AccidentYear = as.numeric(rownames(loss_triangle))
  DevelopmentLag_unique = DevelopmentLag = as.numeric(colnames(loss_triangle))
  DevelopmentLag = DevelopmentLag[-1]
  no_acc = length(AccidentYear); 
  no_dev = length(DevelopmentLag)
  AccidentYear = rep(AccidentYear, each = no_dev)
  DevelopmentLag = rep(DevelopmentLag, times = no_acc)
  
  ldf_triangle = loss_triangle[,-1] / loss_triangle[,-ncol(loss_triangle)]
  
  ldf_transform_fun = function(dev) {
    log(pmax(dev,1)-1+eps_log)
  }
  
  ldf_inv_transform_fun = function(dev) {
    exp(dev)+1-eps_log
  }
  
  
  
  ldf = as_tibble(ldf_triangle) %>%
    mutate(AccidentYear = as.integer(rownames(ldf_triangle))) %>%
    gather(DevelopmentLag, ldf, -AccidentYear) %>%
    filter(!is.na(ldf)) %>%
    mutate(
      DevelopmentLag = as.integer(DevelopmentLag),
      ldf_log = ldf_transform_fun(ldf)
    )
  
  scaler_accident_year = scaler(ldf$AccidentYear)
  scaler_development_lag = scaler(ldf$DevelopmentLag)
  DevelopmentLag_transform = function(x) log(x-min(ldf$DevelopmentLag)+1)
  scaler_ldf = scaler(ldf$ldf, l_bound = scaled_lower_bound)
  scaler_ldf_log = scaler(ldf$ldf_log, l_bound = scaled_lower_bound)
  #y_sd = sd(ldf$ldf)
  if (!take_logs) {
    l_bound = scaler_ldf(l_bound)[1,1]
    mu = c(scaler_ldf(mu))
    
  } else {
    l_bound = scaler_ldf_log(ldf_transform_fun(l_bound))[1,1]
    mu = c(scaler_ldf_log(ldf_transform_fun(mu))) 
  }
  ldf = ldf %>%
    mutate(AccidentYear_scaled = scaler_accident_year(AccidentYear)[,1],
           DevelopmentLag_scaled = scaler_development_lag(DevelopmentLag)[,1], 
           DevelopmentLag_trans = DevelopmentLag_transform(DevelopmentLag),
           y_scaled = scaler_ldf(ldf)[,1],
           y_log_scaled = scaler_ldf_log(ldf_log)[,1]) 
  
  
  # Create inputs for out-of-sample predictions
  ldf_new = tibble(AccidentYear ,DevelopmentLag) %>%
    mutate(AccidentYear_scaled = scaler_accident_year(AccidentYear)[,1],
           DevelopmentLag_scaled = scaler_development_lag(DevelopmentLag)[,1], 
           DevelopmentLag_trans = DevelopmentLag_transform(DevelopmentLag),
           y_scaled = 0,
           y_log_scaled = 0)
  
  formula_x = update(formula_x, . ~ . -1) # remove intercept term
  formula_H = update(formula_H, . ~ . -1) # remove intercept term
  
  stan_list <- list(
    N1 = nrow(ldf),
    N2 = nrow(ldf_new),
    x1 =  model.matrix(formula_x, ldf),
    y1 = pmax(l_bound, c(ldf[,all.vars(update(formula_x, .~0)) ][[1]])),
    dev_lag1 = ldf$DevelopmentLag,
    x2 = model.matrix(formula_x, ldf_new),
    dev_lag2 = ldf_new$DevelopmentLag,
    mu = mu,
    l_bound=l_bound,
    prior_eta=prior_eta,
    prior_rho = prior_rho,
    prior_theta = prior_theta,
    prior_sigma = prior_sigma,
    dimSigma = 11-2 
  )
  
  if (take_logs & (all.vars(formula_x)[1] != "y_log_scaled")){
    warning("if logs are used, LHS of formula_x must be 'y_log_scaled'")
  } 
  
  stan_list$D = ncol(stan_list$x1)
  
  if(!is.null(formula_H)){
    stan_list$H1 = t(model.matrix(formula_H, ldf))
    stan_list$H2 = t(model.matrix(formula_H, ldf_new))
    stan_list$H1 = stan_list$H1 - apply(stan_list$H1,1,min)
    stan_list$H2 = stan_list$H2 - apply(stan_list$H1,1,min)
    stan_list$dimH <- nrow(stan_list$H1)
  }
  
  
  if(!is.null(formula_cat)){
    formula_cat = update(formula_cat,  ~ . -1) # remove intercept term
    stan_list$x_cat1 = c(model.matrix(formula_cat, ldf))
    stan_list$x_cat2 = c(model.matrix(formula_cat, ldf_new))
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
                        control = list(adapt_delta=adapt_delta, max_treedepth =max_treedepth) )
    
    
    if (plot_graphs) {
      print(stan_trace(stanRun,pars = c("rho","sigma"))); 
      print(stan_dens(stanRun, pars = "yStar[90]"))
    }
  } else{
    stanRun <- vb(object = stanMod, data = stan_list )
  }
  
  sigma = rstan::extract(stanRun, pars= "sigma")[[1]]
  sigma = colMeans(sigma)
  
  if (!take_logs) {
    preds = apply(rstan::extract(stanRun, pars= "yStar")[[1]], 2, FUN = scaler_ldf, unscale = T)
  } else {
    preds = apply(rstan::extract(stanRun, pars= "yStar")[[1]], 2, FUN = scaler_ldf_log, unscale = T)
    preds = apply(preds, 2, FUN = ldf_inv_transform_fun)
  }
  
  preds = cbind(ldf_new[, c("AccidentYear","DevelopmentLag")],t(preds))
  preds = as.data.frame(preds) %>%
    arrange(AccidentYear, DevelopmentLag)
  
  preds = preds %>%
    left_join(ldf %>% select(AccidentYear, DevelopmentLag, ldf, ldf_log), by = c("AccidentYear", "DevelopmentLag"))
  #preds[1:15, 1:5]
  
  ldf_train = preds$ldf #makes sure that ldf and preds are in the same order
  preds$ldf = NULL
  preds$ldf_log = NULL
  
  # Calculate the percentage ranks for incremental loss ratios
  
  ldf_increment = loss_square[,-1] / loss_square[,-ncol(loss_square)]
  ldf_increment = c(t(ldf_increment)) #actual incremental ldf factors for full square
  if (take_logs){
    ldf_increment = ldf_transform_fun(ldf_increment)
  }
  # Calculate rates on LDF factors
  rank_val_incremental = vector_rank(v = ldf_increment, m = as.matrix(preds[,-1:-2]))
  
  
  actual_claims = c(t(loss_square[,-1])) #actual cumulative paid losses for full square
  
  # Calculate the cumulative loss ratios
  preds_cum = as.matrix(preds[,-1:-2])
  
  #dim(preds_cum)
  
  # replaces the predicted loss ratio with the actual for the training period and calculate cumulative claims
  preds_cum[!is.na(ldf_train), ] = ldf_increment[!is.na(ldf_train)] #!!!!!  CHECK!!!!
  #preds[1:5, 1:5]
  #preds_cum[1:20, 1:5]
  #rowMeans(preds_cum) - ldf_train
  
  # Reverse log transform for cumulative cases
  if (take_logs){
    preds_cum = ldf_inv_transform_fun(preds_cum) 
  }
  
  for (i in 1:length(AccidentYear)){
    if (DevelopmentLag[i] > 2)  {
      preds_cum[i, ] = preds_cum[i, ] * preds_cum[i-1, ]
    } 
  }
  #preds_cum[1:20,1:5]
  preds_cum = preds_cum *  rep(loss_triangle[,1], each= no_dev) # multiply back initial year claims
  
  #round(rowMeans(preds_cum) - actual_claims,2)
  
  # Calculate the percentage ranks for cumulative loss ratios
  rank_val_cumulative = vector_rank(actual_claims, preds_cum[,-1:-2])
  #plot(rank_val_cumulative)
  
  
  
  ldf_new = ldf_new %>%
    mutate(
      predicted_claims = rowMeans(preds_cum) , #colMeans(preds_cum), 
      lower_claims = apply(preds_cum,1,quantile,probs = 0.05), 
      upper_claims = apply(preds_cum,1,quantile,probs = 0.95),
      actual = actual_claims,
      predicted_ldf_inc = rowMeans(preds[,-1:-2]) , #colMeans(preds_cum), 
      predicted_ldf_median=apply(preds[,-1:-2],1,quantile,probs = 0.5), 
      lower_ldf_inc=apply(preds[,-1:-2],1,quantile,probs = 0.05), 
      upper_ldf_inc=apply(preds[,-1:-2],1,quantile,probs = 0.95),
      actual_ldf = ldf_increment
    )
  
  # Calculate Loss Square
  loss_square_pred = cbind(loss_triangle[,1],t(matrix(ldf_new$predicted_claims, no_dev, no_acc)))
  
  # Calculate mean/median log version
  pred_ldf2 = ldf_new$predicted_ldf_inc # mean LDF
  pred_ldf_med = ldf_new$predicted_ldf_median # mean LDF
  if (take_logs){
    pred_ldf2 = ldf_inv_transform_fun(pred_ldf2)
    pred_ldf_med = ldf_inv_transform_fun(pred_ldf_med)
  }
  
  pred_ldf2[!is.na(ldf_train)] = ldf_train[!is.na(ldf_train)] # Replace training sims with actual ldf
  loss_square_pred_2 = cbind(loss_triangle[,1],t(matrix(pred_ldf2, no_dev, no_acc)))
  loss_square_pred_2 = t(apply(loss_square_pred_2,1,cumprod))
  
  pred_ldf_med[!is.na(ldf_train)] = ldf_train[!is.na(ldf_train)] # Replace training sims with actual ldf
  loss_square_pred_med = cbind(loss_triangle[,1],t(matrix(pred_ldf_med, no_dev, no_acc)))
  loss_square_pred_med = t(apply(loss_square_pred_med,1,cumprod))
  loss_square_pred_2 / loss_square_pred_med
  
  
  # if (!take_logs){
  #   predicted_loss_ratios = ldf_new$predicted_ldf_inc
  #   loss_square_pred = cbind(loss_triangle[,1],t(matrix(predicted_loss_ratios, no_dev,no_acc))* premiums)
  #   loss_square_pred_med = NA
  #   loss_square_pred_2 = NA
  # } else {
  #   predicted_loss_ratios = rowMeans(exp(preds_cum))
  #   predicted_loss_ratios_med = apply(exp(preds_cum),1,median)
  #   predicted_loss_ratios_2 = exp(rowMeans(preds_cum))
  #   loss_square_pred = cbind(loss_triangle[,1],t(matrix(predicted_loss_ratios, no_dev,no_acc))* premiums)
  #   loss_square_pred_med = cbind(loss_triangle[,1],t(matrix(predicted_loss_ratios_med, no_dev,no_acc))* premiums)
  #   loss_square_pred_2 = cbind(loss_triangle[,1],t(matrix(predicted_loss_ratios_2, no_dev,no_acc))* premiums)
  # }
  
  if (plot_graphs){
    gplot = ggplot(ldf_new) + 
      geom_vline(aes(xintercept = 1998-AccidentYear), color = "green", linetype="dashed")+
      geom_line(aes(DevelopmentLag, actual,color="actual"))+ 
      geom_line(aes(DevelopmentLag, predicted_claims,color="pred"))+
      geom_point(aes(DevelopmentLag, actual,color="actual"), size=.5)+ 
      geom_point(aes(DevelopmentLag, predicted_claims,color="pred"), size=.5)+
      geom_ribbon(aes(DevelopmentLag, ymin=lower_claims, ymax=upper_claims), color= "grey", alpha=0.2)+
      facet_wrap(~AccidentYear, scales="free")+ 
      labs(title = "Cumulative Paid Losses: Actual vs. Predicted\nusing LDF Model") 
    print(gplot)
    
    gplot = ggplot(ldf_new) + 
      geom_vline(aes(xintercept = 1998-AccidentYear), color = "green", linetype="dashed")+
      geom_hline(aes(yintercept = 1), color = "orange", linetype="dashed") +
      geom_line(aes(DevelopmentLag, actual_ldf,color="actual"))+ 
      geom_line(aes(DevelopmentLag, predicted_ldf_inc,color="pred"))+
      geom_point(aes(DevelopmentLag, actual_ldf,color="actual"), size=.5)+ 
      geom_point(aes(DevelopmentLag, predicted_ldf_inc,color="pred"), size=.5)+
      geom_ribbon(aes(DevelopmentLag, ymin=lower_ldf_inc, ymax=upper_ldf_inc), color= "grey", alpha=0.2)+
      facet_wrap(~AccidentYear, scales="free")+ 
      labs(title = "Loss Development Factors: Actual vs. Predicted\nusing LDF Model") 
    print(gplot)
    
    print(stan_plot(stanRun,pars= "sigma"))
    
  }
  
  if (!calc_ranks){
    return(loss_square = loss_square_pred)
  } else {
    
    
    errors = error_calcs(premiums = premiums, predicted = loss_square_pred, actual = loss_square)
    
    errors$rmse_mean = rmse_triangle(actual_sq = loss_square/premiums, predicted_sq = loss_square_pred/    premiums,actual_tr = loss_triangle/premiums)
    errors$rmse_med =  rmse_triangle(actual_sq = loss_square/premiums, predicted_sq = loss_square_pred_med/premiums,actual_tr = loss_triangle/premiums)
    
    
    
    return(list(loss_square=loss_square_pred,loss_square=loss_square_pred, loss_square_2=loss_square_pred_2,
                loss_square_med=loss_square_pred_med, 
                rank_val_incremental=rank_val_incremental,rank_val_cumulative= rank_val_cumulative, sigma = sigma, 
                pred_LDF = apply(preds,2, ldf_inv_transform_fun), 
                errors = errors))
  }
} 









#' TWO STEP Model: Calculate the predicted loss ratio and supporting analytics using an incremental loss ratio model
#' 
#' @param loss_triangle the Loss triangle and loss square
#' @param premiums Vector of premiums
#' @param calc_ranks If TRUE then additional statistics are calculated
#' @param formula_x formula for square exponential 
#' @param formula_H formula for linear terms in kernel (not in mean function)
#' @param take_logs If TRUE, then response is logarithm is used as the predictor.  Note, if TRUE, then LHS of formula_x must be y_log_scaled, otherwise use y_scaled
#' @param plot_graphs If TRUE, then sample investigative graphs are plotted.
#' @return If calc_ranks = FALSE, then only the predicted loss_square is returned.  If TRUE, then the percentage ranks, sigma estimates are also returned.
residual_mod_incremental_lr = function(loss_triangle, premiums=NULL, loss_square = NULL, calc_ranks = FALSE, 
                                       formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                       formula_H = y_scaled ~  AccidentYear_scaled + DevelopmentLag_scaled -1 , formula_cat = NULL,
                                       eps = 0.001, plot_graphs = FALSE,
                                       chains = 6,iter = 600, warmup =300, adapt_delta = 0.9,max_treedepth = 10,
                                       prior_eta = 1, prior_rho = c(5.566019, 5.570317), prior_theta = 1, prior_sigma = 1, take_logs = FALSE, 
                                       scaled_lower_bound = NULL,
                                       stanFile = "gp_compound_mean_03.stan", 
                                       stanFileResidual = "residual_01.stan",
                                       use_vb = FALSE){
  
  if(is.null(loss_square)) 
    loss_square = matrix(NA, nrow(loss_triangle), ncol(loss_triangle))
  
  AccidentYear_unique = AccidentYear = as.numeric(rownames(loss_triangle))
  DevelopmentLag_unique = DevelopmentLag = as.numeric(colnames(loss_triangle))
  DevelopmentLag = DevelopmentLag[-1]
  no_acc = length(AccidentYear); 
  no_dev = length(DevelopmentLag)
  AccidentYear = rep(AccidentYear, each = no_dev)
  DevelopmentLag = rep(DevelopmentLag, times = no_acc)
  
  loss_ratio_triangle = (loss_triangle / premiums)
  loss_ratio_triangle = loss_ratio_triangle[,-1] - loss_ratio_triangle[,-ncol(loss_ratio_triangle)]
  loss_ratio_triangle_log = log(apply((loss_triangle / premiums),2,function(x) pmax(x,eps)))
  loss_ratio_triangle_log = loss_ratio_triangle_log[,-1] - loss_ratio_triangle_log[,-ncol(loss_ratio_triangle_log)]
  
  
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
  #y_sd = sd(loss_ratio$loss_ratio)
  if (!take_logs) {
    l_bound = scaler_loss_ratio(0)[1,1]
  } else {
    l_bound = scaler_loss_ratio_log(0)[1,1]
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
    N1 = nrow(loss_ratio),
    N2 = nrow(loss_ratio_new),
    x1 =  model.matrix(formula_x, loss_ratio),
    y1 = pmax(l_bound, c(loss_ratio[,all.vars(update(formula_x, .~0)) ][[1]])),
    dev_lag1 = loss_ratio$DevelopmentLag,
    x2 = model.matrix(formula_x, loss_ratio_new),
    dev_lag2 = loss_ratio_new$DevelopmentLag,
    mu = l_bound,
    l_bound=l_bound,
    prior_eta=prior_eta,
    prior_rho = prior_rho,
    prior_theta = prior_theta,
    prior_sigma = prior_sigma,
    dimSigma = 11-2 
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
  
  stanMod <- stan_model(file=paste0(stan_dir, stanFile))
  
  
  if (!use_vb){
    stanRun <- sampling(object = stanMod, data = stan_list, iter=iter, warmup=warmup, chains=chains,
                        control = list(adapt_delta=adapt_delta, max_treedepth =max_treedepth) )
    #stan_trace(stanRun)
  } else{
    stanRun <- vb(object = stanMod, data = stan_list )
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
  #rank_val_incremental = vector_rank(v = lr_increment, m = as.matrix(preds[,-1:-2]))
  
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
  #rank_val_cumulative = vector_rank(lr_cumulative, preds_cum[,-1:-2])
  #plot(rank_val_cumulative)
  
  
  
  loss_ratio_new = loss_ratio_new %>%
    mutate(
      predicted_lr = rowMeans(preds_cum) , #colMeans(preds_cum), 
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
    loss_square_pred_med = NA
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
      geom_line(aes(DevelopmentLag, actual_inc,color="actual"))+ 
      geom_line(aes(DevelopmentLag, predicted_lr_inc,color="pred"))+
      geom_point(aes(DevelopmentLag, actual_inc,color="actual"), size=.5)+ 
      geom_point(aes(DevelopmentLag, predicted_lr_inc,color="pred"), size=.5)+
      geom_ribbon(aes(DevelopmentLag, ymin=lower_lr_inc, ymax=upper_lr_inc), color= "grey", alpha=0.2)+
      facet_wrap(~AccidentYear)+ 
      labs(title = "Incremental Loss Ratio: Actual vs. Predicted") 
    print(gplot)
    
    print(stan_plot(stanRun,pars= "sigma"))
    
  }
  
  #### Residual Model
  #fStar = rstan::extract(stanRun, pars = "fStar")[[1]]
  # fStar = preds %>%
  #   arrange(AccidentYear, DevelopmentLag ) %>%
  #   select(-AccidentYear, -DevelopmentLag) %>%
  #   as.matrix() %>%
  #   rowMeans()
  
  loss_ratio_new$training_set = loss_ratio_new$AccidentYear + loss_ratio_new$DevelopmentLag <= max(loss_ratio_new$AccidentYear)+1
  stan_list_residual = list(
    N1 = stan_list$N1,
    N2 = stan_list$N2,
    dimSigma = stan_list$dimSigma,
    dev_lag1 = loss_ratio_new$DevelopmentLag[loss_ratio_new$training_set],
    dev_lag2 = loss_ratio_new$DevelopmentLag,
    prior_sigma = prior_sigma,
    f1 = loss_ratio_new$predicted_lr_inc[loss_ratio_new$training_set],
    f2 = loss_ratio_new$predicted_lr_inc,
    y1 = loss_ratio_new$actual_inc[loss_ratio_new$training_set] 
  )
  
  stanModResidual = stan_model(file=paste0(stan_dir, stanFileResidual))
  a=Sys.time()
  stanRunResidual <- sampling(object = stanModResidual, data = stan_list_residual, iter=1200, warmup=200, chains=1,
                              control = list(adapt_delta=adapt_delta, max_treedepth =max_treedepth) )
  b=Sys.time();b-a
  # stanRunResidual
  # plot(stanRunResidual,pars="sigma")
  # plot(stanRunResidual,pars="yStar")
  
  preds_step2 = t(rstan::extract(stanRunResidual,pars = "yStar")[[1]])
  dim(preds_step2)
  preds_step2_cum = preds_step2
  #replaces the predicted loss ratio with the actual for the training period and calculate cumulative claims
  preds_step2_cum[!is.na(loss_ratio_inc_train), ] = loss_ratio_inc_train[!is.na(loss_ratio_inc_train)] 
  #preds_step2_cum[1:20, 1:5]
  #rowMeans(preds_step2_cum) - loss_ratio_inc_train
  
  for (i in 1:length(AccidentYear)){
    if (DevelopmentLag[i] > 2)  {
      preds_step2_cum[i, ] = preds_step2_cum[i, ] + preds_step2_cum[i-1, ]
    } 
  }
  preds_step2_cum[1:5,1:5]
  if (!take_logs){
    preds_step2_cum = preds_step2_cum +  rep(loss_triangle[,1]/premiums, each= no_dev) # add back initial year claims
  } else {
    preds_step2_cum = preds_step2_cum +  rep(log(pmax(eps,loss_triangle[,1]/premiums)), each= no_dev) # add back initial year claims
  }  
  
  # Incremental LR stats
  loss_ratio_new$pred2 = rowMeans(preds_step2)
  loss_ratio_new$lower_2 = apply(preds_step2, 1, quantile, probs=0.05)
  loss_ratio_new$upper_2 = apply(preds_step2, 1, quantile, probs=0.95)
  # Cumulative LR stats
  loss_ratio_new$pred2_cum = rowMeans(preds_step2_cum)
  loss_ratio_new$lower_2_cum = apply(preds_step2_cum, 1, quantile, probs=0.05)
  loss_ratio_new$upper_2_cum = apply(preds_step2_cum, 1, quantile, probs=0.95)
  
  if (plot_graphs){
    loss_ratio_new$DevelopmentLag =as.integer(loss_ratio_new$DevelopmentLag)
    gplot = ggplot(loss_ratio_new) + 
      geom_vline(aes(xintercept = 1998-AccidentYear), color = "green", linetype="dashed")+
      geom_line(aes(DevelopmentLag, actual,color="actual"))+ 
      geom_line(aes(DevelopmentLag, pred2_cum,color="pred"))+
      geom_point(aes(DevelopmentLag, actual,color="actual"), size=.5)+ 
      geom_point(aes(DevelopmentLag, pred2_cum,color="pred"), size=.5)+
      geom_ribbon(aes(DevelopmentLag, ymin=lower_2_cum, ymax=upper_2_cum), color= "grey", alpha=0.2)+
      facet_wrap(~AccidentYear)+ 
      labs(title = "Two Step Model: Cumulative Loss Ratio: Actual vs. Predicted") 
    print(gplot)
    
    gplot = ggplot(loss_ratio_new) + 
      geom_vline(aes(xintercept = 1998-AccidentYear), color = "green", linetype="dashed")+
      geom_hline(aes(yintercept = 0), color = "orange", linetype="dashed")+
      geom_line(aes(DevelopmentLag, actual_inc,color="actual"))+ 
      geom_line(aes(DevelopmentLag, pred2,color="pred"))+
      geom_point(aes(DevelopmentLag, actual_inc,color="actual"), size=.5)+ 
      geom_point(aes(DevelopmentLag, pred2,color="pred"), size=.5)+
      geom_ribbon(aes(DevelopmentLag, ymin=lower_2, ymax=upper_2), color= "grey", alpha=0.2)+
      facet_wrap(~AccidentYear)+ 
      labs(title = "TWO Step Model: Incremental Loss Ratio: Actual vs. Predicted") 
    print(gplot)
    
    print(stan_plot(stanRun,pars= "sigma"))
    
  }
  rank_val_incremental = vector_rank(v = loss_ratio_new$actual_inc, m = preds_step2)
  rank_val_cumulative = vector_rank(loss_ratio_new$actual, preds_step2_cum)

  if (!take_logs){
    loss_square_pred = cbind(loss_triangle[,1],t(matrix(loss_ratio_new$pred2_cum, no_dev,no_acc))* premiums) 
    loss_square_pred_med = NA
    loss_square_pred_2 = NA 
  } else {
    predicted_loss_ratios = rowMeans(exp(preds_step2_cum))
    predicted_loss_ratios_med = apply(exp(preds_step2_cum),1,median)
    predicted_loss_ratios_2 = exp(rowMeans(preds_step2_cum))
    loss_square_pred = cbind(loss_triangle[,1],t(matrix(loss_ratio_new$pred2_cum, no_dev,no_acc))* premiums) 
    loss_square_pred_med = cbind(loss_triangle[,1],t(matrix(predicted_loss_ratios_med, no_dev,no_acc))* premiums) 
    loss_square_pred_2 = cbind(loss_triangle[,1],t(matrix(predicted_loss_ratios_2, no_dev,no_acc))* premiums)
  }
  
  if (!calc_ranks){
    return(loss_square = loss_square_pred)
  } else {
    return(list(loss_square=loss_square_pred,loss_square_med=loss_square_pred_med, loss_square_2=loss_square_pred_2, 
                rank_val_incremental=rank_val_incremental,rank_val_cumulative= rank_val_cumulative, sigma = sigma))
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

stan_incremental_loss_ratio2 = function(loss_triangle, premiums=NULL, loss_square = NULL, calc_ranks = FALSE, verbose=TRUE,
                                        formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                        formula_H = y_scaled ~  AccidentYear_scaled + DevelopmentLag_scaled -1 , formula_cat = NULL,
                                        eps = 0.001, plot_graphs = FALSE,
                                        chains = 6,iter = 600, warmup =300, adapt_delta = 0.9,max_treedepth = 10,
                                        prior_eta = 1, prior_rho = c(5.566019, 5.570317), prior_theta = 1, prior_sigma = 1, 
                                        prior_beta = 1, take_logs = FALSE, 
                                        scaled_lower_bound = NULL, l_bound = 0, mu=9:1/9-.5,
                                        virtual_data = FALSE,
                                        stanFile = "gp_compound_mean_03.stan", use_vb = FALSE){
  

  if(is.null(loss_square)) 
    loss_square = matrix(NA, nrow(loss_triangle), ncol(loss_triangle))
  
  AccidentYear_unique = AccidentYear = as.numeric(rownames(loss_triangle))
  DevelopmentLag_unique = DevelopmentLag = as.numeric(colnames(loss_triangle))
  DevelopmentLag = DevelopmentLag[-1]
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
  loss_ratio_triangle = loss_ratio_triangle[,-1] - loss_ratio_triangle[,-ncol(loss_ratio_triangle)]
  loss_ratio_triangle_log = log(apply((loss_triangle / premiums),2,function(x) pmax(x,eps)))
  loss_ratio_triangle_log = loss_ratio_triangle_log[,-1] - loss_ratio_triangle_log[,-ncol(loss_ratio_triangle_log)]
  
  if (virtual_data){
    loss_ratio_triangle = cbind(loss_ratio_triangle, 0)
    loss_ratio_triangle_log = cbind(loss_ratio_triangle_log, 0)
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
  
  scaler_accident_year = scaler(loss_ratio$AccidentYear)
  scaler_development_lag = scaler(loss_ratio$DevelopmentLag)
  DevelopmentLag_transform = function(x) log(x-min(loss_ratio$DevelopmentLag)+1)
  scaler_loss_ratio = scaler(loss_ratio$loss_ratio, l_bound = scaled_lower_bound)
  scaler_loss_ratio_log = scaler(loss_ratio$loss_ratio_log, l_bound = scaled_lower_bound)
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
    N1 = nrow(loss_ratio),
    N2 = nrow(loss_ratio_new),
    x1 =  model.matrix(formula_x, loss_ratio),
    y1 = pmax(l_bound, c(loss_ratio[,all.vars(update(formula_x, .~0)) ][[1]])),
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
    dimSigma = 11-2 + 1 *(virtual_data == 1) 
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
    stanRun <- sampling(object = stanMod, data = stan_list, iter=iter, warmup=warmup, chains=chains,init=init_list,
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

    return(list(loss_square=loss_square_pred,loss_square_med=loss_square_pred_med, loss_square_2=loss_square_pred_2, 
                rank_val_incremental=rank_val_incremental,rank_val_cumulative= rank_val_cumulative, rank_ultimate= rank_ultimate, sigma = sigma,
                errors = errors, summary_stats = summary_stats, 
                scores = c(crps1=crps1, crps2=crps2, nlpd=nlpd)))
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
stan_incremental_loss_ratio_virtual = function(loss_triangle, premiums=NULL, loss_square = NULL, calc_ranks = FALSE, 
                                               formula_x = y_scaled ~ AccidentYear_scaled + DevelopmentLag_scaled -1,
                                               formula_H = y_scaled ~  AccidentYear_scaled + DevelopmentLag_scaled -1 , formula_cat = NULL,
                                               eps = 0.001, plot_graphs = FALSE,
                                               chains = 6,iter = 600, warmup =300, adapt_delta = 0.9,max_treedepth = 10,
                                               prior_eta = 1, prior_rho = c(5.566019, 5.570317), prior_theta = 1, prior_sigma = 1, 
                                               prior_beta = 1, take_logs = FALSE, 
                                               scaled_lower_bound = NULL, l_bound = 0, mu=9:1/9-.5,
                                               virtual_data = FALSE,virtual_point = NA, 
                                               stanFile = "gp_compound_mean_03.stan", use_vb = FALSE){
  
  if(is.null(loss_square)) 
    loss_square = matrix(NA, nrow(loss_triangle), ncol(loss_triangle))
  
  AccidentYear_unique = AccidentYear = as.numeric(rownames(loss_triangle))
  DevelopmentLag_unique = DevelopmentLag = as.numeric(colnames(loss_triangle))
  DevelopmentLag = DevelopmentLag[-1]
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
  loss_ratio_triangle = loss_ratio_triangle[,-1] - loss_ratio_triangle[,-ncol(loss_ratio_triangle)]
  loss_ratio_triangle_log = log(apply((loss_triangle / premiums),2,function(x) pmax(x,eps)))
  loss_ratio_triangle_log = loss_ratio_triangle_log[,-1] - loss_ratio_triangle_log[,-ncol(loss_ratio_triangle_log)]
  
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
    # y1 = pmax(l_bound, c(loss_ratio[,all.vars(update(formula_x, .~0)) ][[1]]))[actual_data_flag],
    # y1_virtual = pmax(l_bound, c(loss_ratio[,all.vars(update(formula_x, .~0)) ][[1]])),
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
                        control = list(adapt_delta=adapt_delta, max_treedepth =max_treedepth) )
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
    
    
    return(list(loss_square=loss_square_pred,loss_square_med=loss_square_pred_med, loss_square_2=loss_square_pred_2, 
                rank_val_incremental=rank_val_incremental,rank_val_cumulative= rank_val_cumulative, sigma = sigma,
                errors = errors, summary_stats = summary_stats))
  }
}




