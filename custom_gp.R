
####################################################3
##  Customized Gaussian Process Functions
####################################################3


#' Gaussian Process Model
#' 
#' @param x1 design matrix of inputs
#' @param y target values to be modeled
#' @param x2 (Optional) New input values over which to predict new y values.  If ignored, x2 is set to x1.
#' @param prior_eta_sq, prior_theta, prior_sigma_sq Dispersion values for the priors on eta_sq, theta, sigma-sq.   eta_sq = amplitude, theta = lengthscale, sigma_sq = noise
#' @details Fits a GP regression using Greta.  Data is scaled internally. Prior's relate to scaled data.
#' @results A fitted square exponential model
gp_greta = function(x1, y, x2 = NULL,   
                    prior_eta_sq = 1, prior_theta = 1, prior_sigma_sq =1,
                    plot_graphs = TRUE, n_samples = 100, warmup=100) {
  
  x1=as.matrix(x1)
  
  if (is.null(x2)){
    x2 = x1
  } else {
    x2= as.matrix(x2)
  }
  
  y= as.matrix(y)
  
  # Dimensionality of inputs
  D=ncol(x1) 
  N1 = nrow(x1)
  N2 = nrow(x2)
  
  # Scale Data 
  mean_x1 = mean(x1)
  sd_x1=sd(x1)
  scale_x1 = (x1-mean_x1)/sd_x1
  scale_x2 = (x2-mean_x1)/sd_x1
  
  # Calculate distances between predictors
  xd = distance_calc(scale_x1,scale_x2)
  xd1=as_data(xd$xd1)
  xd2=as_data(xd$xd2)
  xd12=as_data(xd$xd12)
  
  ####################################
  ####################################
  ###  Start of Greta Model
  y1=greta_array(y, dim=c(1,N1))
  dim(y1)
  
  eye1 = greta_array(diag(N1), dim=c(N1,N1))
  eye2 = greta_array(diag(N2), dim=c(N2,N2))
  
  #Priors
  eta_sq = lognormal(0,prior_eta_sq,dim = 1)
  sigma_sq =  normal(0,prior_sigma_sq,dim = 1, truncation = c(0, Inf))
  theta = normal(0, prior_theta,dim = D, truncation = c(0, Inf))
  rho_sq = theta ^  - 2
  mu = zeros(dim=N1)
  
  #Calculate Covariance
  dist =  greta_array(xd1[1,,], dim=c(N1,N1)) * rho_sq[1]
  if (D>1){
    for (d in 2:D){
      dist = dist +  greta_array(xd1[d,,], dim=c(N1,N1)) * rho_sq[d];
    }
  }
  K1 = eta_sq * exp(-dist) + eye1*(sigma_sq + 0.000005)
  
  # Distribution
  
  # Cholesky alternative
  if (T) {
    z <- normal(0, 1, dim = N1)
    L <- t(chol(K1))
    distribution(y) = normal(mu + L %*% z,1)
  } else {
    distribution(y) = multivariate_normal(mu, K1)
  }
  ####################################
  ####################################
  
  
  # Define Model
  m= model(eta_sq, theta, rho_sq, sigma_sq)
  
  if (plot_graphs)
    plot(m)
  
  ### Optimization
  if (F) {
    a= Sys.time()
    map_pars = opt(m)
    b= Sys.time(); b-a
    
    map_pars$par
  }
  
  ### MCMC
  a= Sys.time()
  draws=greta::mcmc(m, n_samples = n_samples, warmup=warmup)
  b= Sys.time(); b-a
  draws
  mcmc_trace(draws)
  mcmc_intervals(draws)
  
  sum_draws = summary(draws)
  sum_draws
  #gelman.diag(draws)
  n_eff = effectiveSize(draws)

  
  
  ######################################################################################
  # Make Predictions through R
  ######################################################################################
  
  
  sim_gp = function(y, eta_sq, rho_sq, sigma_sq, D){
    dist1 =  xd$xd1[1,,] * rho_sq[1] 
    dist2 =  xd$xd2[1,,] * rho_sq[1] 
    dist12 =  xd$xd12[1,,] * rho_sq[1]
    if (D>1){
      for (d in 2:D){
        dist1 = dist1 +  xd$xd1[d,,]* rho_sq[d]
        dist2 = dist2 +  xd$xd2[d,,]* rho_sq[d]
        dist12 = dist12 +  xd$xd12[d,,]* rho_sq[d]
      }
    }
    K1 = eta_sq * exp(-dist1) + diag(N1)*(sigma_sq + 0.00005)
    K2 = eta_sq * exp(-dist2) + diag(N2)*(sigma_sq + 0.00005)
    K12 = eta_sq * exp(-dist12)
    
    K12_transpose_div_K1 = t(K12) %*% solve(K1)
    mu_star = K12_transpose_div_K1 %*% y
    Tau = K2 - K12_transpose_div_K1 %*% K12
    Tau = Matrix::forceSymmetric(Tau,uplo="U") %>% as.matrix()
    y_star = c(rmvnorm(1, c(mu_star), Tau, method = "chol"))
    return(y_star)
    
  }

  no_sims = nrow(draws[[1]])
  sims = matrix(0, no_sims, N2)
  rho_sq_sims = as_tibble(draws[[1]]) %>% select(starts_with("rho_sq")) %>% as.matrix()
  
  for (n in 1:no_sims){
    print(n)
    sims[n,] = sim_gp(y = y, eta_sq = draws[[1]][n,"eta_sq"], 
                      rho_sq = rho_sq_sims[n,] , 
                      sigma_sq = draws[[1]][n,"sigma_sq"] , D = D)
  }
  sims
  res = as.data.frame(x2) %>%
    mutate(predicted = colMeans(sims), 
           lower = apply(sims, 2, function(x) quantile(x, probs = 0.05)),
           upper = apply(sims, 2, function(x) quantile(x, probs = 0.95))
    )
  res

}



#' Fit of  a Gaussian Process with a Trend Function.  Uses Fast calculation methodology.
#' 
#' Prediction can be done on out-of-sample datasets
#' 
#' @param data A data frame containing the input training set.  
#' @param data_new The out-of-sample prediction set.  If NULL, data_new is set to data.  No target values are required.
#' @param formula_x The GP regression formula
#' @param forumla_H The right hand side formula for the predictors of the mean function
#' @param prior_eta_sq Prior for eta_sq parameter
#' @param prior_theta Prior for theta parameter
#' @param prior_sigma_sq Prior for sigma_sq parameter
#' @param optimizeMod Default is \code{TRUE}.  If \TRUE}, then an optimized set of parameters are provided.  If \code{FALSE}, then MCMC is conducted
#' @param iter Number of iterations if MCMC is run
#' @param chains Number of chains if MCMC is run
#' @return The MCMC chains or the optimized valus of the parameters of the Gaussian Process model
gp_stan <- function(data, data_new = NULL, formula_x, formula_H = NULL, 
                    prior_eta_sq = 1, prior_theta = 1, prior_sigma_sq = 1,  
                    adapt_delta=0.99, max_treedepth = 13,
                    optimizeMod = TRUE, iter =400, warmup = iter/4, chains = 4,
                    stanFile = "gpFastPredict04.stan"){
  


  if (is.null(data_new)){
    data_new = data
  }
  
  # Establish data to be passed to Stan
  N1 = nrow(data)
  N2 = nrow(data_new)
  
  data_new$blank = 0
  
  stan_list <- list(
    N1 = N1,
    N2 = N2,
    x1 = model.matrix(formula_x, data),
    H1 = model.matrix(formula_H, data), 
    y1 = c(data[,all.vars(update(formula_x, .~0)) ][[1]]),
    x2 = model.matrix(update(formula_x, blank ~ .), data_new),
    H2 = model.matrix(formula_H, data_new), 
    prior_eta_sq=prior_eta_sq,
    prior_theta=prior_theta,
    prior_sigma_sq=prior_sigma_sq
  )
  stan_list$D = ncol(stan_list$x1)
  stan_list$dimH <- ncol(stan_list$H1)
  
  if (ncol(stan_list$H1) != ncol(stan_list$H2)){
    stop("data and data_new have different dimensions for mean function ")
  }
  
  stanMod <- stan_model(file=paste0(stan_dir, stanFile))

  #save(stanMod, file=paste0(stanDirectory, "gpFast01.rda"))
  #load(paste0(stanDirectory, "gpFast01.rda"))
  
  if (!optimizeMod){
    stanRun <- sampling(object = stanMod, data = stan_list, iter=iter, chains=chains, control = list(adapt_delta=adapt_delta, max_treedepth =max_treedepth) )
  } else {
    stanRun <- optimizing(object = stanMod, data = stan_list, verbose=T)
  }
  stanRun
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
gp_stan_cat <- function(data, data_new = NULL, formula_x, formula_H = NULL, formula_cat = NULL, 
                    prior_eta_sq = 1, prior_theta = 1, prior_alpha = 1, prior_sigma_sq = 1, prior_beta = 1,
                    priorDist_theta="lognormal", prior_nu = 4, 
                    adapt_delta=0.99, max_treedepth = 13,
                    optimizeMod = TRUE, iter =400, chains = 4, warmup = 100,
                    stanFile = "gpFastPredict_categorical.stan"){
  
  
  
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
    H1 = model.matrix(formula_H, data), 
    y1 = c(data[,all.vars(update(formula_x, .~0)) ][[1]]),
    x2 = model.matrix(upgp_date(formula_x, blank ~ .), data_new),
    H2 = model.matrix(formula_H, data_new), 
    prior_eta_sq=prior_eta_sq,
    prior_theta=prior_theta,
    prior_beta=prior_beta,
    prior_alpha=prior_alpha,
    prior_sigma_sq=prior_sigma_sq
  )
  
  if(!is.null(formula_cat)){
    formula_cat = update(formula_cat,  ~ . -1) # remove intercept term
    stan_list$x_cat1 = c(model.matrix(formula_cat, data))
    stan_list$x_cat2 = c(model.matrix(formula_cat, data_new))
  }
    
  stan_list$D = ncol(stan_list$x1)
  stan_list$dimH <- ncol(stan_list$H1)
  
  if (ncol(stan_list$H1) != ncol(stan_list$H2)){
    stop("data and data_new have different dimensions for mean function ")
  }
  
  stanMod <- stan_model(file=paste0(stan_dir, stanFile))
  
  #save(stanMod, file=paste0(stanDirectory, "gpFast01.rda"))
  #load(paste0(stanDirectory, "gpFast01.rda"))
  
  if (!optimizeMod){
    stanRun <- sampling(object = stanMod, data = stan_list, iter=iter, warmup=warmup, chains=chains, control = list(adapt_delta=adapt_delta, max_treedepth =max_treedepth) )
  } else {
    stanRun <- optimizing(object = stanMod, data = stan_list, verbose=T)
  }
  stanRun
  plot(stanRun, pars = "yStar")
  
  if(F){
    preds = colMeans(rstan::extract(stanRun, pars = "yStar")[[1]])
    str(preds)
    data_new$preds =  preds
    head(data_new)
    
    ##### Loss Reserve Check
    ggplot(data_new) + 
      geom_line(aes(DevelopmentYear, preds, color = as.factor(AccidentYear))) + 
      facet_wrap(~ as.factor(GRCODE), scales = "free")
    ##### Interest Rate Check
    data_new$preds = scaler_yield(preds, unscale = T)
    ggplot(data_new) + 
      geom_line(aes(term, preds, color = "GP Predicted")) + 
      geom_point(aes(term, yield, color = "Actual")) + 
      facet_wrap(~as.factor(date)) 
    # + 
    #   facet_wrap(~ as.factor(GRCODE), scales = "free")
    
  }
  
  ans = list(draws = stanRun, mu_pred = rstan::extract(stanRun, pars = "yStar")[[1]], tau_pred = NA)
  ans
}


#' Fit of  a Gaussian Process with a Trend Function Using the Kronecker Method
#' 
#' Prediction can be done on out-of-sample datasets
#' 
#' @param data A data frame containing the input training set.  
#' @param data_new The out-of-sample prediction set.  If NULL, data_new is set to data.  No target values are required.
#' @param name_x_1, name_x_2  The column names of the first and second variables of data and data_new using in GP regression
#' @param forumla_H The right hand side formula for the predictors of the mean function
#' @param forumla_cat The right hand side formula for the single categorical predictor.  The designated column must contain real numbers.
#' @param prior_eta_sq Prior for eta_sq parameter
#' @param prior_theta Prior for theta parameter
#' @param prior_sigma_sq Prior for sigma_sq parameter
#' @param optimizeMod Default is \code{TRUE}.  If \TRUE}, then an optimized set of parameters are provided.  If \code{FALSE}, then MCMC is conducted
#' @param iter Number of iterations if MCMC is run
#' @param chains Number of chains if MCMC is run
#' @return The MCMC chains or the optimized valus of the parameters of the Gaussian Process model
gp_stan_kronecker <- function(data, data_new = NULL, name_x_1, name_x_2, 
                              target, formula_H = NULL, formula_cat = NULL, 
                              prior_eta_sq = 1, prior_theta = 1, prior_alpha = 1, prior_sigma_sq = 1, 
                              adapt_delta=0.9, max_treedepth = 10,
                              optimizeMod = TRUE, iter =400, chains = 6, warmup=100,
                              stanFile = "gp_kronecker_01.stan"){

  if (is.null(data_new)){
    data_new = data
  }
  data=as_tibble(data)
  data_new= as_tibble(data_new)
  
  var_1 = sym(name_x_1)
  var_2 = sym(name_x_2)
  
  data = data %>%
    arrange(!! var_1, !! var_2)
  data_new = data_new %>%
    arrange(!! var_1, !! var_2)
  
  #check size
  N_train_1 = length(unique(pull(data, name_x_1)))
  N_train_2 = length(unique(pull(data, name_x_2)))
  N_test_1 = length(unique(pull(data_new, name_x_1)))
  N_test_2 = length(unique(pull(data_new, name_x_2)))
  
  if (nrow(data) != N_train_1 * N_train_2){
    stop("A complete grid of values must be provided for the training set")
  }
  if (nrow(data) != N_test_1 * N_test_2){
    stop("A complete grid of values must be provided for the test set")
  }
  
  # Establish data to be passed to Stan
  
  stan_list <- list(
    N_train_1 = N_train_1,
    N_train_2 = N_train_2,
    N_test_1 = N_test_1,
    N_test_2 = N_test_2,
    x_train_1 = unique(pull(data, name_x_1)),
    x_train_2 = unique(pull(data, name_x_2)),
    x_test_1 = unique(pull(data_new, name_x_1)),
    x_test_2 = unique(pull(data_new, name_x_2)),
    H1 = model.matrix(formula_H, data), 
    H2 = model.matrix(formula_H, data_new), 
    prior_eta_sq=prior_eta_sq,
    prior_theta=prior_theta,
    prior_alpha=prior_alpha,
    prior_sigma_sq=prior_sigma_sq
  )
  stan_list$y = matrix(pull(data, target),nrow = stan_list$N_train_2, ncol = stan_list$N_train_1)
  
  
  # if(!is.null(formula_cat)){
  #   formula_cat = update(formula_cat,  ~ . -1) # remove intercept term
  #   stan_list$x_cat1 = c(model.matrix(formula_cat, data))
  #   stan_list$x_cat2 = c(model.matrix(formula_cat, data_new))
  # }
  
  stan_list$dimH <- ncol(stan_list$H1)
  
  if (ncol(stan_list$H1) != ncol(stan_list$H2)){
    stop("data and data_new have different dimensions for mean function ")
  }
  
  stanMod <- stan_model(file=paste0(stan_dir, stanFile))
  
  #save(stanMod, file=paste0(stanDirectory, "gpFast01.rda"))
  #load(paste0(stanDirectory, "gpFast01.rda"))
  
  if (!optimizeMod){
    stanRun <- sampling(object = stanMod, data = stan_list, iter=iter, warmup=warmup, chains=chains, control = list(adapt_delta=adapt_delta, max_treedepth =max_treedepth) )
  } else {
    stanRun <- optimizing(object = stanMod, data = stan_list, verbose=T)
  }
  stanRun
  #plot(stanRun, pars = "yStar")

  if(F){
    preds = colMeans(rstan::extract(stanRun, pars = "yStar")[[1]])
    str(preds)
    data_new$preds =  preds
    head(data_new)
    
    ##### Loss Reserve Check
    ggplot(data_new) + 
      geom_line(aes(DevelopmentYear, preds, color = as.factor(AccidentYear))) + 
      facet_wrap(~ as.factor(GRCODE), scales = "free")
    ##### Interest Rate Check
    data_new$preds = scaler_yield(preds, unscale = T)
    ggplot(data_new) + 
      geom_line(aes(term, preds, color = "GP Predicted")) + 
      geom_point(aes(term, yield, color = "Actual")) + 
      facet_wrap(~as.factor(date)) 
    # + 
    #   facet_wrap(~ as.factor(GRCODE), scales = "free")
    
  }
  
  #ans = list(draws = stanRun, mu_pred = rstan::extract(stanRun, pars = "yStar")[[1]], tau_pred = NA)
  ans = list(draws = stanRun, mu_pred = NA, tau_pred = NA)
  ans
}





#' Fit of  a Gaussian Process with a Trend Function and a single categorical variable using Greta.
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
gp_greta_cat <- function(data, data_new = NULL, formula_x, formula_H = NULL, formula_cat = NULL,
                        prior_eta_sq = 1, prior_theta = 1, prior_alpha = 1, prior_beta = 1, 
                        prior_sigma_sq = 1, priorDist_theta="lognormal",
                        optimizeMod = TRUE, n_samples= 400, warmup = 50){
  
  
  
  if (is.null(data_new)){
    data_new = data
  }
  
  if(F){
    ggplot(data) + geom_line(aes(DevelopmentLag, ldf_scaled,color = as.factor(AccidentYear))) + 
      facet_wrap(~ as.factor(GRCODE), scales = "free")
    data %>% filter(GRCODE == 13528)
  }
  
  # Establish data to be passed to Stan
  N1 = nrow(data)
  N2 = nrow(data_new)
  
  data_new$blank = 0
  #data$blank = 0
  
  formula_x = update(formula_x, . ~ . -1) # remove intercept term
  formula_cat = update(formula_cat,  ~ . -1) # remove intercept term
  
  stan_list <- list(
    N1 = N1,
    N2 = N2,
    x1 = model.matrix(formula_x, data),
    H1 = model.matrix(formula_H, data), 
    x_cat1 = c(model.matrix(formula_cat, data)), 
    y1 = c(data[,all.vars(update(formula_x, .~0)) ][[1]]),
    x2 = model.matrix(update(formula_x, blank ~ .), data_new),
    H2 = model.matrix(formula_H, data_new), 
    x_cat2 = c(model.matrix(formula_cat, data_new)), 
    prior_eta_sq=prior_eta_sq,
    prior_theta=prior_theta,
    prior_alpha=prior_alpha,
    prior_sigma_sq=prior_sigma_sq
  )
  stan_list$D = ncol(stan_list$x1)
  stan_list$dimH <- ncol(stan_list$H1)
  
  if (ncol(stan_list$H1) != ncol(stan_list$H2)){
    stop("data and data_new have different dimensions for mean function ")
  }
  
  
  ######################################
  ### Greta Model ######################
  ######################################
  
  # Calculate distances between predictors
  print("Calculating distances for covariance matrix...")
  xd = distance_calc(stan_list$x1,stan_list$x2)
  xd1 = xd$xd1
  xd2 = xd$xd2
  xd12 = xd$xd12
  rm(xd)
  
  
  #xd1=as_data(dist_calc(stan_list$x1))
  #xd2=as_data(dist_calc(stan_list$x2))
  #xd12=as_data(dist_calc(stan_list$x1, stan_list$x2))

  # Calculate categorical variable matching matrices
  cat_mat1 = as_data(cat_calc(stan_list$x_cat1))
  cat_mat2 = as_data(cat_calc(stan_list$x_cat2))
  cat_mat12 = as_data(cat_calc(stan_list$x_cat1, stan_list$x_cat2))
  
  
  
  #/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/
  ###  Start of Greta Model
  print("Running Greta model")
  eye1 = greta_array(diag(N1), dim=c(N1,N1))
  eye2 = greta_array(diag(N2), dim=c(N2,N2))
  
  H1 = as_data(stan_list$H1)
  H2 = as_data(stan_list$H2)
  
  y=greta_array(stan_list$y1, dim=c(1,N1))
  dim(y)

  #Priors
  eta_sq = lognormal(0,prior_eta_sq,dim = 1)
  alpha = lognormal(0,prior_alpha,dim = 1)
  sigma_sq =  normal(0,prior_sigma_sq,dim = 1, truncation = c(0, Inf))
  theta = normal(0, prior_theta,dim = stan_list$D, truncation = c(0, Inf))
  beta_par =  normal(0,prior_beta,dim = stan_list$dimH)
  
  rho_sq = theta ^  - 2
  
  mu = t(H1 %*% beta_par )
  

  #Calculate Covariance
  dist =  greta_array(xd1[1,,], dim=c(N1,N1)) * rho_sq[1] #xd1[1,,] * rho_sq[1] 
  if (stan_list$D>1){
    for (d in 2:stan_list$D){
      dist = dist + greta_array(xd1[d,,], dim=c(N1,N1)) * rho_sq[d];
    }
  }
  
  dist = dist + alpha * cat_mat1
  
  K1 = eta_sq * exp(-dist) + eye1*(sigma_sq + 0.000005)
  
  # Distribution
  distribution(y) = multivariate_normal(mu, K1)
  

  # Define Model
  m= model(eta_sq, theta, sigma_sq, alpha, beta_par)
  
  #plot(m)
  
  ### Optimization
  if (F) {
    a= Sys.time()
    map_pars = opt(m)
    b= Sys.time(); b-a
    
    map_pars$par
  }
  
  

  ### MCMC
  #init = m$dag$example_parameters()
  #init = rep(0.1, length(init))
  a= Sys.time()
  #draws=mcmc(m, n_samples = n_samples, warmup=warmup, initial_values = init)
  draws=greta::mcmc(m, n_samples = n_samples, warmup=warmup, one_by_one = TRUE)
  b= Sys.time(); 
  print(paste("Time taken for MCMC:",b-a))

  
  
  ###########################################3
  # Make Predictions through Greta
  print("Making projections")
  
  dist2 =  greta_array(xd2[1,,] , dim=c(N2,N2))  * greta_array(rho_sq[1], dim=c(1))
  dist12 =  greta_array(xd12[1,,], dim=c(N1,N2)) * greta_array(rho_sq[1], dim=c(1))
  if (stan_list$D>1){
    for (d in 2:stan_list$D){
      dist2 = dist2 +  greta_array(xd2[d,,], dim=c(N2,N2)) *  greta_array(rho_sq[d], dim=c(1))
      dist12 = dist12 +  greta_array(xd12[d,,], dim=c(N1,N2)) *  greta_array(rho_sq[d], dim=c(1))
    }
  }
  dist2 = dist2 + alpha * cat_mat2
  dist12 = dist12 + alpha * cat_mat12
  
  # K2 = greta_array(eta_sq * exp(dist2) + diag(N2)*(sigma_sq + 0.000005),dim=c(N2,N2))
  # K12 = greta_array(eta_sq * exp(dist12), dim=c(N1,N2))
  K2 = eta_sq * exp(-dist2) + diag(N2)*(sigma_sq + 0.000005)
  K12 = eta_sq * exp(-dist12)
  
  div_K1 = solve(K1)
  K12_transpose_div_K1 = greta_array(t(K12) %*% div_K1, dim=c(N2, N1))
  #K12_transpose_div_K1 = t(K12) %*% solve(K1)
  
  mu_star = greta_array(H2 %*% beta_par + K12_transpose_div_K1 %*% (t(y) - H1 %*% beta_par) , dim = c(N2,1))
  
  mu_pred = calculate(mu_star, draws)
  
  #Uncomment if you want to sim Tau
  # Tau = greta_array(K2 - K12_transpose_div_K1 %*% K12, dim=c(N2,N2)) +
  #   (H2 - K12_transpose_div_K1 %*% H1) %*% solve(t(H1) %*% div_K1 %*% H1) %*% t(H2 - K12_transpose_div_K1 %*% H1)
  # #Tau_pred = calculate(Tau, draws)
  # y_star = multivariate_normal(t(mu_star), Tau)
  # calculate(y_star)
  
  Tau_pred = list(NA)
  
  ans = list(draws = draws, mu_pred = mu_pred[[1]], tau_pred = Tau_pred[[1]])
  ans
  
}



#' Fit of  a Gaussian Process with a Trend Function  using Greta. Full Version Using Transformed MultivariateNormal  
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
gp_greta_full <- function(data, data_new = NULL, formula_x, formula_H = NULL, max_iterations=100,
                         prior_eta_sq = 1, prior_theta = 1, prior_beta = 1, 
                         prior_sigma_sq = 1, 
                         optimizeMod = TRUE, n_samples= 400, warmup = 50){
  
  
  
  if (is.null(data_new)){
    data_new = data
  }
  data=as_tibble(data)
  data_new=as_tibble(data_new)
  
  data_new$blank = 0
  
  # Establish data to be passed to Stan
  N1 = nrow(data)
  N2 = nrow(data_new)
  
  formula_x = update(formula_x, . ~ . -1) # remove intercept term
  
  stan_list <- list(
    N1 = N1,
    N2 = N2,
    x1 = model.matrix(formula_x, data),
    H1 = model.matrix(formula_H, data), 
    y1 = matrix(pull(data,all.vars(update(formula_x, .~0))),1),
    x2 = model.matrix(update(formula_x, blank ~ .), data_new),
    H2 = model.matrix(formula_H, data_new), 
    prior_eta_sq=prior_eta_sq,
    prior_theta=prior_theta,
    prior_sigma_sq=prior_sigma_sq,
    prior_beta = prior_beta
  )
  stan_list$D = ncol(stan_list$x1)
  stan_list$dimH <- ncol(stan_list$H1)
  
  
  if (ncol(stan_list$H1) != ncol(stan_list$H2)){
    stop("data and data_new have different dimensions for mean function ")
  }
  
  
  ######################################
  ### Greta Model ######################
  ######################################
  
  # Calculate distances between predictors
  print("Calculating distances for covariance matrix...")
  xd = distance_calc(stan_list$x1,stan_list$x2)
  xd1 = as_data(xd$xd1)
  xd2 = as_data(xd$xd2)
  xd12 = as_data(xd$xd12)
  rm(xd)
  
  
  #/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/
  ###  Start of Greta Model
  print("Running Greta model")
  eye1 = greta_array(diag(N1), dim=c(N1,N1))
  eye2 = greta_array(diag(N2), dim=c(N2,N2))
  
  H1 = as_data(stan_list$H1)
  H2 = as_data(stan_list$H2)
  
  y=greta_array(stan_list$y1, dim=c(N1,1))
  dim(y)
  
  #Priors
  eta_sq = lognormal(0,prior_eta_sq,dim = 1)
  sigma_sq =  normal(0,prior_sigma_sq,dim = 1, truncation = c(0, Inf))
  theta = lognormal(0, prior_theta,dim = stan_list$D, truncation = c(0, Inf))
  beta_par =  normal(0,prior_beta,dim = stan_list$dimH)
  
  #rho_sq = theta ^  - 2
  
  mu = H1 %*% beta_par 
  
  
  #Calculate Covariance
  dist =  greta_array(xd1[1,,], dim=c(N1,N1)) / theta[1]^2 #xd1[1,,] * rho_sq[1] 
  if (stan_list$D>1){
    for (d in 2:stan_list$D){
      dist = dist + greta_array(xd1[d,,], dim=c(N1,N1)) / theta[d] ^ 2    
    }
  }
  
  tol = 0.00005
  K1 = eta_sq * exp(dist) + eye1*(tol)
  
  # Distribution
  
  if (T) {
    z <- normal(0, 1, dim = N1)
    L <- t(chol(K1))
    distribution(y) = normal(mu + L %*% z, sigma_sq^.5 ) 
  } else {
    distribution(y) = multivariate_normal(mu, K1 + eye1*(sigma_sq))
  }
  
  
  
  
  # Define Model
  m= model(eta_sq, theta, sigma_sq, beta_par)
  
  #plot(m)
  
  ### Optimization
  if (optimizeMod) {
    a= Sys.time()
    map_pars = opt(m, control=list(learning_rate = 0.8), max_iterations = max_iterations)
    b= Sys.time(); b-a
    
    map_pars$par
    return(map_pars)
  }
  
  
  
  ### MCMC
  #init = m$dag$example_parameters()
  #init = rep(0.1, length(init))
  a= Sys.time()
  #draws=mcmc(m, n_samples = n_samples, warmup=warmup, initial_values = init)
  draws=greta::mcmc(m, n_samples = n_samples, warmup=warmup)
  b= Sys.time(); 
  print(paste("Time taken for MCMC:",b-a))
  
  
  
  ###########################################3
  # Make Predictions through Greta
  print("Making projections")
  
  dist2 =  greta_array(xd2[1,,] , dim=c(N2,N2))  / greta_array(theta[1]^2, dim=c(1))
  dist12 =  greta_array(xd12[1,,], dim=c(N1,N2)) / greta_array(theta[1]^2, dim=c(1))
  if (stan_list$D>1){
    for (d in 2:stan_list$D){
      dist2 = dist2 +  greta_array(xd2[d,,], dim=c(N2,N2)) /  greta_array(theta[d]^2, dim=c(1))
      dist12 = dist12 +  greta_array(xd12[d,,], dim=c(N1,N2)) /  greta_array(theta[d]^2, dim=c(1))
    }
  }
  
  # K2 = greta_array(eta_sq * exp(dist2) + diag(N2)*(sigma_sq + 0.000005),dim=c(N2,N2))
  # K12 = greta_array(eta_sq * exp(dist12), dim=c(N1,N2))
  K2 = eta_sq * exp(dist2) + diag(N2)*(sigma_sq + 0.000005)
  K12 = eta_sq * exp(dist12)
  
  div_K1 = solve(K1 + eye1*(sigma_sq ))
  K12_transpose_div_K1 = greta_array(t(K12) %*% div_K1, dim=c(N2, N1))
  #K12_transpose_div_K1 = t(K12) %*% solve(K1)
  
  #mu_star = greta_array(H2 %*% beta_par + K12_transpose_div_K1 %*% (t(y) - mu) , dim = c(N2,1))
  mu_star = greta_array(H2 %*% beta_par + K12_transpose_div_K1 %*% (y - mu) , dim = c(N2,1))
  
  mu_pred = calculate(mu_star, draws)
  
  if (F){
    opt_list = as.list(map_pars$par)
    mu_pred = calculate(mu_star, list(eta_sq = 1, theta= c(1,1), sigma_sq=0.94, beta_par=-1.12))
    
  }
  
  #Uncomment if you want to sim Tau
  # Tau = greta_array(K2 - K12_transpose_div_K1 %*% K12, dim=c(N2,N2)) +  
  #   (H2 - K12_transpose_div_K1 %*% H1) %*% solve(t(H1) %*% div_K1 %*% H1) %*% t(H2 - K12_transpose_div_K1 %*% H1)
  #Tau_pred = calculate(Tau, draws)
  Tau_pred = list(NA)
  
  ans = list(draws = draws, mu_pred = mu_pred[[1]], tau_pred = Tau_pred[[1]])
  ans
  
}




#' Fit of  a Gaussian Process with a Trend Function  using Greta. Full Version No Transformation  
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
gp_greta_notrans <- function(data, data_new = NULL, formula_x, formula_H = NULL, 
                          prior_eta_sq = 1, prior_theta = 1, prior_beta = 1, 
                          prior_sigma_sq = 1, 
                          optimizeMod = TRUE, max_iterations=100,learning_rate = 0.8, 
                          n_samples= 400, warmup = 50, thin = 1){
  
  
  
  if (is.null(data_new)){
    data_new = data
  }
  data=as_tibble(data)
  data_new=as_tibble(data_new)
  
  data_new$blank = 0
  
  # Establish data to be passed to Stan
  N_train = nrow(data)
  N_test = nrow(data_new)
  
  formula_x = update(formula_x, . ~ . -1) # remove intercept term
  
  x_train = model.matrix(formula_x, data)
  x_test = model.matrix(update(formula_x, blank ~ .), data_new)
  y = greta_array(pull(data,all.vars(update(formula_x, .~0))),dim = c(1, N_train)) 
  
  H1 = as_data(model.matrix(formula_H, data) )
  H2 = as_data(model.matrix(formula_H, data_new))

  D = ncol(x_train)
  dimH <- ncol(H1)
  
  # Calculate distances between predictors
  print("Calculating distances for covariance matrix...")
  xd = distance_calc(x_train,x_test)
  xd1 = xd$xd1
  xd2 = xd$xd2
  xd12 = xd$xd12
  rm(xd)
  
  ######################################
  ### Greta Model ######################
  ######################################
  
  
  
  #/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/
  ###  Start of Greta Model
  print("Running Greta model")
  eye1 = greta_array(diag(N_train), dim=c(N_train,N_train))
  eye2 = greta_array(diag(N_test), dim=c(N_test,N_test))
  

  #Priors
  eta_sq = greta::lognormal(0,prior_eta_sq,dim = 1)
  sigma_sq =  greta::normal(0,prior_sigma_sq,dim = 1, truncation = c(0, Inf))
  theta = greta::lognormal(matrix(0, length(prior_theta)), prior_theta,dim = D, truncation = c(0, Inf))
  beta_par =  greta::normal(0,prior_beta,dim = dimH)
  # eta_sq=(c(1))
  # #sigma_sq=(c(.1))
  # theta=(c(1,1))
  # #beta_par =  c(0)
  
  mu = H1 %*% beta_par 
  
  
  #Calculate Covariance
  dist =  greta_array((xd1[1,,]), dim = c(N_train, N_train)) / theta[1]^2 #xd1[1,,] * rho_sq[1] 
  if (D>1){
    for (d in 2:D){
      dist = dist + greta_array((xd1[d,,]), dim = c(N_train, N_train)) / theta[d] ^ 2    
    }
  }
  
  tol = as_data(0.0005)
  
  # Distribution
  TRANSFORM = FALSE
  if (TRANSFORM) {
    K1 = eta_sq * exp(dist) + diag(tol, N_train, N_train)
    z <- normal(0, 1, dim = N_train)
    L <- t(chol(K1))
    distribution(y) = normal(mu + L %*% z, sigma_sq^.5 ) 
  } else {
    K1 = eta_sq * exp(dist) + (sigma_sq + tol) * eye1
    distribution(y) = multivariate_normal(mu, K1)
  }
  
  # Define Model
  m= model(eta_sq, theta, sigma_sq, beta_par)
  #m= model(beta_par)
  
  #plot(m)
  
  ### Optimization
  if (optimizeMod) {
    a= Sys.time()
    map_pars = opt(m, control=list(learning_rate = learning_rate), max_iterations = max_iterations)
    b= Sys.time(); b-a
    
    map_pars$par
    return(map_pars)
  }
  
  
  
  ### MCMC
  #init = m$dag$example_parameters()
  #init = rep(0.1, length(init))
  a= Sys.time()
  #draws=mcmc(m, n_samples = n_samples, warmup=warmup, initial_values = init)
  draws=greta::mcmc(m, n_samples = n_samples, warmup=warmup, thin = thin)
  b= Sys.time(); 
  print(paste("Time taken for MCMC:",b-a))
  
  
  
  ###########################################3
  # Make Predictions through Greta
  print("Making projections")
  
  dist2 =  greta_array(xd2[1,,] , dim=c(N_test,N_test))  / greta_array(theta[1]^2, dim=c(1))
  dist12 =  greta_array(xd12[1,,], dim=c(N_train,N_test)) / greta_array(theta[1]^2, dim=c(1))
  if (D>1){
    for (d in 2:D){
      dist2 = dist2 +  greta_array(xd2[d,,], dim=c(N_test,N_test)) /  greta_array(theta[d]^2, dim=c(1))
      dist12 = dist12 +  greta_array(xd12[d,,], dim=c(N_train,N_test)) /  greta_array(theta[d]^2, dim=c(1))
    }
  }
  
  # K2 = greta_array(eta_sq * exp(dist2) + diag(N_test)*(sigma_sq + 0.000005),dim=c(N_test,N_test))
  # K12 = greta_array(eta_sq * exp(dist12), dim=c(N_train,N_test))
  K12 = eta_sq * exp(dist12)
  K2 = eta_sq * exp(dist2) + (sigma_sq + tol) * eye2
  
  if (TRANSFORM){
    div_K1 = solve(K1 + diag(sigma_sq, N_train, N_train))
    K12_transpose_div_K1 = greta_array(t(K12) %*% div_K1, dim=c(N_test, N_train))
    mu_star = greta_array(H2 %*% beta_par + K12_transpose_div_K1 %*% (y - mu) , dim = c(N_test,1))
  } else {
    div_K1 = solve(K1)
    K12_transpose_div_K1 = greta_array(t(K12) %*% div_K1, dim=c(N_test, N_train))
    mu_star = greta_array(H2 %*% beta_par + K12_transpose_div_K1 %*% (t(y) - mu) , dim = c(N_test,1))
  }

  mu_pred = calculate(mu_star, draws)
  
  if (F){
    opt_list = as.list(map_pars$par)
    mu_pred = calculate(mu_star, list(eta_sq = 1, theta= c(1,1), sigma_sq=0.94, beta_par=-1.12))
    
  }
  
  #Uncomment if you want to sim Tau
  # Tau = greta_array(K2 - K12_transpose_div_K1 %*% K12, dim=c(N_test,N_test)) +  
  #   (H2 - K12_transpose_div_K1 %*% H1) %*% solve(t(H1) %*% div_K1 %*% H1) %*% t(H2 - K12_transpose_div_K1 %*% H1)
  #Tau_pred = calculate(Tau, draws)
  Tau_pred = list(NA)
  
  ans = list(draws = draws, mu_pred = mu_pred[[1]], tau_pred = Tau_pred[[1]])
  ans
  
}




#' Fit of  a Gaussian Process with  Square Exponential Kernel
#' 
#' Prediction can be done on out-of-sample datasets
#' 
#' @param data A data frame containing the input training set.  
#' @param data_new The out-of-sample prediction set.  If NULL, data_new is set to data.  No target values are required.
#' @param formula_x The GP regression formula
# @param forumla_H The right hand side formula for the predictors of the mean function 
# @param forumla_cat The right hand side formula for the single categorical predictor.  The designated column must contain real numbers.
#' @param prior_eta_sq Prior for eta_sq parameter
#' @param prior_rho Prior for theta parameter
#' @param prior_sigma_sq Prior for sigma_sq parameter
#' @param optimizeMod Default is \code{TRUE}.  If \TRUE}, then an optimized set of parameters are provided.  If \code{FALSE}, then MCMC is conducted
#' @param iter Number of iterations if MCMC is run
#' @param chains Number of chains if MCMC is run
#' @return The MCMC chains or the optimized valus of the parameters of the Gaussian Process model
gp_stan_se <- function(data, data_new = NULL, formula_x, formula_H = NULL, include_dev_lag = FALSE, include_premiums = FALSE,
                       formula_cat = NULL, 
                       prior_eta_sq = 1, prior_eta = 1, prior_rho = c(5,5), prior_theta = 1, 
                       prior_alpha = 1, prior_sigma = 1, prior_lambda = 1, prior_beta = 1,
                        adapt_delta=0.99, max_treedepth = 13,
                        optimizeMod = TRUE, iter =400, chains = 4, warmup = 100,
                        stanFile = "gp_se_03.stan"){
  

  if (is.null(data_new)){
    data_new = data
  }

  data=as_tibble(data)
  data_new= as_tibble(data_new)
  
  # Establish data to be passed to Stan
  data_new$blank = 0
  formula_x = update(formula_x, . ~ . -1) # remove intercept term
  
  stan_list <- list(
    N1 = nrow(data),
    N2 = nrow(data_new),
    x1 = model.matrix(formula_x, data),
    y1 = c(data[,all.vars(update(formula_x, .~0)) ][[1]]),
    x2 = model.matrix(update(formula_x, blank ~ .), data_new),
    prior_eta_sq=prior_eta_sq,
    prior_eta=prior_eta,
    prior_rho=prior_rho,
    prior_theta = prior_theta,
    prior_sigma=prior_sigma,
    prior_lambda = prior_lambda
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
  
  stan_list$D = ncol(stan_list$x1)
  
  stanMod <- stan_model(file=paste0(stan_dir, stanFile))
  
  #save(stanMod, file=paste0(stanDirectory, "gpFast01.rda"))
  #load(paste0(stanDirectory, "gpFast01.rda"))
  
  if (!optimizeMod){
    stanRun <- sampling(object = stanMod, data = stan_list, iter=iter, warmup=warmup, chains=chains, 
                        control = list(adapt_delta=adapt_delta, max_treedepth =max_treedepth) )
  } else {
    stanRun <- optimizing(object = stanMod, data = stan_list, verbose=T)
  }

  ans = list(draws=stanRun)
}


if(F){
  
  # Get data for a company
  
  data = load_losses_subset()
  loss = loss_table(data, "comauto", company_code = 353)
  loss$ldf = dev_factors(loss$paid_tr) # Get loss development factors
  loss = losses_long(loss) # make into a "long" format table
  loss$ldf = loss$ldf %>% rename(ldf=loss) #rename column heading
  premiums = get_premium(data, "comauto", company_code = 353)
  
  # Scale data
  scaler_accident_year = scaler((loss$paid_sq$AccidentYear))
  scaler_development_lag = scaler(loss$paid_sq$DevelopmentLag)
  scaler_ldf = scaler(loss$ldf$ldf)
  
  scaler_all = function(m, unscale = FALSE){
    browser()
    m = m %>%
      mutate(AccidentYear = scaler_accident_year(AccidentYear, unscale=unscale),
             DevelopmentLag = scaler_development_lag(DevelopmentLag, unscale=unscale))
    if ("ldf" %in% colnames(m)){
      m = m %>%
        mutate(ldf = scaler_ldf(ldf, unscale=unscale))
    }
    m
  }
  
  loss_scaled = map(loss, scaler_all)
  
  
  
}
