
#' Create a function to scale and unscale data
#' 
#' @param x_reference A reference matrix
#' @param scaled_sd The  targeted size of the standard deviation of the scaled data.  If the sd is 0, then the sd is set to 1 to ensure that no NaN are created
#' @param l_bound Shifts the scaled values so that the lower bound is l_bound.  If l_bound is NULL, this is ignored.
#' @return a function that calculates the scaled values on each column of matrix x if unscale = FALSE, or unscales each column of matrix x, if unscale = TRUE 
scaler = function(x_reference, scaled_sd = 1, l_bound = NULL){
  
  x_reference = as.matrix(x_reference)
  
  means = apply(x_reference,2,mean)
  sds = apply(x_reference,2,sd) / scaled_sd
  sds = ifelse(sds == 0,1,sds)
  if (!is.null(l_bound)){
    min_vals = (apply(x_reference,2,min) - means)/sds
  } else {
    min_vals = 0
    l_bound = 0
  }
  rm(x_reference)
  
  FUN = function(x, unscale = FALSE){
    x = as.matrix(x)
    if (ncol(x) != length(means)){
      stop("scaler and x size mismatch")
    }
    
    if (!unscale){
      return(t((t(x) - means)/sds - min_vals + l_bound) )
    } else {
      return(t((t(x)+min_vals-l_bound)*sds + means))
    }
  }
  
  FUN
}

# Remove NA data from a Vector
remove_na = function(x){
  x[!is.na[x]]
}

#' Calculates the shape and rate parameters for the Inverse Gamma to produce the desired tails
#' 
#' @param tail_size the size of the tail on each of the leftand right hand size of the distribution
#' @param lower_ls, upper_ls The lower and upper lengthscales
#' @param guess initial guess value for the parameters
rho_prior_tune = function(tail_size = 0.01, lower_ls, upper_ls, guess = c(10,10)){
  
  error_rate = function(pars, tail_size = 0.01, lower_ls, upper_ls){
    shape = pars[1]
    rate = pars[2]
    delta_1 =pinvgamma(q = lower_ls, shape = shape, rate) - tail_size
    delta_2 =1-pinvgamma(q = upper_ls, shape = shape, rate) - tail_size
    ans = delta_1^2 + delta_2^2
    ans
  }
  res = optim(par = guess, fn = error_rate, tail_size = tail_size, lower_ls = lower_ls, upper_ls = upper_ls)
  res
}
# rho_prior_tune(.01, 0.44, 3.55)
# rho_prior_tune(.001, 0.44, 3.55)
# rho_prior_tune(.005, 0.44, 2.5)$par

#' One-hot envode a vector
#' 
#' @param var Vector to be one-hot-encoded.The vector must be a factor
one_hot = function(var){
  N= length(var)
  vals = levels(var)
  N_vals  = length(vals)
  m= matrix(NA, N, N_vals)
  for (i in 1:N_vals){
    m[,i] = ifelse(var == vals[i],1,0)
  }
  m
}


# Creates a character index based on a given set of vectors
make_index = function( ...){
  nameVecs <- list(...) %>%
    map(as.character)
  ans = nameVecs[[1]]
  if (length(nameVecs) > 1 ){
    for (i in 2: length(nameVecs)){
      ans <- paste(ans, nameVecs[[i]], sep="_")
    }
  }
  ans
}

# Creates a character index based on a given set of vectors, where all possible combinations of the vectors are provided
make_index_comb = function( ...){
  nameVecs <- do.call(expand.grid,list(...)) %>%
    map(as.character)
  ans = nameVecs[[1]]
  if (length(nameVecs) > 1 ){
    for (i in 2: length(nameVecs)){
      ans <- paste(ans, nameVecs[[i]], sep="_")
    }
  }
  ans
}

#' Append an item to a list or vector
#' 
#' @param l A list of vector
#' @param a The item to append to the list
list_append = function(l, a){
  if (is.null(l)){
    l=a
  } else {
    N = length(l)
    l[[N+1]] = a
  }
  l
}


# Get incremental ranks
#' @param v vector
#' @param m matrix. m and v must have the same number of rows.
#' @result the percentage rank of each element in v relative to each row in m.  
vector_rank = function(v,m){
  if (length(v) != nrow(m)){
    stop("m and v must have the same number of rows")
  }
  N= length(v)
  M=ncol(m)
  rank_val = rep(0,N)
  for (i in 1:N){
    rank_val[i] = (sum(v[i]>m[i,]) + 0.5*(sum(v[i]==m[i,])))/M
  }
  rank_val
}

  

#' Calculation of Distances
#' 
#' @param x1, x2 input matrices
#' @return a list of distance difference to be fed into Gaussian Process model=
distance_calc = function(x1, x2=NULL){
  x1 = as.matrix(x1)
  if (is.null(x2))
    x2=x1
  x2=as.matrix(x2)
  D = ncol(x1)
  if (D != ncol(x2)){
    stop("x1 and x2 must have the same number of columns")
  }
  N1= nrow(x1)
  N2= nrow(x2)
  
  xd1 = array(0, dim = c(D, N1, N1))
  xd2 = array(0, dim = c(D, N2, N2))
  xd12 = array(0, dim = c(D, N1, N2))
  
  square_diff = function(a,b){
    outer(a,b, function(a,b) (a-b)^2/2)
  }
  
  for(d in 1:D) {
    xd1[d,,] = square_diff(x1[,d], x1[,d])
    xd2[d,,] = square_diff(x2[,d], x2[,d])
    xd12[d,,] = square_diff(x1[,d], x2[,d])
    
    #################
    # for (i in 1:N1){
    #   xd1[d, i, i] = 0;
    #   if (i < N1){
    #     for (j in (i+1):N1) {
    #       xd1[d, i, j] = -(x1[i,d]-x1[j,d])^2 / 2;
    #       xd1[d, j, i] = xd1[d, i, j];
    #     }  
    #   }
    # }
    
    # for (i in 1:N2){
    #   xd2[d, i, i] = 0;
    #   if (i < N2){
    #     for (j in (i+1):N2) {
    #       xd2[d, i, j] = -(x2[i,d]-x2[j,d])^2 / 2;
    #       xd2[d, j, i] = xd2[d, i, j];
    #     }
    #   }
    #   
    # }
    # 
    # for (i in 1:N1){
    #   #xd12[d, i, i] = 0;
    #   for (j in 1:N2) {
    #     xd12[d, i, j] = -(x1[i,d]-x2[j,d])^2 / 2;
    #   }
    # }
  }
  
  return(list(xd1=xd1, xd2=xd2, xd12=xd12))
}



#' Calculation of Linear factors for kernel
#' 
#' @param h1, h2 input matrices
#' @return a list of linear factors to be fed into Gaussian Process model=
linear_calc = function(h1, h2=NULL){
  h1 = as.matrix(h1)
  if (is.null(h2))
    h2=h1
  h2=as.matrix(h2)
  
  min_vals = apply(h1,2,min)
  h1=t(t(h1)-min_vals)
  h2=t(t(h2)-min_vals)
  
  
  D = ncol(h1)
  if (D != ncol(h2)){
    stop("h1 and h2 must have the same number of columns")
  }
  N1= nrow(h1)
  N2= nrow(h2)
  
  hl1 = array(0, dim = c(D, N1, N1))
  hl2 = array(0, dim = c(D, N2, N2))
  hl12 = array(0, dim = c(D, N1, N2))
  
  linear_val = function(a,b){
    matrix(a,ncol=1) %*% matrix(b, nrow=1)
  } 
  
  for(d in 1:D) {
    hl1[d,,] = linear_val(h1[,d], h1[,d])
    hl2[d,,] = linear_val(h2[,d], h2[,d])
    hl12[d,,] = linear_val(h1[,d], h2[,d])
  }
    
  return(list(hl1=hl1, hl2=hl2, hl12=hl12))
}



#' Calculation of Distances
#' 
#' @param x1, x2 input matrices
#' @return a list of distance difference to be fed into Gaussian Process model=
distance_calc_old = function(x1, x2=NULL){
  x1 = as.matrix(x1)
  if (is.null(x2))
    x2=x1
  x2=as.matrix(x2)
  D = ncol(x1)
  if (D != ncol(x2)){
    stop("x1 and x2 must have the same number of columns")
  }
  N1= nrow(x1)
  N2= nrow(x2)
  
  xd1 = array(0, dim = c(D, N1, N1))
  xd2 = array(0, dim = c(D, N2, N2))
  xd12 = array(0, dim = c(D, N1, N2))
  
  for(d in 1:D) {
    for (i in 1:N1){
      xd1[d, i, i] = 0;
      if (i < N1){
        for (j in (i+1):N1) {
          xd1[d, i, j] = (x1[i,d]-x1[j,d])^2 / 2
          xd1[d, j, i] = xd1[d, i, j];
        }  
      }
      
    }
  
  
    for (i in 1:N2){
      xd2[d, i, i] = 0;
      if (i < N2){
        for (j in (i+1):N2) {
          xd2[d, i, j] = (x2[i,d]-x2[j,d])^2 / 2
          xd2[d, j, i] = xd2[d, i, j];
        }
      }
      
    }
  
    for (i in 1:N1){
      for (j in 1:N2) {
        xd12[d, i, j] = (x1[i,d]-x2[j,d])^2 / 2;
      }
    }
  }
  
  return(list(xd1=xd1, xd2=xd2, xd12=xd12))
}



# Calculate a the 0.5 x square difference between  two vectorss
dist_calc = function(x1, x2=NULL){
  if (is.null(x2)){
    x2=x1
  }
  outer(x1,x2, FUN = function(a,b) 0.5 * (a-b) ^2)
}




# Calculate a matrix where elements are 1 if x!=y, 0 elsewhere 
cat_calc = function(x1, x2=NULL){
  if (is.null(x2)){
    x2=x1
  }
  outer(x1,x2, FUN = function(a,b) (a!=b) * 1)
}


#' Warps data using the beta function
#' 
#' @param x A Vector
#' @param alpha, beta Parameters of the beta used to warp the data
beta_warp = function(x, alpha, beta){
  pbeta(q = (x-min(x))/(max(x)-min(x)),alpha,beta)
}
