///  Gaussian Process Regression with Stan
///  Utilizes Cholesky Decomposition
///  Compound kernel
///  Allows for specification of a fixed prior mean
///  Heteroskedastic sigma
///  Hurdle Model at the lower bound
///  Vector prior mean
///  Multi company Model

functions {
  

  
  // mu, sigma = parameters of normal distribution
  // l_bound = left truncation value
  real normal_lb_rng(real mu, real sigma, real l_bound) {
    real lb = normal_cdf(l_bound, mu, sigma);  // cdf with lower bound
    real u;
    if (lb == 1){
      return (l_bound);
    } else {
      u = uniform_rng(lb, 1);
    }
    return (sigma * inv_Phi(u)) + mu;  // inverse cdf 
  }
  
  real hurdle_rng(real mu, real sigma, real l_bound) {
    real hurdle = normal_cdf(l_bound, mu, sigma);  // cdf with lower bound
    real u;
    u = uniform_rng(0, 1);
    if (u <= hurdle){
      return (l_bound);
    } else {
      return (sigma * inv_Phi(u)) + mu;  // inverse cdf 
    }
  }
  
  // Add a fixed constant to the add_diagonal for a square matrix
  matrix add_diagonal(matrix K, real delta){
    int N = rows(K);
    matrix[N,N] KK = K;
    
    for (i in 1:N){
      KK[i,i] = K[i,i] + delta;
    }
    
    return KK;
  } 

  // Add the sigma value to the diagonal of the covariance matrix
  matrix add_sigma(matrix K, vector sigma, int[] dev_lag){
    int N = rows(K);
    int dimSigma = rows(sigma);
    matrix[N,N] KK = K;
    
    for (i in 1:N){
      KK[i,i] = K[i,i] + sigma[dimSigma+1-dev_lag[i]];
    }
    
    return KK;
  } 

  
  // Square exponential kernel
  // no_input_vectors = 1 or 2: Use 1 if x1 and x2 are identical and are on the main diagonal (i.e.  for K, K**) 
  //                    otherwise use 2 (i.e. for K*)
  matrix cov_exp_quad_se(vector[] x1, vector[] x2, real eta,  vector rho, real delta, int no_input_vectors) {
    
    int N1 = size(x1);
    int N2 = size(x2);
    matrix[N1, N2] K;
    real eta_sq = eta ^2;
    
    
    if( no_input_vectors == 1){
      for (i in 1:(N1-1)) {
        K[i, i] = eta_sq   + delta;
        for (j in (i+1):N2) {
          K[i, j] = eta_sq * exp(-0.5 * dot_self((x1[i] - x2[j]) ./ rho));
          K[j, i] = K[i, j];
        }
      }
      K[N1, N1] = eta_sq   + delta;
    } else {
      for (i in 1:N1) {
        for (j in 1:N2) {
          K[i, j] = eta_sq * exp(-0.5 * dot_self((x1[i] - x2[j]) ./ rho));
        }
      }    
    }
    // if (delta >0) {
    //   K = add_diagonal(K, delta);
    // }
    return K;
  }
  
  // Linear kernel
  matrix cov_linear(vector h1, vector h2, real theta, real delta, int no_input_vectors) {
    
    int N1 = rows(h1);
    int N2 = rows(h2);
    
    matrix[N1, N2] K;
    
    if (no_input_vectors == 1){
      for (i in 1:N1) {
        for (j in i:N2) {
          K[i, j] = theta * h1[i] * h2[j];
          K[j, i] = K[i, j];
        }
      }
    } else {
      for (i in 1:N1) {
        for (j in 1:N2) {
          K[i, j] = theta * h1[i] * h2[j];
        }
      }
    }
    if (delta >0) {
      K = add_diagonal(K, delta);
    }
    return K;
  }
  
    // Calculates the categorical factor kernel. 
  // Calculates a matrix where each value is 1 if c1!=c2, 0 otherwise
  // c1, c2 categorical values.  
  matrix cov_categorical(int[] c1, int[] c2, real alpha){
    int n1 = size(c1);
    int n2 = size(c2);
    matrix[n1, n2] cat_mat = rep_matrix(-alpha, n1, n2);
    
    for (i1 in 1:n1){
      for (i2 in 1:n2){
        if(c1[i1] == c2[i2]){
          cat_mat[i1,i2] = 0;
        }
      }
    }
    
    return exp(cat_mat);
  }

  
  
  // Binary kernel (alternative method of including categorical variable, using a one-hot encoding structure  not used)
  // NOT USED
  matrix cov_binary(vector[] c1, vector[] c2, real alpha, real delta, int no_input_vectors) {
    
    int N1 = size(c1);
    int N2 = size(c2);
    matrix[N1, N2] K;
    
    if(no_input_vectors == 1){
      for (i in 1:N1) {
        for (j in i:N2) {
          K[i, j] = exp( - alpha * (1-(c1[i]' * c2[j])));
          K[j, i] = K[i, j];
        }
      }
    } else {
      for (i in 1:N1) {
        for (j in 1:N2) {
          K[i, j] = exp( - alpha * (1-(c1[i]' * c2[j])));
        }
      }    
    }
    if (delta >0) {
      K = add_diagonal(K, delta);
    }
    return K;
  }

  
  // Random Walk kernel
  matrix cov_random_walk(vector h1, vector h2, real theta, real delta, int no_input_vectors) {
    
    int N1 = rows(h1);
    int N2 = rows(h2);
    
    matrix[N1, N2] K;
    
    if (no_input_vectors == 1){
      for (i in 1:N1) {
        for (j in i:N2) {
          K[i, j] = theta * fmin(h1[i], h2[j]);
          K[j, i] = K[i, j];
        }
      }
    } else {
        for (i in 1:N1) {
          for (j in 1:N2) {
            K[i, j] = theta * fmin(h1[i], h2[j]);
          }
        }
    }
    
    if (delta >0) {
      K = add_diagonal(K, delta);
    }
    return K;
  }
  
  
  // General function to calculate an advance kernel
  // In this case, it just case the square exponential kernel
  matrix cov_fun(vector[] x1, vector[] x2, int[] x_cat1, int[] x_cat2, vector[] h1, vector[] h2,   
    real eta, vector rho, vector theta, real alpha, real delta, int no_input_vectors){
      matrix[size(x1), size(x2)] res;
      int N_h = size(h1);

      res = cov_exp_quad_se(x1, x2, eta,  rho, delta,no_input_vectors);
      
      for (i in 1:N_h){
        res = res +  cov_linear(h1[i], h2[i], theta[i],0,no_input_vectors);
      }
      
      res = res .* cov_categorical(x_cat1, x_cat2, alpha);

    return res;
  }
  
  
  // prediction function
  vector gp_pred_rng(vector[] x1, vector[] x2, int[] x_cat1, int[] x_cat2, vector[] h1, vector[] h2, int[] dev_lag,
                     vector y1,  vector mu_vec1, vector mu_vec2,
                     real eta, vector rho, vector sigma, vector theta, real alpha,
                     real delta) {
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] f2;
    {
      matrix[N1, N1] K =   add_sigma(cov_fun(x1, x1, x_cat1, x_cat1, h1, h1, eta, rho, theta,  alpha, delta, 1), 
                                      square(sigma),dev_lag);
      matrix[N1, N1] L_K = cholesky_decompose(K);

      vector[N1] L_K_div_y1 = mdivide_left_tri_low(L_K, y1-mu_vec1);
      vector[N1] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_K)';
      matrix[N1, N2] k_x1_x2 = cov_fun(x1, x2, x_cat1, x_cat2, h1, h2, eta, rho, theta, alpha, 0, 2);
      vector[N2] f2_mu = mu_vec2 + (k_x1_x2' * K_div_y1);
      matrix[N1, N2] v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      matrix[N2, N2] cov_f2 = cov_fun(x2, x2, x_cat2, x_cat2, h2, h2, eta, rho, theta, alpha, delta,1) -v_pred' * v_pred;
      f2 = multi_normal_rng(f2_mu, cov_f2);
    }
    return f2;
  }
}

data {
  int<lower=1> N1;  // historic data
  int<lower=1> N2;  // predicted points
  int<lower=1> D; // Dimension of input values
  int<lower=1> dimH; // dimension of basis functions
  int<lower=1> dimSigma; // dimension sigma

  // Historic Data
  vector[D] x1[N1];  // input values (N1 x D matrix from R)
  int x_cat1[N1];  // input values containing a categorical value
  
  vector[N1] H1[dimH];  // input values (N1 x dimH matrix from R)
  int dev_lag1[N1];  // development lag
  vector[N1] y1;  // output values 

  // Inputs for prediction
  vector[D] x2[N2];  // input values (N x D matrix from R)
  int x_cat2[N2];  // input values containing a categorical value
  vector[N2] H2[dimH];  // input values (N1 x dimH matrix from R)
  int dev_lag2[N2];  // development lag
 
 // inputs for priors on parameters
 vector[dimSigma] mu;
 real l_bound;
 real<lower=0> prior_eta;
 real prior_rho[2];
 real<lower=0> prior_sigma;
 real<lower=0> prior_theta;
 real<lower=0> prior_alpha;
 

}

transformed data{
  real delta = 0.000001;
  vector[N1] H1_pos[dimH];  
  vector[N2] H2_pos[dimH];  
  vector[N1] mu_vec1;
  vector[N2] mu_vec2;
  int below_hurdle[N1];
  int sigma_index[N1];

  for (n in 1:N1){
    sigma_index[n] = dimSigma+1-dev_lag1[n];
    if(y1[n] <= l_bound){
      below_hurdle[n] = 1;
    } else {
      below_hurdle[n] = 0;
    }
  }
  
  for (i in 1:dimH){
    H1_pos[i] = H1[i] - min(H1[i]);
    H2_pos[i] = H2[i] - min(H2[i]);
    
  }
  
  for (i in 1:N1){
    mu_vec1[i] = mu[dev_lag1[i]];
  }
  for (i in 1:N2){
    mu_vec2[i] = mu[dev_lag2[i]];
  }


}  

parameters {
  vector<lower=0>[D] rho;
  real<lower=0> eta;
  positive_ordered[dimSigma] sigma;
  vector<lower=0>[dimH] theta;
  real<lower=0> alpha;
  vector[N1] f;
}

transformed parameters {
  vector<lower=0>[D] rho_sq;
  real<lower=0> eta_sq;
  

  
  
  eta_sq = eta ^2;
  for (i in 1:D){
    rho_sq[i] = rho[i]^2;
  }
  
}

model {
  
  // Mean field function
  matrix[N1, N1] cov =   cov_fun(x1, x1, x_cat1, x_cat1, H1_pos, H1_pos,  eta, rho, theta,  alpha, delta, 1);
  matrix[N1, N1] L_cov;
  vector[N1] hurdle;
  L_cov = cholesky_decompose(cov);
  
  
  f ~ multi_normal_cholesky(mu_vec1, L_cov);
  

  // Priors
  rho ~ inv_gamma(prior_rho[1], prior_rho[2]);
  eta ~ normal(0, prior_eta);
  sigma ~ normal(0, prior_sigma);
  theta ~ normal(0, prior_theta);
  alpha ~ normal(0, prior_alpha);

  // Sampling Statement
  hurdle  = Phi((l_bound - f)./ sigma[sigma_index]); // hurdle probability
  below_hurdle ~ bernoulli(hurdle); // likelihood for hurdle
  
  // likelihood for normal distribution truncated to left at hurdle
  for (n in 1:N1){
    if(y1[n] > l_bound){
      target +=  normal_lpdf(y1[n] | f[n], sigma[sigma_index[n]]) - log(1-hurdle[n]);  //normal_lccdf(l_bound | f[n], sigma[sigma_index[n]]);
      //y1[n] ~ normal(f[n], sigma[dimSigma+1-dev_lag1[n]]) T[l_bound, ];
    }
  }
}

generated quantities {
  vector[N2] fStar = gp_pred_rng(x1, x2, x_cat1, x_cat2, H1_pos, H2_pos, dev_lag1, y1, mu_vec1, mu_vec2, 
                                 eta, rho, sigma, theta, alpha, delta);
  vector[N2] yStar;
  for (n in 1:N2)
    yStar[n] = hurdle_rng(fStar[n], sigma[dimSigma+1-dev_lag2[n]], l_bound);
}
