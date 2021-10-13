// This version changes the likelihood to be (y(T) - R) ~ N(0, I).
// Error free
functions{
  //define transposed rotation matrix
  matrix rotation(real theta){
    matrix[2,2] rot;
    rot[1,1] = cos(theta);
    rot[1,2] = -sin(theta);
    rot[2,1] = -rot[1,2];
    rot[2,2] = rot[1,1];
    return rot;
  }
  
    matrix var_exp(row_vector[] x, real alpha, real rho, int N){
    matrix[N, N] K;
    real alpha_sq = square(alpha);

    for (i in 1: (N-1)){
      K[i,i] = alpha_sq;
      for( j in (i +1) :N){
        K[i, j] = alpha_sq * exp(-rho * distance(x[i], x[j]));
        K[j, i] = K[i, j];
      }
    }
    K[N,N] = alpha_sq;
    return K;
  }
  
    matrix cov_exp(row_vector[] x1, row_vector[] x2, real alpha, real rho, int N, int N_ref){
    matrix[N,N_ref] K;
    real alpha_sq = square(alpha);
    for (i in 1: N){
      for (j in 1 : N_ref){
        K[i,j] = alpha_sq * exp(-rho * distance(x1[i], x2[j]));
      }
    }
    return K;
  }
  
  row_vector[] m_to_rv(matrix M){
    int Ncol = cols(M);
    int Nrow = rows(M);
    row_vector[Nrow] rv[Ncol];
    
    for (c in 1: Ncol) rv[c] = to_row_vector(col(M, c));
    
    return rv;
  }
  real regularizer(vector reg){
    real reg_sum = sqrt(.5 * ((reg[1]^2) /3 + (reg[2]^2) /3 + reg[3]^2 + (reg[4]^2)/3 + (reg[5]^2) /3 + reg[6]^2 +
                        .5 * (reg[1]*reg[2] + reg[4]*reg[5]) + reg[1]*reg[3] + reg[2]*reg[3] + reg[4]*reg[6]+ reg[5]*reg[6] -
                        2.0/3.0 *(reg[1]+reg[5]) -.5*(reg[2] + reg[4])-(reg[3]+reg[6]) + 2.0/3.0));
    return reg_sum;
  }
  
  matrix Tmat(matrix s, real rot, vector translation){
    matrix[3,3] mat = diag_matrix(rep_vector(1, 3));
    
    mat[1:2, 1:2] = s * rotation(rot);
    mat[3, 1] = translation[1];
    mat[3, 2] = translation[2];
    
    return mat;
  }
  
    matrix Tmat_rv(matrix s, real rot, vector translation){
    matrix[3,3] mat = diag_matrix(rep_vector(1, 3));
    row_vector[2] c = to_row_vector(translation) * rotation(rot)*s;
    
    mat[1:2, 1:2] = rotation(rot)*s;
    
    mat[3, 1] = c[1];
    mat[3, 2] = c[2];
    return mat;
  }
}

 data{
      int<lower = 1> N;
      int NROW;
      int<lower = 1> N_ref;
      int N_rv;
      // int<lower = 1> Nsubj;//number of subject
      vector[N_ref] ref; // reference map
      vector[N] X; // activation map
      vector[N_rv] ref_rv;
      vector[N] X_rv;

      matrix[2, N] coord;
      matrix[2, N_ref] coord_ref;
      matrix[2, N_rv] coord_rv;
      
      // lower and upper bound of coord_ref
      real<lower=0> lim_var;
      real<lower=0> lim_theta[2];
      real<lower=0> lim_log_sigma[2];
      real limRota;

      // real<lower=0> nu;
      real<lower=0> rho;
      real b0;
      real tx0;
      real ty0;
      real lsx0;
      real lsy0;
      real r0;
      real phi0;
      real<lower = 0> lbd_b;
      real<lower = 0> lbd_T;
      
      real b0_rv;
      real tx0_rv;
      real ty0_rv;
      real lsx0_rv;
      real lsy0_rv;
      real r0_rv;
      real<lower = 0> lbd_b_rv;
      real<lower = 0> lbd_T_rv;
}

transformed data{
        vector[N] ones = rep_vector(1, N);
        real delta = 1e-9;
        // int N = NROW *NCOL;
        real lim_x_theta= lim_theta[1];
        real lim_y_theta = lim_theta[2];
        real lim_x_sigma= lim_log_sigma[1];
        real lim_y_sigma= lim_log_sigma[2];  
        real<lower = 0, upper = lim_var> alpha = 1;
        // real lim_b = b0 * .6;
        matrix[3,3] Tmat0 = Tmat(diag_matrix([exp(lsx0), exp(lsy0)]'), r0, [tx0, ty0]');
        matrix[3,3] Tmat0_rv = Tmat_rv(diag_matrix([exp(lsx0_rv), exp(lsy0_rv)]'), r0_rv, [tx0_rv, ty0_rv]');

        matrix[N,N] Sig= var_exp(m_to_rv(coord/NROW), alpha, rho, N) + diag_matrix(rep_vector(delta, N));// variance of reference
        matrix[N, N] L = cholesky_decompose(Sig); //Sig = L*L'
        vector[N] Sig_div_1 = mdivide_left_tri_low(L, ones);
        real oneT_SigInv_one = dot_self(Sig_div_1);
        vector[N] Sig_div_X = mdivide_left_tri_low(L, X);//L^{-1} * (ref-mu0) ~ normal(0,1)
        real mu0 = 1/oneT_SigInv_one * Sig_div_1' * Sig_div_X;
        vector[N] Sig_div_y = mdivide_left_tri_low(L, X-mu0);//L^{-1} * (ref-mu0) ~ normal(0,1)
        
        vector[N] Sig_div_X_rv = mdivide_left_tri_low(L, X_rv);//L^{-1} * (ref-mu0) ~ normal(0,1)
        real mu0_rv = 1/oneT_SigInv_one * Sig_div_1' * Sig_div_X_rv;
        vector[N] Sig_div_y_rv = mdivide_left_tri_low(L, X_rv-mu0_rv);//L^{-1} * (ref-mu0) ~ normal(0,1)
        
        Sig_div_y = mdivide_right_tri_low(Sig_div_y', L)'; //inv{L'*L} * (ref - mu0)
        Sig_div_y_rv = mdivide_right_tri_low(Sig_div_y_rv', L)'; 

}

parameters{
      real<lower = -lim_x_theta + tx0, upper = lim_x_theta + tx0> theta_x;
      real<lower = -lim_y_theta + ty0, upper = lim_y_theta + ty0> theta_y;
      
      //scaling parameters
      real<lower = -lim_x_sigma + lsx0, upper = lim_x_sigma + lsx0> log_sigma_x;
      real<lower = -lim_y_sigma + lsy0, upper = lim_y_sigma + lsy0> log_sigma_y;

      //rotation parameters
      real<lower =  r0 -limRota, upper = r0 + limRota> rota;
      
      // parameters for covariance matrix
      real<lower = 0, upper = 5> phi; //person specific error
      real<lower = 0>  b;
      
      real<lower = -lim_x_theta + tx0, upper = lim_x_theta + tx0> theta_x_rv;
      real<lower = -lim_y_theta + ty0, upper = lim_y_theta + ty0> theta_y_rv;
      
      //scaling parameters
      real<lower = -lim_x_sigma + lsx0, upper = lim_x_sigma + lsx0> log_sigma_x_rv;
      real<lower = -lim_y_sigma + lsy0, upper = lim_y_sigma + lsy0> log_sigma_y_rv;

      //rotation parameters
      real<lower =  r0 -limRota, upper = r0 + limRota> rota_rv;
      
      // parameters for covariance matrix
      real<lower = 0>  b_rv;

}

transformed parameters{
      // row_vector[2] coord_t[N_ref];  //concatenate transformed coordinates
      
      vector[N_ref] mu_kriged; //mean of the kriging 
      vector[N_rv] mu_kriged_rv;
      vector[N_ref] mu_kr; //mean of the kriging 
      vector[N_rv] mu_kr_rv; //mean of the kriging 

      matrix[N,N_ref] k_org_tf; // each column is covariance of transformed data to the original ones.
      vector[N_ref] ref_keep_transf = ref / sd(ref);
      
      matrix[N, N_rv] k_org_tf_rv; // each column is covariance of transformed data to the original ones.
      vector[N_rv] ref_keep_transf_rv = ref_rv / sd(ref);
      
      real log_lik[N_ref];
      real log_lik_rv[N_rv];

      matrix[2,2] scaling = diag_matrix([exp(log_sigma_x), exp(log_sigma_y)]');
      matrix[2,2] scaling_rv = diag_matrix([exp(log_sigma_x_rv), exp(log_sigma_y_rv)]');
      
      matrix[3, 3]regmat = Tmat(scaling, rota, [theta_x/NROW, theta_y/NROW]') - Tmat0;
      vector[6] reg = to_vector(regmat[1:3, 1:2]);
      real reg_sum = regularizer(reg);

      matrix[3, 3]regmat_rv = Tmat_rv(scaling_rv, rota_rv, [theta_x_rv/NROW, theta_y_rv/NROW]') - Tmat0_rv;
      vector[6] reg_rv = to_vector(regmat_rv[1:3, 1:2]);
      real reg_sum_rv = regularizer(reg_rv);
      
      matrix[2, N_ref] coord_t = scaling * rotation(rota) * coord_ref/NROW;
      matrix[2, N_rv] coord_t_rv = append_row(coord_rv[1]/NROW + theta_x_rv/NROW, coord_rv[2]/NROW + theta_y_rv/NROW);
      
      coord_t = append_row(coord_t[1] + theta_x/NROW, coord_t[2] + theta_y/NROW);
      coord_t_rv = rotation(rota_rv) * scaling_rv * coord_t_rv;
      
      // prepare for kriging: distance from transformed data to the original
      //l^th element in the transformed data 
      //m^th element in the original

      // kriging 
      k_org_tf = cov_exp(m_to_rv(coord/NROW), m_to_rv(coord_t), alpha, rho, N, N_ref);
      mu_kr  = k_org_tf' * Sig_div_y + mu0;
      mu_kriged = mu_kr / sd(ref) * b;
      
      k_org_tf_rv = cov_exp(m_to_rv(coord/NROW), m_to_rv(coord_t_rv), alpha, rho, N, N_rv);
      mu_kr_rv  = k_org_tf_rv' * Sig_div_y_rv + mu0_rv;
      mu_kriged_rv = mu_kr_rv / sd(ref) * b_rv;

      for (i in 1:N_ref) log_lik[i] = normal_lpdf(mu_kriged[i]| ref_keep_transf[i], phi);
      for (i in 1:N_rv) log_lik_rv[i] = normal_lpdf(mu_kriged_rv[i]| ref_keep_transf_rv[i], phi);
}

model{
      b ~ lognormal(log(b0), phi / sqrt(lbd_b * N_ref));
      reg_sum ~ normal(0, phi /  sqrt(lbd_T * N_ref));
      
      b_rv ~ lognormal(log(b0_rv), phi / sqrt(lbd_b_rv * N_ref));
      reg_sum_rv ~ normal(0, phi /  sqrt(lbd_T_rv * N_ref));
      
      target += sum(log_lik) + sum(log_lik_rv) - 0.5*log(phi);

    }
      
generated quantities{
      real<lower = 0> sigma_x;
      real<lower = 0> sigma_y;
      real<lower = 0> sigma_x_rv;
      real<lower = 0> sigma_y_rv;
      
      sigma_x = exp(log_sigma_x);
      sigma_y = exp(log_sigma_y);
      sigma_x_rv = exp(log_sigma_x_rv);
      sigma_y_rv = exp(log_sigma_y_rv);
}
