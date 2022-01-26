// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;

// helper function for detecting NAs in the data supplied from R
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// helper function to make sparse SPDE precision matrix
// Inputs:
//    logkappa: log(kappa) parameter value
//    logtau: log(tau) parameter value
//    M0, M1, M2: these sparse matrices are output from:
//     R::INLA::inla.spde2.matern()$param.inla$M*
template<class Type>
SparseMatrix<Type> spde_Q(Type logkappa, Type logtau, SparseMatrix<Type> M0,
                          SparseMatrix<Type> M1, SparseMatrix<Type> M2) {
  SparseMatrix<Type> Q;
  Type kappa2 = exp(2. * logkappa);
  Type kappa4 = kappa2*kappa2;
  Q = pow(exp(logtau), 2.)  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
  return Q;
}


// helper function to use the same penalized complexity prior on
//  matern params that is used in INLA

template<class Type>
Type dPCPriSPDE(Type logtau, Type logkappa,
                Type matern_par_a, Type matern_par_b,
                Type matern_par_c, Type matern_par_d,
                //vector<Type> matern_pri(4),
                int give_log=0)
{
  
  // matern_pri = c(a, b, c, d): P(range < a) = b; P(sigma > c) = d
  
  Type penalty; // prior contribution to jnll
  
  Type d = 2.;  // dimension
  Type lambda1 = -log(matern_par_b) * pow(matern_par_a, d/2.);
  Type lambda2 = -log(matern_par_d) / matern_par_c;
  Type range   = sqrt(8.0) / exp(logkappa);
  Type sigma   = 1.0 / sqrt(4.0 * 3.14159265359 * exp(2.0 * logtau) *
    exp(2.0 * logkappa));
  
  penalty = (-d/2. - 1.) * log(range) - lambda1 * pow(range, -d/2.) -
    lambda2 * sigma;
  // Note: (rho, sigma) --> (x=log kappa, y=log tau) -->
  //  transforms: rho = sqrt(8)/e^x & sigma = 1/(sqrt(4pi)*e^x*e^y)
  //  --> Jacobian: |J| propto e^(-y -2x)
  Type jacobian = - logtau - 2.0*logkappa;
  penalty += jacobian;
  
  if(give_log)return penalty; else return exp(penalty);
}

template<class Type>
Type robustMix(vector<Type> tmpW, vector<Type> tmpUrb, int nLen){
  Type maxVal = tmpUrb(0)+log(tmpW(0));
  for(int k = 0; k < nLen; ++k){
    if(maxVal < tmpUrb(k)+log(tmpW(k)))
      maxVal = tmpUrb(k)+log(tmpW(k));
  }
  Type tmpVal = 0.0;
  for(int k = 0; k < nLen; ++k){
    tmpVal += exp(tmpUrb(k)+log(tmpW(k))-maxVal);
  }
  tmpVal = maxVal + log(tmpVal);
  
  return(tmpVal);
}

///////////////////////////
// the main function     //
// to calculate the jnll //
///////////////////////////
template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // ~~~~~~~~~------------------------------------------------------~~
  // FIRST, we define params/values/data that will be passed in from R
  // ~~~~~~~~~~~------------------------------------------------------
  
  // normalization flag - used for speed-up
  DATA_INTEGER( flag1 ); // flag == 0 => no data contribution added to jnll
  // flag for choosing Gaussian or Binomial
  DATA_INTEGER( flag2 ); // flag2 == 0 => Gaussian
  
  // Indices
  DATA_INTEGER( num_iUrban );   // Number of urban data points in space
  DATA_INTEGER( num_iRural );   // Number of rural data points in space
  DATA_INTEGER( n_integrationPointsUrban );   // number of integration points for each observation
  DATA_INTEGER( n_integrationPointsRural );
  DATA_INTEGER( num_s );   // Number of mesh points in space mesh
  
  // Data (all except for X_ij is a vector of length num_i)
  DATA_VECTOR( y_iUrban );   // obs per binomial experiment at point i (clust)
  DATA_VECTOR( y_iRural );
  DATA_VECTOR( n_iUrban );   // Trials per cluster
  DATA_VECTOR( n_iRural );
  DATA_MATRIX( X_alphaUrban );  // 'design matrix' for just int
  DATA_MATRIX( X_alphaRural );
  DATA_MATRIX( wUrban ); // nObsUrban x nIntegrationPointsUrban weight matrix
  DATA_MATRIX( wRural ); // nObsRural x nIntegrationPointsRural weight matrix
  
  // SPDE objects
  DATA_SPARSE_MATRIX( M0 );
  DATA_SPARSE_MATRIX( M1 );
  DATA_SPARSE_MATRIX( M2 );
  DATA_SPARSE_MATRIX( AprojUrban );
  DATA_SPARSE_MATRIX( AprojRural );
  
  // Options
  DATA_VECTOR( options );
  // options[0] == 1 : use normalization trick
  // options[1] == 1 : adreport transformed params
  
  // Prior specifications
  DATA_VECTOR( alpha_pri );
  DATA_VECTOR( matern_pri);
  // matern_pri = c(a, b, c, d): P(range < a) = b; P(sigma > c) = d
  Type matern_par_a = matern_pri[0]; // range limit:    rho0
  Type matern_par_b = matern_pri[1]; // range prob:     alpha_rho
  Type matern_par_c = matern_pri[2]; // field sd limit: sigma0
  Type matern_par_d = matern_pri[3]; // field sd prob:  alpha_sigma
  
  // Fixed effects
  PARAMETER( alpha ); // Intercept
  // Log of INLA tau param (precision of space covariance matrix)
  PARAMETER( log_tau );
  // Log of INLA kappa (related to spatial correlation and range)
  PARAMETER( log_kappa );
  // Log of nugget standard deviation
  PARAMETER( log_nug_std );
  // Random effects for each spatial mesh vertex
  PARAMETER_VECTOR( Epsilon_s );
  
  // ~~~~~~~~~------------------------------------------------~~
  // SECOND, we define all other objects that we need internally
  // ~~~~~~~~~------------------------------------------------~~
  
  // objective function -- joint negative log-likelihood/posterior
  Type jnll = 0;
  
  // Make spatial precision matrix
  SparseMatrix<Type> Q_ss = spde_Q(log_kappa, log_tau, M0, M1, M2);
  
  // Transform some of our parameters
  Type sp_range = sqrt(8.0) / exp(log_kappa);
  Type sp_sigma = 1.0 / sqrt(4.0 * 3.14159265359 *
    exp(2.0 * log_tau) * exp(2.0 * log_kappa));
  
  // Define objects for derived values
  vector<Type> fe_iUrban(num_iUrban); // main effect: alpha
  vector<Type> fe_iRural(num_iRural);
  // Logit estimated prob for each cluster i
  //vector<Type> latent_field_iUrban(num_iUrban);
  //vector<Type> latent_field_iRural(num_iRural);
  // value of gmrf at data points
  vector<Type> projepsilon_iUrban(num_iUrban*n_integrationPointsUrban);
  vector<Type> projepsilon_iRural(num_iRural*n_integrationPointsRural);
  
  // fixed effects is just alpha in this example
  fe_iUrban = X_alphaUrban * Type(alpha); // initialize
  fe_iRural = X_alphaRural * Type(alpha);
  
  // Project GP approx from mesh points to data points
  // projepsilon_iUrban = AprojUrban * Epsilon_s.matrix();
  // projepsilon_iRural = AprojRural * Epsilon_s.matrix();
  projepsilon_iUrban = AprojUrban * Epsilon_s;
  projepsilon_iRural = AprojRural * Epsilon_s;
  
  // ~~~~~~~~~------------------------------------------------~~-
  // THIRD, we calculate the contribution to the likelihood from:
  // 1) GP field first, for 'flag' normalization purposes
  // 2) priors
  // 3) GP field
  // ~~~~~~~~~------------------------------------------------~~-
  
  /////////
  // (1) //
  /////////
  // the random effects. we do this first so to do the
  //   normalization outside of every optimization step
  // NOTE: likelihoods from namespace 'density' already return NEGATIVE
  //       log-liks so we add other likelihoods return positive log-liks
  if(options[0] == 1){
    // then we are not calculating the normalizing constant in the inner opt
    // that norm constant means taking an expensive determinant of Q_ss
    jnll += GMRF(Q_ss, false)(Epsilon_s);
    // return without data ll contrib to avoid unneccesary log(det(Q)) calcs
    if (flag1 == 0) return jnll;
  } else{
    jnll += GMRF(Q_ss)(Epsilon_s);
  }
  
  /////////
  // (2) //
  /////////
  // Prior contributions to joint likelihood (if options[1]==1)
  
  // add in priors for spde gp
  jnll -= dPCPriSPDE(log_tau, log_kappa,
                     matern_par_a, matern_par_b, matern_par_c, matern_par_d,
                     true);
  
  // prior for intercept
  jnll -= dnorm(alpha, alpha_pri[0], alpha_pri[1], true); // N(mean, sd)
  
  /////////
  // (3) //
  /////////
  // jnll contribution from each datapoint i
  
  // Precompute standard deviation
  Type nugStd = exp(log_nug_std);
  
  int thisIndex;
  Type thisLatentField;
  Type thisWeight;
  Type thislik;
  for (int i = 0; i < num_iUrban; i++){
    vector<Type> tmpUrb(n_integrationPointsUrban);
    vector<Type> tmpW(n_integrationPointsUrban);
    for (int j = 0; j < n_integrationPointsUrban; j++){
      thisIndex = i + num_iUrban * j;
      
      // latent field estimate at each obs
      thisLatentField = fe_iUrban(i) + projepsilon_iUrban(thisIndex);
      
      // and add data contribution to jnll
      if(!isNA(y_iUrban(i))){
        // get integration weight
        thisWeight = wUrban(i,j);
        
        if ( flag2 ==0 ){
          // Use Gaussian
          tmpW(j) = thisWeight;
          tmpUrb(j) = dnorm(y_iUrban(i), thisLatentField, nugStd, true);
          //thislik += thisWeight*dnorm(y_iUrban(i), thisLatentField, nugStd, false);
        }else{
          // Use dbinom_robust function, which takes the logit probability
          tmpW(j) = thisWeight;
          tmpUrb(j) = dbinom_robust( y_iUrban(i), n_iUrban(i), thisLatentField, true);
          //thislik += thisWeight*dbinom_robust( y_iUrban(i), n_iUrban(i), thisLatentField, false);
        }
        
        
      } // !isNA
      if(!isNA(y_iUrban(i))){
        thislik = robustMix(tmpW, tmpUrb, n_integrationPointsUrban);
        //thislik = tmpUrb[0];
      }else{
        thislik = 0.0;
      }
      
    } // for( j )
    
    jnll -= thislik;
  } // for( i )
  
  for (int i = 0; i < num_iRural; i++){
    vector<Type> tmpRur(n_integrationPointsRural);
    vector<Type> tmpW(n_integrationPointsRural);
    
    for (int j = 0; j < n_integrationPointsRural; j++){
      thisIndex = i + num_iRural * j;
      
      // latent field estimate at each obs
      thisLatentField = fe_iRural(i) + projepsilon_iRural(thisIndex);
      
      // and add data contribution to jnll
      if(!isNA(y_iRural(i))){
        // get integration weight
        thisWeight = wRural(i,j);
        
        if ( flag2 == 0 ){
          // Use Gaussian
          tmpW(j) = thisWeight;
          tmpRur(j) = dnorm(y_iRural(i), thisLatentField, nugStd, true);
          //thislik += thisWeight*dnorm(y_iRural(i), thisLatentField, nugStd, false);
        }else{
          // Use dbinom_robust function, which takes the logit probability
          tmpW(j) = thisWeight;
          tmpRur(j) = dbinom_robust( y_iRural(i), n_iRural(i), thisLatentField, true);
          //thislik += thisWeight*dbinom_robust( y_iRural(i), n_iRural(i), thisLatentField, false);
        }
        
      } // !isNA
      
    } // for( j )
    if(!isNA(y_iRural(i))){
      thislik = robustMix(tmpW, tmpRur, n_integrationPointsRural);
      //thislik = tmpRur[0];
    }else{
      thislik = 0.0;
    }
    
    jnll -= thislik;
  } // for( i )
  
  
  // ~~~~~~~~~~~
  // ADREPORT: used to return estimates and cov for transforms?
  // ~~~~~~~~~~~
  if(options[1]==1){
    ADREPORT(sp_range);
    ADREPORT(sp_sigma);
  }
  
  return jnll;
  
}

















