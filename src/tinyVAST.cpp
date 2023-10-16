#define TMB_LIB_INIT R_init_tinyVAST
#include <TMB.hpp>   //Links in the TMB libraries
template<class Type>
Type objective_function<Type>::operator() (){
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  // Data
  DATA_VECTOR(Y);        // The response
  DATA_MATRIX(X);        // Design matrix for fixed covariates
  DATA_MATRIX(Z);        // Design matrix for splines
  DATA_SPARSE_MATRIX(S); // Sparse penalization matrix
  DATA_IVECTOR(Sdims);   // Dimensions of blockwise components of S

  // Spatial objects
  DATA_STRUCT(spatial_list, spde_t);
  DATA_SPARSE_MATRIX(A_is);

  // Prediction options
  DATA_MATRIX(Xpred);       // Design matrix for fixed covariates
  DATA_UPDATE(Xpred);       //
  DATA_MATRIX(Zpred);       // Design matrix for splines
  DATA_UPDATE(Zpred);       //
  DATA_MATRIX(Apred_is);       // Design matrix for SPDE projection (must be dense for DATA_UPDATE)
  DATA_UPDATE(Apred_is);       //

  // Params
  PARAMETER(log_kappa);
  PARAMETER(log_tau);
  PARAMETER_VECTOR(beta); // Fixed covariate parameters
  PARAMETER_VECTOR(gamma); // Spline regression parameters
  PARAMETER_VECTOR(omega); // SPDE random effecs
  PARAMETER_VECTOR(log_lambda); //Penalization parameters
  PARAMETER(log_sigma);   

  // Transformations
  Type sigma = exp(log_sigma);
  vector<Type> lambda = exp(log_lambda);
  Eigen::SparseMatrix<Type> Q_spatial = Q_spde(spatial_list, exp(log_kappa));
  //Calculate the objective function
  Type nll = 0;
  nll += GMRF(Q_spatial)(omega);

  // Likelihood of spline components
  int k=0;   // Counter
  for(int i=0; i<Sdims.size(); i++){
    int m_i = Sdims(i);
    vector<Type> gamma_i = gamma.segment(k,m_i);       // Recover gamma_i
    SparseMatrix<Type> S_i = S.block(k,k,m_i,m_i);  // Recover S_i
    nll -= Type(0.5)*m_i*log_lambda(i) - 0.5*lambda(i)*GMRF(S_i).Quadform(gamma_i);
    k += m_i;
  }
  
  // Linear predictor and data
  vector<Type> mu = X*beta + Z*gamma + A_is*omega/exp(log_tau);
  for(int i=0; i<Y.size(); i++){
    nll -= dnorm(Y(i), mu(i), sigma, true);
  }

  // Predictions
 vector<Type> mu_pred = Xpred*beta + Zpred*gamma + Apred_is*omega/exp(log_tau);

  // Reporting
  ADREPORT(beta);
  REPORT( mu )
  REPORT( mu_pred )

  return nll;
}
