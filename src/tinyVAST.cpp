#define TMB_LIB_INIT R_init_tinyVAST
#include <TMB.hpp>   //Links in the TMB libraries
template<class Type>
Type objective_function<Type>::operator() (){
  using namespace density;
  using namespace Eigen; //Needed for utilisation of sparse structures

  // Data
  DATA_VECTOR(Y);        // The response
  DATA_MATRIX(X);        // Design matrix for fixed covariates
  DATA_MATRIX(Z);        // Design matrix for splines
  DATA_SPARSE_MATRIX(S); // Sparse penalization matrix
  DATA_IVECTOR(Sdims);   // Dimensions of blockwise components of S

  // Prediction options
  DATA_MATRIX(Xpred);       //Design matrix for fixed covariates
  DATA_UPDATE(Xpred);       //Design matrix for fixed covariates
  DATA_MATRIX(Zpred);       //Design matrix for splines
  DATA_UPDATE(Zpred);       //Design matrix for splines

  // Params
  PARAMETER_VECTOR(beta); // Fixed covariate parameters
  PARAMETER_VECTOR(gamma); // Spline regression parameters
  PARAMETER_VECTOR(log_lambda); //Penalization parameters
  PARAMETER(log_sigma);   

  // Transformations
  Type sigma = exp(log_sigma);
  vector<Type> lambda = exp(log_lambda);

  //Calculate the objective function--
  Type nll=0;

  // Likelihood of spline components
  int k=0;   // Counter
  for(int i=0;i<Sdims.size();i++){
    int m_i = Sdims(i);
    vector<Type> gamma_i = gamma.segment(k,m_i);       // Recover betai
    SparseMatrix<Type> S_i = S.block(k,k,m_i,m_i);  // Recover Si
    nll -= Type(0.5)*m_i*log_lambda(i) - 0.5*lambda(i)*GMRF(S_i).Quadform(gamma_i);
    k += m_i;
  }
  
  // Linear predictor and data
  vector<Type> mu = X*beta + Z*gamma;
  for(int i=0; i<Y.size(); i++){
    nll -= dnorm(Y(i), mu(i), sigma, true);
  }

  // Predictions
  vector<Type> mu_Xpred = Xpred*beta;
  vector<Type> mu_Zpred = Zpred*gamma;
  vector<Type> mu_pred = mu_Xpred + mu_Zpred;

  // Reporting
  ADREPORT(beta);
  REPORT( mu )
  REPORT( mu_pred )

  return nll;
}
