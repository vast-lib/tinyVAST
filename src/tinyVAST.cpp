#define TMB_LIB_INIT R_init_tinyVAST
#include <TMB.hpp>   //Links in the TMB libraries
template<class Type>
Type objective_function<Type>::operator() (){
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  // Data
  DATA_INTEGER( n_t );
  DATA_INTEGER( n_c );
  DATA_VECTOR(Y);        // The response
  DATA_MATRIX(X);        // Design matrix for fixed covariates
  DATA_MATRIX(Z);        // Design matrix for splines
  DATA_IVECTOR( t_i );
  DATA_IVECTOR( c_i );
  DATA_SPARSE_MATRIX(S); // Sparse penalization matrix
  DATA_IVECTOR(Sdims);   // Dimensions of blockwise components of S

  // Spatial objects
  DATA_STRUCT(spatial_list, spde_t);
  DATA_SPARSE_MATRIX(A_is);

  // SEM objects
  DATA_IMATRIX( RAM );
  DATA_VECTOR( RAMstart );

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
  PARAMETER_VECTOR(alpha); // Fixed covariate parameters
  PARAMETER_VECTOR(gamma); // Spline regression parameters
  PARAMETER_VECTOR(omega); // SPDE random effecs
  PARAMETER_VECTOR(beta_z); // SEM coefficients
  PARAMETER_VECTOR(log_lambda); //Penalization parameters
  PARAMETER(log_sigma);   
  PARAMETER_VECTOR( delta0_c );
  PARAMETER_ARRAY( x_tc );

  // Assemble precision
  int n_k = n_t * n_c;      // data
  Eigen::SparseMatrix<Type> Q_kk( n_k, n_k );
  // SEM
  Eigen::SparseMatrix<Type> Linv_kk(n_k, n_k);
  Eigen::SparseMatrix<Type> Rho_kk(n_k, n_k);
  Eigen::SparseMatrix<Type> Gammainv_kk(n_k, n_k);
  Eigen::SparseMatrix<Type> I_kk( n_k, n_k );
  Rho_kk.setZero();
  Gammainv_kk.setZero();
  I_kk.setIdentity();
  Type tmp;
  for(int r=0; r<RAM.rows(); r++){
    // Extract estimated or fixed value
    if(RAM(r,3)>=1){
      tmp = beta_z(RAM(r,3)-1);
    }else{
      tmp = RAMstart(r);
    }
    if(RAM(r,0)==1) Rho_kk.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = tmp;
    if(RAM(r,0)==2) Gammainv_kk.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = 1 / tmp;
  }
  Linv_kk = Gammainv_kk * ( I_kk - Rho_kk );
  Q_kk = Linv_kk.transpose() * Linv_kk;

  // Calculate effect of initial condition -- SPARSE version
  vector<Type> delta_k( n_k );
  delta_k.setZero();
  if( delta0_c.size() > 0 ){
    // Compute delta_k
    matrix<Type> delta0_k1( n_k, 1 );
    delta0_k1.setZero();
    int k;
    for(int c=0; c<n_c; c++){
      k = c * n_t;
      delta0_k1(k,0) = delta0_c(c);
    }

    // SPARSE version
    // See C:\Users\James.Thorson\Desktop\Work files\AFSC\2023-06 -- Sparse inverse-product\Kasper example\lu.cpp
    Eigen::SparseMatrix<Type> IminusRho_kk = I_kk - Rho_kk;
    Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > lu;
    lu.compute(IminusRho_kk);
    matrix<Type> x = lu.solve(delta0_k1);

    REPORT( delta0_k1 );
    delta_k = x.array();
  }
  REPORT( delta_k );

  // Transformations
  Type sigma = exp(log_sigma);
  vector<Type> lambda = exp(log_lambda);
  Eigen::SparseMatrix<Type> Q_spatial = Q_spde(spatial_list, exp(log_kappa));
  //Calculate the objective function
  Type nll = 0;
  nll += GMRF(Q_spatial)(omega);
  nll += GMRF(Q_kk)( x_tc - delta_k );

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
  vector<Type> mu = X*alpha + Z*gamma + A_is*omega/exp(log_tau);
  for(int i=0; i<Y.size(); i++){
    if( (n_t>0) & (n_c>0) ){
      mu(i) += x_tc( t_i(i), c_i(i) );
    }
    nll -= dnorm(Y(i), mu(i), sigma, true);
  }

  // Predictions
 vector<Type> mu_pred = Xpred*alpha + Zpred*gamma + Apred_is*omega/exp(log_tau);

  // Reporting
  ADREPORT(alpha);
  REPORT( mu );
  REPORT( mu_pred );
  REPORT( Q_kk );

  return nll;
}
