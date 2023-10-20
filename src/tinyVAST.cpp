#define TMB_LIB_INIT R_init_tinyVAST
#include <TMB.hpp>   //Links in the TMB libraries

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() (){
  //using namespace R_inla;
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
  DATA_STRUCT(spatial_list, R_inla::spde_t);
  DATA_IMATRIX(Aistc);
  DATA_VECTOR(Ax);

  // SEM objects
  DATA_IMATRIX( RAM );
  DATA_VECTOR( RAMstart );

  // Prediction options
  DATA_MATRIX(predX);       // Design matrix for fixed covariates
  //DATA_UPDATE(Xpred);       //
  DATA_MATRIX(predZ);       // Design matrix for splines
  //DATA_UPDATE(Zpred);       //
  DATA_IMATRIX(predAistc);       // Design matrix for SPDE projection (must be dense for DATA_UPDATE)
  DATA_VECTOR( predAx );
  //DATA_UPDATE(Apred);       //
  DATA_IVECTOR( predt );
  //DATA_UPDATE( tpred );
  DATA_IVECTOR( predc );
  //DATA_UPDATE( cpred );

  // Params
  PARAMETER(log_kappa);
  //PARAMETER(log_tau);
  PARAMETER_VECTOR(alpha); // Fixed covariate parameters
  PARAMETER_VECTOR(gamma); // Spline regression parameters
  //PARAMETER_VECTOR(omega); // SPDE random effecs
  PARAMETER_VECTOR(beta_z); // SEM coefficients
  PARAMETER_VECTOR(log_lambda); //Penalization parameters
  PARAMETER(log_sigma);
  PARAMETER_VECTOR( delta0_c );
  //PARAMETER_ARRAY( x_tc );
  PARAMETER_ARRAY( epsilon_stc );

  Type nll = 0;
  // Assemble precision
  int n_k = n_t * n_c;      // data
  int k;
  int n_s = epsilon_stc.rows();
  Eigen::SparseMatrix<Type> Q_kk( n_k, n_k );
  Type log_tau = log( 1.0 / (exp(log_kappa) * sqrt(4.0*M_PI)) );
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
  array<Type> epsilon_sk( n_s, n_k );
  for( int s=0; s<n_s; s++ ){
  for( int t=0; t<n_t; t++ ){
  for( int c=0; c<n_c; c++ ){
    k = c*n_t + t;
    epsilon_sk(s,k) = epsilon_stc(s,t,c);
  }}}
  Eigen::SparseMatrix<Type> Q_spatial = R_inla::Q_spde(spatial_list, exp(log_kappa));

  // distributions
  if( (n_s * n_k) > 0 ){
    nll += SEPARABLE( GMRF(Q_kk), GMRF(Q_spatial) )( epsilon_sk );
  }

  // Likelihood of spline components
  k = 0;   // Counter
  for(int i=0; i<Sdims.size(); i++){
    int m_i = Sdims(i);
    vector<Type> gamma_i = gamma.segment(k,m_i);       // Recover gamma_i
    SparseMatrix<Type> S_i = S.block(k,k,m_i,m_i);  // Recover S_i
    nll -= Type(0.5)*m_i*log_lambda(i) - 0.5*lambda(i)*GMRF(S_i).Quadform(gamma_i);
    k += m_i;
  }

  // Linear predictor
  //vector<Type> mu = X*alpha + Z*gamma + A_is*omega/exp(log_tau);
  vector<Type> mu = X*alpha + Z*gamma;
  vector<Type> devresid( Y.size() );
  // Sparse product of matrix A and 3D array epsilon_stc
  for( int row=0; row<Aistc.rows(); row++ ){
    if( (!isNA(Aistc(row,2))) & (!isNA(Aistc(row,3))) ){   // Ignoring epsilon_stc when `sem` not passed to fit(.)
      mu(Aistc(row,0)) += Ax(row) * epsilon_stc(Aistc(row,1),Aistc(row,2),Aistc(row,3)) / exp(log_tau);
    }
  }
  for(int i=0; i<Y.size(); i++){
    if( (n_k>0) ){     // (!isNA(c_i(i))) & (!isNA(t_i(i))) &
      k = c_i(i)*n_t + t_i(i);
      mu(i) -= delta_k(k);
    }
    nll -= dnorm(Y(i), mu(i), sigma, true);
    devresid(i) = Y(i) - mu(i);
  }

  // Predictions
  vector<Type> mu_pred = predX*alpha + predZ*gamma;
  // Sparse product of matrix A and 3D array epsilon_stc
  for( int row=0; row<predAistc.rows(); row++ ){
    if( (!isNA(predAistc(row,2))) & (!isNA(predAistc(row,3))) ){   // Ignoring epsilon_stc when `sem` not passed to fit(.)
      mu_pred(predAistc(row,0)) += predAx(row) * epsilon_stc(predAistc(row,1),predAistc(row,2),predAistc(row,3)) / exp(log_tau);
    }
  }
  for( int row=0; row<predAistc.rows(); row++ ){
    if( (n_k>0) ){     // (!isNA(c_i(i))) & (!isNA(t_i(i))) &
      k = predAistc(row,3)*n_t + predAistc(row,2);
      mu_pred(predAistc(row,0)) -= delta_k(k);
    }
  }

  // Reporting
  Type range = pow(8.0, 0.5) / exp( log_kappa );
  REPORT( range );
  ADREPORT(alpha);
  REPORT( mu );
  REPORT( mu_pred );
  REPORT( Q_kk );
  REPORT( Q_spatial );
  //REPORT( Q_joint );
  REPORT( devresid );

  return nll;
}
