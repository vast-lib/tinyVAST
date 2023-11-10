#define TMB_LIB_INIT R_init_tinyVAST
#include <TMB.hpp>   //Links in the TMB libraries

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Sparse array * matrix
// NAs in IVECTOR or IMATRIX get converted to -2147483648, so isNA doesn't work as intended
template<class Type>
vector<Type> multiply_3d_sparse( matrix<int> A, vector<Type> weight, array<Type> x, int n_i ){
  vector<Type> out( n_i );
  out.setZero();
  for( int z=0; z<A.rows(); z++ ){
    out(A(z,0)) += weight(z) * x(A(z,1),A(z,2),A(z,3));
  }
  return out;
}

// get sign of double, only for REPORT use
template<class Type>
Type sign(Type x){
  return x / pow(pow(x,2),0.5);
}

// Deviance for the Tweedie
// https://en.wikipedia.org/wiki/Tweedie_distribution#Properties
template<class Type>
Type devresid_tweedie( Type y,
                       Type mu,
                       Type p ){

  Type c1 = pow( y, 2.0-p ) / (1.0-p) / (2.0-p);
  Type c2 = y * pow( mu, 1.0-p ) / (1.0-p);
  Type c3 = pow( mu, 2.0-p ) / (2.0-p);
  Type deviance = 2 * (c1 - c2 + c3 );
  Type devresid = sign( y - mu ) * pow( deviance, 0.5 );
  return devresid;
}

// Coding principles
// 1. Don't telescope model features using `map` (which is bug-prone), but
//    instead use 0-length vectors/matrices
// 2. Use Roman for data, Greek for parameters
// 3. Name objects with [var]_[indices] so that number of indices indicates
//    dimensionality
// 4. spatial_graph always has dimension 1+, but DSEM functionality is removed if
//    sem=NULL such that nrow(RAM)=0

template<class Type>
Type objective_function<Type>::operator() (){
  //using namespace R_inla;  // Not loaded globally, but used below
  using namespace density;
  using namespace Eigen;

  // Settings
  DATA_IMATRIX( f_ez );      // family/link
  DATA_IVECTOR( e_i );      // Asssociate each sample with a family/link
  DATA_IMATRIX( Edims_ez );  // Start (in C++ indexing) and Length, of log_sigma for each family/link

  // Data
  DATA_VECTOR( y_i );        // The response
  DATA_MATRIX( X_ij );        // Design matrix for fixed covariates
  DATA_MATRIX( Z_ik );        // Design matrix for splines
  DATA_IVECTOR( t_i );
  DATA_IVECTOR( c_i );
  DATA_SPARSE_MATRIX( S_kk ); // Sparse penalization matrix
  DATA_IVECTOR( Sdims );   // Dimensions of blockwise components of S_kk

  // Spatial objects
  DATA_INTEGER( spatial_method_code );   // Switch to string: https://kaskr.github.io/adcomp/compois_8cpp-example.html#a2
  DATA_IMATRIX( Aistc_zz );    // NAs get converted to -2147483648
  DATA_VECTOR( Axi_z );

  // SEM objects
  DATA_IMATRIX( RAM );
  DATA_VECTOR( RAMstart );

  // Prediction options
  DATA_MATRIX( X_gj );       // Design matrix for fixed covariates
  DATA_MATRIX( Z_gk );       // Design matrix for splines
  DATA_IMATRIX( Agstc_zz );       // Design matrix for SPDE projection (must be dense for DATA_UPDATE)
  DATA_VECTOR( Axg_z );
  DATA_IVECTOR( t_g );
  DATA_IVECTOR( c_g );
  DATA_IVECTOR( e_g );

  // Expansion options
  DATA_MATRIX( W_gz );            // Covariates for expansion
    // W_gz.col(0): Area
    // W_gz.col(1): Extra covariate, e.g., coordinates
  DATA_IMATRIX( V_gz );            // Settings for expansion
    // E_gz.col(0) : Expansion Type ( 0=sum(D*Area);  1=weighted.mean(W_gz.col(1),w=D*Area))
    // E_gz.col(1) : Prior row for bivariate weighting

  // Params
  PARAMETER( log_kappa );
  PARAMETER_VECTOR( alpha_j ); // Fixed covariate parameters
  PARAMETER_VECTOR( gamma_k ); // Spline regression parameters
  //PARAMETER_VECTOR(omega); // SPDE random effecs
  PARAMETER_VECTOR( beta_z ); // SEM coefficients
  PARAMETER_VECTOR( log_lambda ); //Penalization parameters
  PARAMETER_VECTOR( log_sigma );
  PARAMETER_VECTOR( delta0_c );
  PARAMETER_ARRAY( epsilon_stc );
  PARAMETER_VECTOR( eps );     // manual epsilon bias-correction, empty to turn off

  Type nll = 0;

  // Assemble precision
  int n_s = epsilon_stc.dim(0);
  int n_t = epsilon_stc.dim(1);
  int n_c = epsilon_stc.dim(2);
  int n_h = n_t * n_c;      // data
  Eigen::SparseMatrix<Type> Q_hh( n_h, n_h );
  // SEM
  Eigen::SparseMatrix<Type> Linv_hh(n_h, n_h);
  Eigen::SparseMatrix<Type> Rho_hh(n_h, n_h);
  Eigen::SparseMatrix<Type> Gammainv_hh(n_h, n_h);
  Eigen::SparseMatrix<Type> I_hh( n_h, n_h );
  Rho_hh.setZero();
  Gammainv_hh.setZero();
  I_hh.setIdentity();
  Type tmp;
  REPORT( n_h );
  for(int r=0; r<RAM.rows(); r++){
    // Extract estimated or fixed value
    if(RAM(r,3)>=1){
      tmp = beta_z(RAM(r,3)-1);
    }else{
      tmp = RAMstart(r);
    }
    if(RAM(r,0)==1) Rho_hh.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = tmp;
    if(RAM(r,0)==2) Gammainv_hh.coeffRef( RAM(r,1)-1, RAM(r,2)-1 ) = 1 / tmp;
  }
  Linv_hh = Gammainv_hh * ( I_hh - Rho_hh );
  Q_hh = Linv_hh.transpose() * Linv_hh;

  // Calculate effect of initial condition -- SPARSE version
  vector<Type> delta_h( n_h );
  delta_h.setZero();
  int h;
  if( delta0_c.size() > 0 ){
    // Compute delta_k
    matrix<Type> delta0_h1( n_h, 1 );
    delta0_h1.setZero();
    for(int c=0; c<n_c; c++){
      h = c * n_t;
      delta0_h1(h,0) = delta0_c(c);
    }

    // Sparse product
    Eigen::SparseMatrix<Type> IminusRho_hh = I_hh - Rho_hh;
    Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > lu;
    lu.compute(IminusRho_hh);
    matrix<Type> x = lu.solve(delta0_h1);

    // Resize
    delta_h = delta0_h1.array();
    REPORT( delta_h );
  }

  // Transformations
  vector<Type> lambda = exp(log_lambda);
  array<Type> epsilon_sh( n_s, n_h );
  for( int s=0; s<n_s; s++ ){
  for( int t=0; t<n_t; t++ ){
  for( int c=0; c<n_c; c++ ){
    h = c*n_t + t;
    epsilon_sh(s,h) = epsilon_stc(s,t,c);
  }}}

  // Spatial distribution
  Type log_tau = 0;
  Eigen::SparseMatrix<Type> Q_ss;
  if( (spatial_method_code==1) | (spatial_method_code==3) | (spatial_method_code==4) ){
    // Using INLA
    DATA_STRUCT(spatial_list, R_inla::spde_t);
    Q_ss = R_inla::Q_spde(spatial_list, exp(log_kappa));
    log_tau = log( 1.0 / (exp(log_kappa) * sqrt(4.0*M_PI)) );
  }else if( spatial_method_code==2 ){
    /// Using SAR
    DATA_SPARSE_MATRIX( Adj );
    Eigen::SparseMatrix<Type> I_ss( Adj.rows(), Adj.rows() );
    Eigen::SparseMatrix<Type> Lspatial_ss( Adj.rows(), Adj.rows() );
    I_ss.setIdentity();
    Lspatial_ss = ( I_ss - exp(log_kappa)*Adj );
    Q_ss = Lspatial_ss.transpose() * Lspatial_ss;
    log_tau = 0.0;
  }

  // GMRF:  non-separable time-variable, with separable spatace
  if( (n_s>0) & (n_h>0) ){ // PARALLEL_REGION
    // Including this line with Makevars below seems to cause a crash:
    // PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS)
    // PKG_CXXFLAGS=$(SHLIB_OPENMP_CXXFLAGS)
    nll += SEPARABLE( GMRF(Q_hh), GMRF(Q_ss) )( epsilon_sh );
  }

  // Distribution for spline components
  int k = 0;   // Counter
  for(int z=0; z<Sdims.size(); z++){   // PARALLEL_REGION
    int m_z = Sdims(z);
    vector<Type> gamma_segment = gamma_k.segment(k,m_z);       // Recover gamma_segment
    SparseMatrix<Type> S_block = S_kk.block(k,k,m_z,m_z);  // Recover S_i
    nll -= Type(0.5)*m_z*log_lambda(z) - 0.5*lambda(z)*GMRF(S_block).Quadform(gamma_segment);
    k += m_z;
  }

  // Linear predictor
  vector<Type> p_i = X_ij*alpha_j + Z_ik*gamma_k;
  p_i += multiply_3d_sparse( Aistc_zz, Axi_z, epsilon_stc, p_i.size() ) / exp(log_tau);
  for( int i=0; i<p_i.size(); i++ ){
    if( (n_h>0) ){     // (!isNA(c_i(i))) & (!isNA(t_i(i))) &
      h = c_i(i)*n_t + t_i(i);
      p_i(i) -= delta_h(h);
    }
  }
  // Likelihood
  vector<Type> devresid_i( y_i.size() );
  vector<Type> mu_i( y_i.size() );
  for( int i=0; i<y_i.size(); i++ ) {       // PARALLEL_REGION
      //vector<int> log_sigma_segment = Edims_ez.row(e_i(i));
    vector<Type> log_sigma_segment = log_sigma.segment( Edims_ez(e_i(i),0), Edims_ez(e_i(i),1) );
    switch( f_ez(e_i(i),0) ){
      case 0:  // Normal distribution
        mu_i(i) = p_i(i);
        nll -= dnorm( y_i(i), mu_i(i), exp(log_sigma_segment(0)), true );
        devresid_i(i) = y_i(i) - p_i(i);
        break;
      case 1: // Tweedie
        mu_i(i) = exp(p_i(i));
        nll -= dtweedie( y_i(i), mu_i(i), exp(log_sigma_segment(0)), 1.0 + invlogit(log_sigma_segment(1)), true );
        devresid_i(i) = devresid_tweedie( y_i(i), mu_i(i), 1.0 + invlogit(log_sigma_segment(1)) );
        break;
      default:
        error("Family not implemented.");
    }
  }

  // Predictions
  vector<Type> p_g = X_gj*alpha_j + Z_gk*gamma_k;
  vector<Type> mu_g( p_g.size() );
  p_g += multiply_3d_sparse( Agstc_zz, Axg_z, epsilon_stc, p_g.size() ) / exp(log_tau);
  for( int g=0; g<p_g.size(); g++ ){
    if( (n_h>0) ){     // (!isNA(c_i(i))) & (!isNA(t_i(i))) &
      h = c_g(g)*n_t + t_g(g);
      p_g(g) -= delta_h(h);
    }
  }
  for( int g=0; g<p_g.size(); g++ ){
    switch( f_ez(e_g(g),0) ){
      case 0: // Normal distribution
        mu_g(g) = p_g(g);
        break;
      case 1: // Tweedie
        mu_g(g) = exp(p_g(g));
        break;
    }
  }

  // Expansion
  if( (W_gz.rows()==mu_g.size()) & (V_gz.rows()==mu_g.size()) ){
    // First sweep
    vector<Type> phi0_g( mu_g.size() );
    for( int g=0; g<mu_g.size(); g++ ){
      if( (V_gz(g,0)==0) | (V_gz(g,0)==1) | (V_gz(g,0)==2) | (V_gz(g,0)==3) ){
        // Area-weighted average
        phi0_g(g) = mu_g(g) * W_gz(g,0);
      }
    }
    Type sumphi0 = sum(phi0_g);
    REPORT( phi0_g );

    // Second sweep for covariate or density-weighted averages
    vector<Type> phi_g( mu_g.size() );
    phi_g.setZero();
    for( int g=0; g<mu_g.size(); g++ ){
      if( V_gz(g,0)==0 ){
        // Exclude from 2nd-sweep calculation
        phi_g(g) = 0;
      }
      if( V_gz(g,0)==1 ){
        // Default:  equal to 1st swepp
        phi_g(g) = phi0_g(g);
      }
      if( V_gz(g,0)==2 ){
        // density-weighted average of W_gz(g,1)
        phi_g(g) = (phi0_g(g) / sumphi0) * W_gz(g,1);
      }
      if( (V_gz(g,0)==3) & (V_gz(g,1)>=0) & (V_gz(g,1)<=g) ){
        // density-weighted average of prediction
        phi_g(g) = (phi0_g(V_gz(g,1)) / sumphi0) * mu_g(g);
      }
    }
    //Type Metric = sum(phi_g);
    Type Metric = newton::Tag( sum(phi_g) ); // Set lowrank tag on Metric = sum(exp(x))
    REPORT( phi_g );
    REPORT( Metric );
    ADREPORT( Metric );

    if( eps.size() == 1 ){
      nll += Metric * eps(0);
    }
  }

  // Reporting
  Type range = pow(8.0, 0.5) / exp( log_kappa );
  REPORT( range );
  REPORT( p_i );
  REPORT( mu_i );
  if(p_g.size()>0){
    REPORT( p_g );
    REPORT( mu_g );
  }
  //REPORT( Q_hh );
  //REPORT( Q_ss );
  REPORT( devresid_i );

  return nll;
}
