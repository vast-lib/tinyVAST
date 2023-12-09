#define TMB_LIB_INIT R_init_tinyVAST
#include <TMB.hpp>   //Links in the TMB libraries

// Needed for returning SparseMatrix for Ornstein-Uhlenbeck network correlations
template<class Type>
Eigen::SparseMatrix<Type> Q_network( Type log_theta,
                                     int n_s,
                                     vector<int> parent_s,
                                     vector<int> child_s,
                                     vector<Type> dist_s ){

  Eigen::SparseMatrix<Type> Q( n_s, n_s );
  Type theta = exp( log_theta );
  for(int s=0; s<n_s; s++){
    Q.coeffRef( s, s ) = Type(1.0);
  }
  for(int s=1; s<parent_s.size(); s++){
    if( exp(-dist_s(s))!=0 ){
      Q.coeffRef( parent_s(s), child_s(s) ) = -exp(-theta*dist_s(s)) / (1-exp(-2*theta*dist_s(s)));
      Q.coeffRef( child_s(s), parent_s(s) ) = Q.coeffRef( parent_s(s), child_s(s) );
      Q.coeffRef( parent_s(s), parent_s(s) ) += exp(-2*theta*dist_s(s)) / (1-exp(-2*theta*dist_s(s)));
      Q.coeffRef( child_s(s), child_s(s) ) += exp(-2*theta*dist_s(s)) / (1-exp(-2*theta*dist_s(s)));
    }
  }
  return Q;
}

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
template<class Type>
vector<Type> multiply_2d_sparse( matrix<int> A, vector<Type> weight, array<Type> x, int n_i ){
  vector<Type> out( n_i );
  out.setZero();
  for( int z=0; z<A.rows(); z++ ){
    out(A(z,0)) += weight(z) * x(A(z,1),A(z,2));
  }
  return out;
}

// get sign of double, only for REPORT use
template<class Type>
Type sign(Type x){
  return x / pow(pow(x,2),0.5);
}

// dlnorm
template<class Type>
Type dlnorm( Type x,
             Type meanlog,
             Type sdlog,
             int give_log=0){

  //return 1/(sqrt(2*M_PI)*sd) * exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
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
//    instead use 0-length vectors/matrices.  However, still retaining map = list()
//    to allow future uses.
// 2. Use Roman for data, Greek for parameters
// 3. Name objects with [var]_[indices] so that number of indices indicates
//    dimensionality
// 4. spatial_graph always has dimension 1+, but DSEM functionality is removed if
//    sem=NULL such that nrow(ram_dsem)=0

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
  DATA_VECTOR( offset_i );
  DATA_SPARSE_MATRIX( S_kk ); // Sparse penalization matrix
  DATA_IVECTOR( Sdims );   // Dimensions of blockwise components of S_kk

  // Spatial objects
  DATA_IVECTOR( spatial_options );   //
  // spatial_options(0)==1: SPDE;  spatial_options(0)==2: SAR;  spatial_options(0)==3: Off;  spatial_options(0)==4: stream-network
  DATA_IMATRIX( Aepsilon_zz );    // NAs get converted to -2147483648
  DATA_VECTOR( Aepsilon_z );
  DATA_IMATRIX( Aomega_zz );    // NAs get converted to -2147483648
  DATA_VECTOR( Aomega_z );

  // DSEM objects
  DATA_IMATRIX( ram_dsem );
  DATA_VECTOR( ram_dsem_start );

  // SEM objects
  DATA_IMATRIX( ram_sem );
  DATA_VECTOR( ram_sem_start );

  // Prediction options
  DATA_MATRIX( X_gj );       // Design matrix for fixed covariates
  DATA_MATRIX( Z_gk );       // Design matrix for splines
  DATA_IMATRIX( AepsilonG_zz );       // Design matrix for SPDE projection (must be dense for DATA_UPDATE)
  DATA_VECTOR( AepsilonG_z );
  DATA_IMATRIX( AomegaG_zz );       // Design matrix for SPDE projection (must be dense for DATA_UPDATE)
  DATA_VECTOR( AomegaG_z );
  DATA_IVECTOR( t_g );
  DATA_IVECTOR( c_g );
  DATA_VECTOR( offset_g );
  DATA_IVECTOR( e_g );

  // Expansion options
  DATA_MATRIX( W_gz );            // Covariates for expansion
    // W_gz.col(0): Area
    // W_gz.col(1): Extra covariate, e.g., coordinates
  DATA_IMATRIX( V_gz );            // Settings for expansion
    // E_gz.col(0) : Expansion Type ( 0=sum(D*Area);  1=weighted.mean(W_gz.col(1),w=D*Area))
    // E_gz.col(1) : Prior row for bivariate weighting

  // Params
  PARAMETER_VECTOR( alpha_j ); // Fixed covariate parameters
  PARAMETER_VECTOR( gamma_k ); // Spline regression parameters
  PARAMETER_VECTOR( beta_z ); // DSEM coefficients
  PARAMETER_VECTOR( theta_z ); // SEM coefficients
  PARAMETER_VECTOR( log_lambda ); //Penalization parameters
  PARAMETER_VECTOR( log_sigma );
  PARAMETER_VECTOR( delta0_c );
  PARAMETER_ARRAY( epsilon_stc );
  PARAMETER_ARRAY( omega_sc );
  PARAMETER_VECTOR( eps );     // manual epsilon bias-correction, empty to turn off

  // Globals
  Type nll = 0;
  Type tmp;

  // Assemble precision
  int n_s = epsilon_stc.dim(0);
  int n_t = epsilon_stc.dim(1);
  int n_c = epsilon_stc.dim(2);
  int n_h = n_t * n_c;      // data
  int h;

  // DSEM
  Eigen::SparseMatrix<Type> Q_hh( n_h, n_h );
  Eigen::SparseMatrix<Type> Linv_hh(n_h, n_h);
  Eigen::SparseMatrix<Type> Rho_hh(n_h, n_h);
  Eigen::SparseMatrix<Type> Gammainv_hh(n_h, n_h);
  Eigen::SparseMatrix<Type> Gamma_hh(n_h, n_h);
  Eigen::SparseMatrix<Type> I_hh( n_h, n_h );
  Rho_hh.setZero();
  Gammainv_hh.setZero();
  Gamma_hh.setZero();
  I_hh.setIdentity();
  for(int r=0; r<ram_dsem.rows(); r++){
    // Extract estimated or fixed value
    if(ram_dsem(r,3)>=1){
      tmp = beta_z(ram_dsem(r,3)-1);
    }else{
      tmp = ram_dsem_start(r);
    }
    if(ram_dsem(r,0)==1) Rho_hh.coeffRef( ram_dsem(r,1)-1, ram_dsem(r,2)-1 ) = tmp;
    if(ram_dsem(r,0)==2){
      Gammainv_hh.coeffRef( ram_dsem(r,1)-1, ram_dsem(r,2)-1 ) = 1 / tmp;
      Gamma_hh.coeffRef( ram_dsem(r,1)-1, ram_dsem(r,2)-1 ) = tmp;
    }
  }
  REPORT( Gamma_hh );
  REPORT( Rho_hh );

  // SEM
  Eigen::SparseMatrix<Type> Q_cc( n_c, n_c );
  Eigen::SparseMatrix<Type> Linv_cc(n_c, n_c);
  Eigen::SparseMatrix<Type> Rho_cc(n_c, n_c);
  Eigen::SparseMatrix<Type> Gammainv_cc(n_c, n_c);
  Eigen::SparseMatrix<Type> Gamma_cc(n_c, n_c);
  Eigen::SparseMatrix<Type> I_cc( n_c, n_c );
  Rho_cc.setZero();
  Gammainv_cc.setZero();
  Gamma_cc.setZero();
  I_cc.setIdentity();
  for(int r=0; r<ram_sem.rows(); r++){
    // Extract estimated or fixed value
    if(ram_sem(r,3)>=1){
      tmp = theta_z(ram_sem(r,3)-1);
    }else{
      tmp = ram_sem_start(r);
    }
    if(ram_sem(r,0)==1) Rho_cc.coeffRef( ram_sem(r,1)-1, ram_sem(r,2)-1 ) = tmp;
    if(ram_sem(r,0)==2){
      Gammainv_cc.coeffRef( ram_sem(r,1)-1, ram_sem(r,2)-1 ) = 1 / tmp;
      Gamma_cc.coeffRef( ram_sem(r,1)-1, ram_sem(r,2)-1 ) = tmp;
    }
  }
  REPORT( Gamma_cc );
  REPORT( Rho_cc );

  // Calculate effect of initial condition -- SPARSE version
  // Where does x go later?
  vector<Type> delta_h( n_h );
  delta_h.setZero();
  if( delta0_c.size() > 0 ){
    error("delta0 not currently working.");
    //// Compute delta_k
    //matrix<Type> delta0_h1( n_h, 1 );
    //delta0_h1.setZero();
    //for(int c=0; c<n_c; c++){
    //  h = c * n_t;
    //  delta0_h1(h,0) = delta0_c(c);
    //}
    //
    //// Sparse product
    //matrix<Type> x = inverseIminusRho_hh.solve(delta0_h1);
    //
    //// Resize
    //delta_h = delta0_h1.array();
    //REPORT( delta_h );
  }

  // Spatial distribution
  Type log_tau = 0;
  Eigen::SparseMatrix<Type> Q_ss;
  if( (spatial_options(0)==1) ){
    // Using INLA
    PARAMETER( log_kappa );
    DATA_STRUCT(spatial_list, R_inla::spde_t);
    Q_ss = R_inla::Q_spde(spatial_list, exp(log_kappa));
    log_tau = log( 1.0 / (exp(log_kappa) * sqrt(4.0*M_PI)) );
    Type range = pow(8.0, 0.5) / exp( log_kappa );
    REPORT( range );
  }else if( spatial_options(0)==2 ){
    /// Using SAR
    PARAMETER( log_kappa );
    DATA_SPARSE_MATRIX( Adj );
    Eigen::SparseMatrix<Type> I_ss( Adj.rows(), Adj.rows() );
    Eigen::SparseMatrix<Type> Lspatial_ss( Adj.rows(), Adj.rows() );
    I_ss.setIdentity();
    Lspatial_ss = ( I_ss - exp(log_kappa)*Adj );
    Q_ss = Lspatial_ss.transpose() * Lspatial_ss;
    log_tau = 0.0;
  }else if( spatial_options(0)==3 ){
    // Off, but using INLA inputs
    DATA_STRUCT(spatial_list, R_inla::spde_t);
    Q_ss = R_inla::Q_spde(spatial_list, Type(1.0));
    log_tau = Type(0.0);
  }else if( spatial_options(0)==4 ){
    // stream-network
    PARAMETER( log_kappa );
    DATA_IMATRIX( graph_sz );
    DATA_VECTOR( dist_s );
    Q_ss = Q_network( log_kappa, n_s, graph_sz.col(0), graph_sz.col(1), dist_s );  // Q_network( log_theta, n_s, parent_s, child_s, dist_s )
    log_tau = 0.0;
    REPORT( Q_ss );
  }

  // Space-variable interaction
  if( omega_sc.size()>0 ){ // PARALLEL_REGION
    if( spatial_options(1) == 0 ){
      // Separable precision
      Linv_cc = Gammainv_cc * ( I_cc - Rho_cc );
      Q_cc = Linv_cc.transpose() * Linv_cc;
      REPORT( Q_cc );

      // GMRF for SEM:  separable variable-space
      nll += SEPARABLE( GMRF(Q_cc), GMRF(Q_ss) )( omega_sc );
      // Including this line with Makevars below seems to cause a crash:
      // PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS)
      // PKG_CXXFLAGS=$(SHLIB_OPENMP_CXXFLAGS)
    }else{
      // Rank-deficient (projection) method
      Eigen::SparseMatrix<Type> I_cc( n_c, n_c );
      I_cc.setIdentity();
      nll += SEPARABLE( GMRF(I_cc), GMRF(Q_ss) )( omega_sc );

      // Sparse inverse-product
      Eigen::SparseMatrix<Type> IminusRho_cc = I_cc - Rho_cc;
      Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > inverseIminusRho_cc;
      inverseIminusRho_cc.compute(IminusRho_cc);

      // (I-Rho)^{-1} * Gamma * Epsilon
      matrix<Type> omega2_cs = Gamma_cc * omega_sc.matrix().transpose();
      matrix<Type> omega3_cs = inverseIminusRho_cc.solve(omega2_cs);
      omega_sc = omega3_cs.transpose();
      REPORT( omega_sc );
    }
  }

  // Space-time-variable interaction
  if( epsilon_stc.size()>0 ){ // PARALLEL_REGION
    // Reshape for either spatial_options
    array<Type> epsilon_hs( n_h, n_s );
    for( int s=0; s<n_s; s++ ){
    for( int t=0; t<n_t; t++ ){
    for( int c=0; c<n_c; c++ ){
      h = c*n_t + t;
      epsilon_hs(h,s) = epsilon_stc(s,t,c);
    }}}

    if( spatial_options(1) == 0 ){
      // Separable precision
      Linv_hh = Gammainv_hh * ( I_hh - Rho_hh );
      Q_hh = Linv_hh.transpose() * Linv_hh;

      // GMRF for DSEM:  non-separable time-variable, with separable space
      nll += SEPARABLE( GMRF(Q_ss), GMRF(Q_hh) )( epsilon_hs );
      REPORT( Q_hh );
    }else{
      // Rank-deficient (projection) method
      Eigen::SparseMatrix<Type> I_hh( n_h, n_h );
      I_hh.setIdentity();
      nll += SEPARABLE( GMRF(Q_ss), GMRF(I_hh) )( epsilon_hs );

      // Sparse inverse-product
      Eigen::SparseMatrix<Type> IminusRho_hh = I_hh - Rho_hh;
      Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > inverseIminusRho_hh;
      inverseIminusRho_hh.compute(IminusRho_hh);

      // (I-Rho)^{-1} * Gamma * Epsilon
      matrix<Type> e2_hs = Gamma_hh * epsilon_hs.matrix();
      matrix<Type> e3_hs = inverseIminusRho_hh.solve(e2_hs);

      // Transformations
      for( int s=0; s<n_s; s++ ){
      for( int t=0; t<n_t; t++ ){
      for( int c=0; c<n_c; c++ ){
        h = c*n_t + t;
        epsilon_stc(s,t,c) = e3_hs(h,s);
      }}}
      REPORT( epsilon_stc );
    }
  }

  // Distribution for spline components
  int k = 0;   // Counter
  for(int z=0; z<Sdims.size(); z++){   // PARALLEL_REGION
    int m_z = Sdims(z);
    vector<Type> gamma_segment = gamma_k.segment(k,m_z);       // Recover gamma_segment
    SparseMatrix<Type> S_block = S_kk.block(k,k,m_z,m_z);  // Recover S_i
    nll -= Type(0.5)*m_z*log_lambda(z) - 0.5*exp(log_lambda(z))*GMRF(S_block).Quadform(gamma_segment);
    k += m_z;
  }

  // Linear predictor
  vector<Type> p_i = X_ij*alpha_j + Z_ik*gamma_k + offset_i;
  p_i += multiply_3d_sparse( Aepsilon_zz, Aepsilon_z, epsilon_stc, p_i.size() ) / exp(log_tau);
  p_i += multiply_2d_sparse( Aomega_zz, Aomega_z, omega_sc, p_i.size() ) / exp(log_tau);
  for( int i=0; i<p_i.size(); i++ ){
    if( (n_h>0) ){     // (!isNA(c_i(i))) & (!isNA(t_i(i))) &
      h = c_i(i)*n_t + t_i(i);
      p_i(i) -= delta_h(h);
    }
  }
  // Likelihood
  vector<Type> devresid_i( y_i.size() );
  vector<Type> mu_i( y_i.size() );
  vector<Type> logmu_i( y_i.size() );
  for( int i=0; i<y_i.size(); i++ ) {       // PARALLEL_REGION
    vector<Type> log_sigma_segment = log_sigma.segment( Edims_ez(e_i(i),0), Edims_ez(e_i(i),1) );
    // Link function
    switch( f_ez(e_i(i),1) ){
      case 0:  // identity-link
        mu_i(i) = p_i(i);
        logmu_i(i) = log( p_i(i) );
        break;
      case 1: // log-link
        mu_i(i) = exp(p_i(i));
        logmu_i(i) = p_i(i);
        break;
      default:
        error("Link not implemented.");
    }
    if( !R_IsNA(asDouble(y_i(i))) ){
      // Distribution
      switch( f_ez(e_i(i),0) ){
        case 0:  // Normal distribution
          nll -= dnorm( y_i(i), mu_i(i), exp(log_sigma_segment(0)), true );
          devresid_i(i) = y_i(i) - p_i(i);
          break;
        case 1: // Tweedie
          nll -= dtweedie( y_i(i), mu_i(i), exp(log_sigma_segment(0)), 1.0 + invlogit(log_sigma_segment(1)), true );
          devresid_i(i) = devresid_tweedie( y_i(i), mu_i(i), 1.0 + invlogit(log_sigma_segment(1)) );
          break;
        case 2: // lognormal
          nll -= dlnorm( y_i(i), logmu_i(i) - 0.5*exp(2.0*log_sigma_segment(0)), exp(log_sigma_segment(0)), true );
          devresid_i(i) = log(y_i(i)) - ( logmu_i(i) - 0.5*exp(2.0*log_sigma_segment(0)) );
          break;
          break;
        case 3: // Poisson
          nll -= dpois( y_i(i), mu_i(i), true );
          // devresid_i(i) = MUST ADD;
          break;
        default:
          error("Distribution not implemented.");
      }
    }
  }

  // Predictions
  vector<Type> palpha_g = X_gj*alpha_j;
  vector<Type> pgamma_g = Z_gk*gamma_k;
  vector<Type> pepsilon_g = multiply_3d_sparse( AepsilonG_zz, AepsilonG_z, epsilon_stc, palpha_g.size() ) / exp(log_tau);
  vector<Type> pomega_g = multiply_2d_sparse( AomegaG_zz, AomegaG_z, omega_sc, palpha_g.size() ) / exp(log_tau);
  vector<Type> p_g = palpha_g + pgamma_g + offset_g + pepsilon_g + pomega_g;
  vector<Type> mu_g( p_g.size() );
  for( int g=0; g<p_g.size(); g++ ){
    if( (n_h>0) ){     // (!isNA(c_i(i))) & (!isNA(t_i(i))) &
      h = c_g(g)*n_t + t_g(g);
      p_g(g) -= delta_h(h);
    }
  }
  for( int g=0; g<p_g.size(); g++ ){
    switch( f_ez(e_g(g),1) ){
      case 0: // identity-link
        mu_g(g) = p_g(g);
        break;
      case 1: // log-link
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
  REPORT( p_i );
  REPORT( mu_i );
  if(p_g.size()>0){
    REPORT( p_g );
    REPORT( mu_g );
    REPORT( palpha_g );
    REPORT( pgamma_g );
    REPORT( pepsilon_g );
    REPORT( pomega_g );
    ADREPORT( p_g );
  }
  //REPORT( Q_ss );
  REPORT( devresid_i );
  REPORT( nll );

  return nll;
}
