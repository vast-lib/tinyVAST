#define TMB_LIB_INIT R_init_tinyVAST
#include <TMB.hpp>   //Links in the TMB libraries
#include "utils.h"

// Coding principles
// 1. Don't telescope model features using `map` (which is bug-prone), but
//    instead use 0-length vectors/matrices.  However, still retaining map = list()
//    to allow future uses.
// 2. Use Roman for data, Greek for parameters
// 3. Name objects with [var]_[indices] so that number of indices indicates
//    dimensionality
// 4. spatial_graph always has dimension 1+, but DSEM functionality is removed if
//    sem=NULL such that nrow(ram_spacetime_term)=0

template<class Type>
Type objective_function<Type>::operator() (){
  //using namespace R_inla;  // Not loaded globally, but used below
  using namespace density;
  using namespace Eigen;
  using namespace tinyVAST;

  // Settings
  DATA_IMATRIX( family_ez );      // family
  DATA_IMATRIX( link_ez );        // family
  DATA_IVECTOR( components_e );   // family
  DATA_IVECTOR( poislink_e );     // family
  DATA_IVECTOR( e_i );      // Asssociate each sample with a family/link
  DATA_IMATRIX( Edims_ez );  // Start (in C++ indexing) and Length, of log_sigma for each family/link

  // Data
  DATA_VECTOR( y_i );        // The response
  DATA_MATRIX( X_ij );        // Design matrix for fixed covariates
  DATA_MATRIX( Z_ik );        // Design matrix for splines
  DATA_MATRIX( W_il );        // Design matrix for SVCs
  DATA_MATRIX( X2_ij );        // Design matrix for fixed covariates
  DATA_MATRIX( Z2_ik );        // Design matrix for splines
  DATA_MATRIX( W2_il );        // Design matrix for SVCs
  DATA_IVECTOR( t_i );
  DATA_IVECTOR( c_i );
  DATA_VECTOR( offset_i );
  DATA_VECTOR( weights_i );
  DATA_UPDATE( weights_i );  // Warning: cannot use branching (weights_i > 0) when using DATA_UPDATE
  DATA_VECTOR( size_i );  // Ignored unless family = binomial
  DATA_SPARSE_MATRIX( S_kk ); // Sparse penalization matrix
  DATA_IVECTOR( Sdims );   // Dimensions of blockwise components of S_kk
  DATA_IVECTOR( Sblock );
  DATA_SPARSE_MATRIX( S2_kk ); // Sparse penalization matrix
  DATA_IVECTOR( S2dims );   // Dimensions of blockwise components of S_kk
  DATA_IVECTOR( S2block );

  // Spatial objects
  DATA_IVECTOR( model_options );   //
  // model_options(0)==1: SPDE;  model_options(0)==2: SAR;  model_options(0)==3: Off;  model_options(0)==4: stream-network
  // model_options(1)==0: use GMRF(Q) to evaluate density;  model_options(1)==1: use GMRF(I) and project by Q^{-0.5} to evaluate density
  // model_options(2)==0: no RSS; model_options(2)==1: yes RSR
  // model_options(3)==0: no extra reporting; model_options(3)==1: yes extra reporting
  // model_options(4)==0: no SE for _g; model_options(4)==1: SE for p_g;  model_options(4)==2: SE for mu_g
  DATA_IMATRIX( Aepsilon_zz );    // NAs get converted to -2147483648
  DATA_VECTOR( Aepsilon_z );
  DATA_IMATRIX( Aomega_zz );    // NAs get converted to -2147483648
  DATA_VECTOR( Aomega_z );
  DATA_SPARSE_MATRIX( A_is );    // Used for SVC
  DATA_STRUCT( nngp_data, nngp_data_t );

  // DSEM objects
  DATA_IMATRIX( ram_space_term );
  DATA_VECTOR( ram_space_term_start );
  DATA_IMATRIX( ram_time_term );
  DATA_VECTOR( ram_time_term_start );
  DATA_IMATRIX( ram_spacetime_term );
  DATA_VECTOR( ram_spacetime_term_start );
  DATA_IMATRIX( ram2_space_term );
  DATA_VECTOR( ram2_space_term_start );
  DATA_IMATRIX( ram2_time_term );
  DATA_VECTOR( ram2_time_term_start );
  DATA_IMATRIX( ram2_spacetime_term );
  DATA_VECTOR( ram2_spacetime_term_start );

  // Prediction options
  DATA_MATRIX( X_gj );       // Design matrix for fixed covariates
  DATA_MATRIX( Z_gk );       // Design matrix for splines
  DATA_MATRIX( W_gl );
  DATA_MATRIX( X2_gj );       // Design matrix for fixed covariates
  DATA_MATRIX( Z2_gk );       // Design matrix for splines
  DATA_MATRIX( W2_gl );
  DATA_IMATRIX( AepsilonG_zz );       // Design matrix for SPDE projection (must be dense for DATA_UPDATE)
  DATA_VECTOR( AepsilonG_z );
  DATA_IMATRIX( AomegaG_zz );       // Design matrix for SPDE projection (must be dense for DATA_UPDATE)
  DATA_VECTOR( AomegaG_z );
  DATA_SPARSE_MATRIX( A_gs );
  DATA_IVECTOR( t_g );
  DATA_IVECTOR( c_g );
  DATA_VECTOR( offset_g );
  DATA_IVECTOR( e_g );

  // Expansion options
  DATA_MATRIX( W_gz );            // Covariates for expansion
    // W_gz.col(0): Area
    // W_gz.col(1): Extra covariate, e.g., coordinates
  DATA_IMATRIX( V_gz );            // Settings for expansion
    // V_gz.col(0) : Expansion Type ( 0=sum(D*Area);  1=weighted.mean(W_gz.col(1),w=D*Area))
    // V_gz.col(1) : Prior row for bivariate weighting
    // V_gz.col(2) : Blocks for derived quantities

  // Params
  PARAMETER_VECTOR( alpha_j ); // Fixed covariate parameters
  PARAMETER_VECTOR( gamma_k ); // Spline regression parameters
  PARAMETER_VECTOR( beta_z ); // DSEM coefficients
  PARAMETER_VECTOR( theta_z ); // SEM coefficients
  PARAMETER_VECTOR( nu_z ); // SEM coefficients
  PARAMETER_VECTOR( log_lambda ); //Penalization parameters
  PARAMETER_VECTOR( log_sigmaxi_l ); // SVC logSD

  PARAMETER_VECTOR( alpha2_j ); // Fixed covariate parameters
  PARAMETER_VECTOR( gamma2_k ); // Spline regression parameters
  PARAMETER_VECTOR( beta2_z ); // DSEM coefficients
  PARAMETER_VECTOR( theta2_z ); // SEM coefficients
  PARAMETER_VECTOR( nu2_z ); // SEM coefficients
  PARAMETER_VECTOR( log_lambda2 ); //Penalization parameters
  PARAMETER_VECTOR( log_sigmaxi2_l ); // delta_SVC logSD

  PARAMETER_VECTOR( log_sigma );
  PARAMETER_ARRAY( epsilon_stc );
  PARAMETER_ARRAY( omega_sc );
  PARAMETER_ARRAY( delta_tc );
  PARAMETER_ARRAY( xi_sl );

  PARAMETER_ARRAY( epsilon2_stc );
  PARAMETER_ARRAY( omega2_sc );
  PARAMETER_ARRAY( delta2_tc );
  PARAMETER_ARRAY( xi2_sl );

  PARAMETER_VECTOR( eps );     // manual epsilon bias-correction, empty to turn off

  // Globals
  Type nll = 0.0;

  // dimensions .. define minimally, to avoid later conflicts between object dimensions
  int n_s = epsilon_stc.dim(0);
  int n_i = y_i.size();
  int n_g = X_gj.rows();

  // Spatial distribution
  PARAMETER( log_kappa );
  PARAMETER_VECTOR( ln_H_input );

  // Anisotropy elements
  matrix<Type> H( ln_H_input.size(), ln_H_input.size() );
  H.setZero();
  H(0,0) = exp(ln_H_input(0));
  H(1,0) = ln_H_input(1);
  H(0,1) = ln_H_input(1);
  H(1,1) = (1.0 + ln_H_input(1)*ln_H_input(1)) / exp(ln_H_input(0));
  for( int i=2; i<ln_H_input.size(); i++ ){
    //H(i,i) = exp(ln_H_input(i));
    // Allowing negative values allows vertex_formula = 0 + x + I(x+y) to be the saem as vertex_formula = 0 + x + I(x-y)
    H(i,i) = ln_H_input(i);
  }

  REPORT( H );
  // Spatial settings
  Type log_tau = 0.0;
  Type range = NAN;
  Eigen::SparseMatrix<Type> Q_ss;
  // Using INLA with geometric anisotropy
  if( model_options(0)==1 ){
    DATA_STRUCT( spatial_list, R_inla::spde_aniso_t );
    // Build precision
    Q_ss = R_inla::Q_spde( spatial_list, exp(log_kappa), H );
    log_tau = log( 1.0 / (exp(log_kappa) * sqrt(4.0*M_PI)) );
    range = pow(8.0, 0.5) / exp( log_kappa );
    REPORT( range );
  }
  /// Using SAR
  if( model_options(0)==2 ){
    DATA_SPARSE_MATRIX( Adj );
    Eigen::SparseMatrix<Type> I_ss( Adj.rows(), Adj.rows() );
    Eigen::SparseMatrix<Type> Lspatial_ss( Adj.rows(), Adj.rows() );
    I_ss.setIdentity();
    Lspatial_ss = ( I_ss - exp(log_kappa)*Adj );
    log_tau = 0.0;
    Q_ss = Lspatial_ss.transpose() * Lspatial_ss;
  }
  // Off, but using INLA inputs
  if( model_options(0)==3 ){
    DATA_STRUCT(spatial_list, R_inla::spde_t);
    Q_ss = R_inla::Q_spde(spatial_list, Type(1.0));
    log_tau = Type(0.0);
  }
  // stream-network
  if( model_options(0)==4 ){
    DATA_IMATRIX( graph_sz );
    DATA_VECTOR( dist_s );
    Q_ss = Q_network2( log_kappa, n_s, graph_sz.col(0), graph_sz.col(1), dist_s );  // Q_network( log_theta, n_s, parent_s, child_s, dist_s )
    // DATA_SPARSE_MATRIX( Dist_ss );
    //Q_ss = Q_network( log_kappa, n_s, Dist_ss );
    log_tau = (log(2.0) + log_kappa) / 2;
    // diag(solve(Q)) = sigma^2 / (2 * theta)  ->  sigma^2 = 2*theta -> tau = 1/sqrt(2*theta)
    REPORT( Q_ss );
    //Eigen::SparseMatrix<Type> Q2_ss = Q_network( log_kappa, n_s, graph_sz.col(0), graph_sz.col(1), dist_s );  // Q_network( log_theta, n_s, parent_s, child_s, dist_s )
    //REPORT( Q2_ss );
  }
  // Using a SAR with geometric anisotropy
  if( model_options(0)==5 ){
    DATA_IVECTOR( i_z );
    DATA_IVECTOR( j_z );
    DATA_MATRIX( delta_z2 );
    log_tau = Type(0.0);
    // rho = inverse_cloglog( log_theta )
    Type rho = exp(-1.0 * exp(log_kappa));
    REPORT( rho );
    //Q_ss = Q_SAR( rho, H, n_s, i_z, j_z, delta_z2 );
    Q_ss = Q_SAR( exp(log_kappa), H, n_s, i_z, j_z, delta_z2 );
  }
  // Using SPDE with covariate-based anisotropy and geometric anisotropy
  if( model_options(0)==6 ){
    DATA_STRUCT( spatial_list, spde_covariates_t );
    //DATA_MATRIX( V_zk );
    PARAMETER_VECTOR( triangle_k );
    // Build precision
    //Eigen::SparseMatrix<Type> G1 = G_spde_covariates( spatial_list, H, V_zk, triangle_k );
    Eigen::SparseMatrix<Type> G1 = G_spde_covariates( spatial_list, H, triangle_k );
    REPORT( G1 );
    Q_ss = exp(4.0*log_kappa)*spatial_list.G0 + Type(2.0)*exp(2.0*log_kappa)*G1 + G1*spatial_list.G0_inv*G1;
    log_tau = log( 1.0 / (exp(log_kappa) * sqrt(4.0*M_PI)) );
    REPORT( Q_ss );
  }
  // Using NNGP
  if( model_options(0)==7 ){
    // DATA_INTEGER(n);      == n_s
    //DATA_IVECTOR(nn_index_flat);
    //DATA_IVECTOR(nn_start);
    //DATA_IVECTOR(nn_len);
    //DATA_VECTOR(dist_to_nn_flat);
    //DATA_VECTOR(dist_within_nn_flat);
    //DATA_IVECTOR(gp_order);
    //range = exp( log_kappa );
    //REPORT( range );
  }
  REPORT( log_tau );

  // spacetime_term
  Eigen::SparseMatrix<Type> Rho_hh = make_ram( ram_spacetime_term, ram_spacetime_term_start, beta_z, epsilon_stc.dim(1)*epsilon_stc.dim(2), int(0) );
  Eigen::SparseMatrix<Type> Gammainv_hh = make_ram( ram_spacetime_term, ram_spacetime_term_start, beta_z, epsilon_stc.dim(1)*epsilon_stc.dim(2), int(1) );
  Eigen::SparseMatrix<Type> Gamma_hh = make_ram( ram_spacetime_term, ram_spacetime_term_start, beta_z, epsilon_stc.dim(1)*epsilon_stc.dim(2), int(2) );

  // time_term
  Eigen::SparseMatrix<Type> Rho_time_hh = make_ram( ram_time_term, ram_time_term_start, nu_z, delta_tc.dim(0)*delta_tc.dim(1), int(0) );
  Eigen::SparseMatrix<Type> Gammainv_time_hh = make_ram( ram_time_term, ram_time_term_start, nu_z, delta_tc.dim(0)*delta_tc.dim(1), int(1) );
  Eigen::SparseMatrix<Type> Gamma_time_hh = make_ram( ram_time_term, ram_time_term_start, nu_z, delta_tc.dim(0)*delta_tc.dim(1), int(2) );

  // space_term
  Eigen::SparseMatrix<Type> Rho_cc = make_ram( ram_space_term, ram_space_term_start, theta_z, omega_sc.dim(1), int(0) );
  Eigen::SparseMatrix<Type> Gammainv_cc = make_ram( ram_space_term, ram_space_term_start, theta_z, omega_sc.dim(1), int(1) );
  Eigen::SparseMatrix<Type> Gamma_cc = make_ram( ram_space_term, ram_space_term_start, theta_z, omega_sc.dim(1), int(2) );

  // Delta spacetime_term
  Eigen::SparseMatrix<Type> Rho2_hh = make_ram( ram2_spacetime_term, ram2_spacetime_term_start, beta2_z, epsilon2_stc.dim(1)*epsilon2_stc.dim(2), int(0) );
  Eigen::SparseMatrix<Type> Gammainv2_hh = make_ram( ram2_spacetime_term, ram2_spacetime_term_start, beta2_z, epsilon2_stc.dim(1)*epsilon2_stc.dim(2), int(1) );
  Eigen::SparseMatrix<Type> Gamma2_hh = make_ram( ram2_spacetime_term, ram2_spacetime_term_start, beta2_z, epsilon2_stc.dim(1)*epsilon2_stc.dim(2), int(2) );

  // Delta time_term
  Eigen::SparseMatrix<Type> Rho2_time_hh = make_ram( ram2_time_term, ram2_time_term_start, nu2_z, delta2_tc.dim(0)*delta2_tc.dim(1), int(0) );
  Eigen::SparseMatrix<Type> Gammainv2_time_hh = make_ram( ram2_time_term, ram2_time_term_start, nu2_z, delta2_tc.dim(0)*delta2_tc.dim(1), int(1) );
  Eigen::SparseMatrix<Type> Gamma2_time_hh = make_ram( ram2_time_term, ram2_time_term_start, nu2_z, delta2_tc.dim(0)*delta2_tc.dim(1), int(2) );

  // Delta space_term
  Eigen::SparseMatrix<Type> Rho2_cc = make_ram( ram2_space_term, ram2_space_term_start, theta2_z, omega2_sc.dim(1), int(0) );
  Eigen::SparseMatrix<Type> Gammainv2_cc = make_ram( ram2_space_term, ram2_space_term_start, theta2_z, omega2_sc.dim(1), int(1) );
  Eigen::SparseMatrix<Type> Gamma2_cc = make_ram( ram2_space_term, ram2_space_term_start, theta2_z, omega2_sc.dim(1), int(2) );

  // space_term
  omega_sc = omega_distribution( omega_sc, model_options, Rho_cc,
                                   Gamma_cc, Gammainv_cc, Q_ss,
                                   exp(log_kappa), nngp_data, nll );          // nngp_data,
  omega2_sc = omega_distribution( omega2_sc, model_options, Rho2_cc,
                                   Gamma2_cc, Gammainv2_cc, Q_ss,
                                   exp(log_kappa), nngp_data, nll );           // nngp_data,

  // spacetime_term
  epsilon_stc = epsilon_distribution( epsilon_stc, model_options, Rho_hh,
                                     Gamma_hh, Gammainv_hh, Q_ss,
                                     exp(log_kappa), nngp_data, nll );
  epsilon2_stc = epsilon_distribution( epsilon2_stc, model_options, Rho2_hh,
                                       Gamma2_hh, Gammainv2_hh, Q_ss,
                                       exp(log_kappa), nngp_data, nll );

  // time_term
  delta_tc = delta_distribution( delta_tc, model_options, Rho_time_hh,
                                     Gamma_time_hh, Gammainv_time_hh, nll );
  delta2_tc = delta_distribution( delta2_tc, model_options, Rho2_time_hh,
                                       Gamma2_time_hh, Gammainv2_time_hh, nll );


  // Distribution for spline components
  nll += gamma_distribution( gamma_k, Sdims, Sblock, S_kk, log_lambda );
  nll += gamma_distribution( gamma2_k, S2dims, S2block, S2_kk, log_lambda2 );

  // Distribution for SVC components
  nll += xi_distribution( model_options, xi_sl, log_sigmaxi_l,
                          Q_ss, exp(log_kappa), nngp_data );
  nll += xi_distribution( model_options, xi2_sl, log_sigmaxi2_l,
                          Q_ss, exp(log_kappa), nngp_data );

  // Linear predictor .. keep partial effects for RSR adjustment below
  vector<Type> p_i( n_i );
  vector<Type> palpha1_i( n_i );
  vector<Type> pgamma1_i( n_i );
  palpha1_i.setZero();
  pgamma1_i.setZero();
  p_i = offset_i;
  if(alpha_j.size() > 0){ palpha1_i = X_ij*alpha_j; }
  if(gamma_k.size() > 0){ pgamma1_i = Z_ik*gamma_k; }
  vector<Type> pepsilon1_i = multiply_epsilon( Aepsilon_zz, Aepsilon_z, epsilon_stc, p_i.size() ) / exp(log_tau);
  vector<Type> pomega1_i = multiply_omega( Aomega_zz, Aomega_z, omega_sc, p_i.size() ) / exp(log_tau);
  vector<Type> pxi1_i = multiply_xi( A_is, xi_sl, W_il ) / exp(log_tau);
  vector<Type> pdelta1_i = multiply_delta( delta_tc, t_i, c_i, n_i );
  p_i += palpha1_i;
  p_i += pgamma1_i;
  p_i += pepsilon1_i;
  p_i += pomega1_i;
  p_i += pxi1_i;
  p_i += pdelta1_i;

  // 2nd linear predictor
  vector<Type> p2_i( n_i );
  vector<Type> palpha2_i( n_i );
  vector<Type> pgamma2_i( n_i );
  palpha2_i.setZero();
  pgamma2_i.setZero();
  p2_i.setZero();
  if(alpha2_j.size() > 0){ palpha2_i = X2_ij*alpha2_j; }
  if(gamma2_k.size() > 0){ pgamma2_i = Z2_ik*gamma2_k; }
  vector<Type> pepsilon2_i = multiply_epsilon( Aepsilon_zz, Aepsilon_z, epsilon2_stc, p_i.size() ) / exp(log_tau);
  vector<Type> pomega2_i = multiply_omega( Aomega_zz, Aomega_z, omega2_sc, p_i.size() ) / exp(log_tau);
  vector<Type> pxi2_i = multiply_xi( A_is, xi2_sl, W2_il ) / exp(log_tau);
  vector<Type> pdelta2_i = multiply_delta( delta2_tc, t_i, c_i, n_i );
  p2_i += palpha2_i;
  p2_i += pgamma2_i;
  p2_i += pepsilon2_i;
  p2_i += pomega2_i;
  p2_i += pxi2_i;
  p2_i += pdelta2_i;

  // Likelihood
  // relative_deviance != devresid^2 for hurdle model
  vector<Type> mu_i( n_i );
  vector<Type> devresid_i( n_i );
  vector<Type> negloglik_i( n_i );
  Type devresid = 0.0;
  Type deviance = 0.0;
  Type dev;
  Type nll_tmp;
  for( int i=0; i<n_i; i++ ) {       // PARALLEL_REGION
    vector<Type> log_sigma_segment = log_sigma.segment( Edims_ez(e_i(i),0), Edims_ez(e_i(i),1) );
    // Link function
    if( components_e(e_i(i))==1 ){
      mu_i(i) = one_predictor_likelihood( y_i(i), p_i(i), size_i(i), link_ez(e_i(i),0), family_ez(e_i(i),0), log_sigma_segment, nll_tmp, devresid, this );
      negloglik_i(i) = nll_tmp;
      nll += weights_i(i) * nll_tmp;
      deviance += pow( devresid, 2.0 );
      devresid_i(i) = devresid;
    }
    if( components_e(e_i(i))==2 ){
      mu_i(i) = two_predictor_likelihood( y_i(i), p_i(i), p2_i(i), size_i(i), link_ez.row(e_i(i)), family_ez.row(e_i(i)), log_sigma_segment, poislink_e(e_i(i)), nll_tmp, dev, this );
      negloglik_i(i) = nll_tmp;
      nll += weights_i(i) * nll_tmp;
      deviance += dev;
      devresid_i(i) = NAN;
    }
  }

  // Restricted spatial regression correction
  // t(A_is) %*% A_is = diag( 1s and 0s )
  // SO:  could adjust pomega1_i + pepsilon1_i + pxi1_i + pdelta1_i too
  // e.g. p_omegaprime1_i = pomega1_i + Z (X^T X)^-1 Z^T pomega1_i, where Z = A^T X
  if( model_options(2) == 1 ){
    matrix<Type> covX_jj = X_ij.transpose() * X_ij;
    matrix<Type> precisionX_jj = atomic::matinv(covX_jj);
    vector<Type> alphaprime_j = alpha_j + (precisionX_jj * (X_ij.transpose() * (pomega1_i + pepsilon1_i + pxi1_i + pdelta1_i).matrix())).array();
    matrix<Type> covX2_jj = X2_ij.transpose() * X2_ij;
    matrix<Type> precisionX2_jj = atomic::matinv(covX2_jj);
    vector<Type> alphaprime2_j = alpha2_j + (precisionX2_jj * (X2_ij.transpose() * (pomega2_i + pepsilon2_i + pxi2_i + pdelta2_i).matrix())).array();
    REPORT( alphaprime_j );
    REPORT( alphaprime2_j );
    ADREPORT( alphaprime_j );
    ADREPORT( alphaprime2_j );
  }

  // Predictions
  if( n_g > 0 ){
    vector<Type> palpha1_g = X_gj*alpha_j;
    vector<Type> pgamma1_g = Z_gk*gamma_k;
    vector<Type> pepsilon1_g = multiply_epsilon( AepsilonG_zz, AepsilonG_z, epsilon_stc, palpha1_g.size() ) / exp(log_tau);
    vector<Type> pomega1_g = multiply_omega( AomegaG_zz, AomegaG_z, omega_sc, palpha1_g.size() ) / exp(log_tau);
    vector<Type> pxi1_g = multiply_xi( A_gs, xi_sl, W_gl ) / exp(log_tau);
    vector<Type> pdelta1_g = multiply_delta( delta_tc, t_g, c_g, n_g );
    vector<Type> p1_g = palpha1_g + pgamma1_g+ pepsilon1_g + pomega1_g + pdelta1_g + pxi1_g + offset_g ;
    // Second linear predictor
    vector<Type> palpha2_g = X2_gj*alpha2_j;
    vector<Type> pgamma2_g = Z2_gk*gamma2_k;
    vector<Type> pepsilon2_g = multiply_epsilon( AepsilonG_zz, AepsilonG_z, epsilon2_stc, palpha2_g.size() ) / exp(log_tau);
    vector<Type> pomega2_g = multiply_omega( AomegaG_zz, AomegaG_z, omega2_sc, palpha2_g.size() ) / exp(log_tau);
    vector<Type> pxi2_g = multiply_xi( A_gs, xi2_sl, W2_gl ) / exp(log_tau);
    vector<Type> pdelta2_g = multiply_delta( delta2_tc, t_g, c_g, n_g );
    vector<Type> p2_g = palpha2_g + pgamma2_g + pepsilon2_g + pomega2_g + pdelta2_g + pxi2_g;
    // Combined
    vector<Type> p_g = p1_g + p2_g;
    vector<Type> mu_g( p_g.size() );
    for( int g=0; g<p1_g.size(); g++ ){
      switch( link_ez(e_g(g),0) ){
        case identity_link:
          mu_g(g) = p1_g(g);
          break;
        case log_link:
          mu_g(g) = exp(p1_g(g));
          break;
        case logit_link:
          mu_g(g) = invlogit(p1_g(g));
          break;
        case cloglog_link:
          mu_g(g) = Type(1.0) - exp( -1*exp(p1_g(g)) );
          break;
        default:
          error("Link not implemented.");
      }
      if( components_e(e_g(g))==2 ){
        // second link
        switch( link_ez(e_g(g),1) ){
          case identity_link:
            mu_g(g) *= p2_g(g);
            break;
          case log_link:
            mu_g(g) *= exp(p2_g(g));
            break;
          // case logit_link: // Logit
          default:
            error("Link not implemented.");
        }
      }
    }

    // Expansion
    if( (W_gz.rows()==mu_g.size()) && (V_gz.rows()==mu_g.size()) ){
      // First sweep
      vector<Type> phi0_g( mu_g.size() );
      for( int g=0; g<mu_g.size(); g++ ){
        if( (V_gz(g,0)==0) || (V_gz(g,0)==1) || (V_gz(g,0)==2) || (V_gz(g,0)==3) || (V_gz(g,0)==4) ){
          // Area-weighted average
          phi0_g(g) = mu_g(g) * W_gz(g,0);
        }
      }
      Type sumphi0 = sum(phi0_g);
      REPORT( phi0_g );

      // Second sweep for covariate or density-weighted averages
      int n_blocks = V_gz.col(2).maxCoeff();
      vector<Type> phi_g( mu_g.size() );
      vector<Type> sumphi_b( n_blocks );
      sumphi_b.setZero();
      phi_g.setZero();
      for( int g=0; g<mu_g.size(); g++ ){
        if( V_gz(g,0)==0 ){
          // Exclude from 2nd-sweep calculation
          phi_g(g) = 0.0;
        }
        if( V_gz(g,0)==1 ){
          // Default: equal to 1st sweep
          phi_g(g) = phi0_g(g);
        }
        if( V_gz(g,0)==2 ){
          // density-weighted average of W_gz(g,1)
          phi_g(g) = (phi0_g(g) / sumphi0) * W_gz(g,1);
        }
        if( (V_gz(g,0)==3) && (V_gz(g,1)>=0) && (V_gz(g,1)<=g) ){
          // density-weighted average of prediction
          phi_g(g) = (phi0_g(V_gz(g,1)) / sumphi0) * mu_g(g);
        }
        if( (V_gz(g,0)==4) && (V_gz(g,1)>=0) && (V_gz(g,1)<=g) ){
          // density-weighted total of prediction
          phi_g(g) = phi0_g(V_gz(g,1)) * mu_g(g);
        }

        // Combine them in blocks
        sumphi_b( V_gz(g,2) - 1 ) += phi_g(g);
      }

      //Type Metric = sum(phi_g);
      vector<Type> Metric( n_blocks );
      for( int b = 0; b < n_blocks; b++ ){
        Metric(b) = newton::Tag( sumphi_b(b) ); // Set lowrank tag on Metric = sum(exp(x))
      }
      // eps.size()=0 when apply_epsilon = FALSE
      if( eps.size() == Metric.size() ){
        nll += (Metric * eps).sum();
      }
      REPORT( phi_g );
      REPORT( Metric );
      ADREPORT( Metric );
    }

    REPORT( p1_g );
    REPORT( palpha1_g );
    REPORT( pgamma1_g );
    REPORT( pepsilon1_g );
    REPORT( pomega1_g );
    REPORT( pdelta1_g );
    REPORT( pxi1_g );
    REPORT( p2_g );
    REPORT( palpha2_g );
    REPORT( pgamma2_g );
    REPORT( pepsilon2_g );
    REPORT( pomega2_g );
    REPORT( pdelta2_g );
    REPORT( pxi2_g );
    REPORT( p_g );
    REPORT( mu_g );
    if(model_options(4) == 1) ADREPORT( p1_g );   // Keeping original order, in case someone reloads
    if(model_options(4) == 2) ADREPORT( p2_g );
    if(model_options(4) == 3) ADREPORT( p_g );
    if(model_options(4) == 4) ADREPORT( mu_g );
  }

  // Reporting
  REPORT( p_i );
  REPORT( p2_i );
  REPORT( mu_i );                      // Needed for `residuals.tinyVAST`
  REPORT( Q_ss );                      // Needed for forming precision of epsilon_sct in `project(.)`
  REPORT( devresid_i );
  REPORT( deviance );
  REPORT( nll );
  REPORT( negloglik_i );
  SIMULATE{
    REPORT(y_i);
  }

  //
  if( model_options(3) == 1 ){
    REPORT( Rho_hh );
    REPORT( Gamma_hh );
    REPORT( Gammainv_hh );
    REPORT( Rho_cc );
    REPORT( Gamma_cc );
    REPORT( Gammainv_cc );
    REPORT( Rho_time_hh );
    REPORT( Gamma_time_hh );
    REPORT( Gammainv_time_hh );
    REPORT( Rho2_hh );
    REPORT( Gamma2_hh );
    REPORT( Gammainv2_hh );
    REPORT( Rho2_cc );
    REPORT( Gamma2_cc );
    REPORT( Gammainv2_cc );
    REPORT( Rho2_time_hh );
    REPORT( Gamma2_time_hh );
    REPORT( Gammainv2_time_hh );
    REPORT( palpha1_i );
    REPORT( pgamma1_i );
    REPORT( pepsilon1_i );
    REPORT( pomega1_i );
    REPORT( pxi1_i );
    REPORT( pdelta1_i );
    REPORT( palpha2_i );
    REPORT( pgamma2_i );
    REPORT( pepsilon2_i );
    REPORT( pomega2_i );
    REPORT( pxi1_i );
    REPORT( pdelta2_i );
  }

  return nll;
}
