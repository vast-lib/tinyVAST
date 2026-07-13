namespace tinyVAST {

enum valid_family {
  gaussian_family  = 0,
  tweedie_family   = 1,
  lognormal_family = 2,
  poisson_family   = 3,
  bernoulli_family = 4,
  binomial_family = 4,
  gamma_family = 5,
  nbinom1_family = 6,
  nbinom2_family = 7,
  student_family = 8
};

enum valid_link {
  identity_link = 0,
  log_link      = 1,
  logit_link    = 2,
  cloglog_link = 3
};

// Old algorithmic constructor for tail-down exponential stream network, copied from VAST
//template<class Type>
//Eigen::SparseMatrix<Type> Q_network( Type log_theta,
//                                     int n_s,
//                                     vector<int> parent_s,
//                                     vector<int> child_s,
//                                     vector<Type> dist_s ){
//
//  Eigen::SparseMatrix<Type> Q( n_s, n_s );
//  Type theta = exp( log_theta );
//  for(int s=0; s<n_s; s++){
//    Q.coeffRef( s, s ) = Type(1.0);
//  }
//  for(int s=0; s<parent_s.size(); s++){
//  //for(int s=1; s<parent_s.size(); s++){
//    if( exp(-dist_s(s))!= Type(0.) ){
//      Type tmp = -exp(-theta*dist_s(s)) / (1-exp(-2*theta*dist_s(s)));
//      Type tmp2 = exp(-2*theta*dist_s(s)) / (1-exp(-2*theta*dist_s(s)));
//      Q.coeffRef( parent_s(s), child_s(s) ) = tmp;
//      Q.coeffRef( child_s(s), parent_s(s) ) = tmp;
//      Q.coeffRef( parent_s(s), parent_s(s) ) += tmp2;
//      Q.coeffRef( child_s(s), child_s(s) ) += tmp2;
//    }
//  }
//  return Q;
//}

template<class Type>
Eigen::SparseMatrix<Type> vectorsToSparseMatrix(
    const vector<int> &i,
    const vector<int> &j,
    const vector<Type> &x,
    int N ){

  // https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#title35
  using namespace Eigen;
  std::vector<Triplet<Type>> tripletList;
  for( int index = 0; index < i.size(); index ++ ){
    tripletList.push_back(Triplet<Type>( i[index], j[index], x[index]) );
  }
  SparseMatrix<Type> M( N, N );
  M.setFromTriplets( tripletList.begin(), tripletList.end() );

  //Eigen::SparseMatrix<Type> M( N, N );
  //for(int index=0; index<N; index++){
  //  M.coeffRef( i[index], j[index] ) += x[index];
  //}

  return M;
}

/** \brief Object containing all elements of an anisotropic SPDE object, i.e. eqn (20) in Lindgren et al. */
// Modified from: https://github.com/kaskr/adcomp/blob/master/TMB/inst/include/tmbutils/R_inla.hpp
template<class Type>
struct spde_covariates_t{
  int n_s;
  //int n_tri;
  vector<Type> Tri_Area;
  matrix<Type> E0;
  matrix<Type> E1;
  matrix<Type> E2;
  matrix<Type> V_zk;
  matrix<int>  TV;
  Eigen::SparseMatrix<Type> G0;
  Eigen::SparseMatrix<Type> G0_inv;

  // Define output
  spde_covariates_t(SEXP x){  /* x = List passed from R */
    n_s = CppAD::Integer(asVector<Type>(getListElement(x,"n_s"))[0]);
    //n_tri = CppAD::Integer(asVector<Type>(getListElement(x,"n_tri"))[0]);
    Tri_Area = asVector<Type>(getListElement(x,"Tri_Area"));
    E0 = asMatrix<Type>(getListElement(x,"E0"));
    E1 = asMatrix<Type>(getListElement(x,"E1"));
    E2 = asMatrix<Type>(getListElement(x,"E2"));
    V_zk = asMatrix<Type>(getListElement(x,"V_zk"));
    TV = asMatrix<int>(getListElement(x,"TV"));
    G0 = tmbutils::asSparseMatrix<Type>(getListElement(x,"G0"));
    G0_inv = tmbutils::asSparseMatrix<Type>(getListElement(x,"G0_inv"));
  }
};
// make stiffness G1 from covariates in columns of E0, E1, E2
template<class Type>
Eigen::SparseMatrix<Type> G_spde_covariates(
    const spde_covariates_t<Type> &spde,
    const matrix<Type> &H_jj,
    //matrix<Type> V_zk,
    const vector<Type> &triangle_k ){

  // TO ADD EDGE COVARIATES:
  // 0.  add n_e = edge0_tj.rows()
  // 1.  Add TE, where dim(TE) = n_t \times 3, listing the edge-number for each triangle-edge of n_e edges
  // 2.  Add A_ek design matrix, where dim(A_ek) = n_e \times n_k (edge-covariates)

  // Extract objects
  vector<Type> area_t = spde.Tri_Area;
  matrix<Type> edge0_tj = spde.E0;
  matrix<Type> edge1_tj = spde.E1;
  matrix<Type> edge2_tj = spde.E2;
  matrix<int> s_tv = spde.TV;         // dim = n_t \times n_v
  matrix<Type> triangle_t1 = spde.V_zk * triangle_k.matrix();
  int n_s = spde.n_s;
  int n_t = s_tv.rows();
  int n_v = s_tv.cols();              // number of vertices per triangle
  int n_j = edge0_tj.cols();          // number of covariates (default = 2)

  // Output
  Eigen::SparseMatrix<Type> G_ss(n_s, n_s);

  // Objects to assemble triangle contributions
  matrix<Type> edges_vj( n_v, n_j );
  matrix<Type> Gt_vv( n_v, n_v );

  // Calculate adjugate of H
  //matrix<Type> adjH_jj = H_jj.inverse() * H_jj.determinant();  // CAUSES CRASH
  matrix<Type> adjH_jj = H_jj.adjoint();

  // Assemble stiffness G
  for( int t=0; t<n_t; t++ ){
    edges_vj.row(0) = edge0_tj.row(t);
    edges_vj.row(1) = edge1_tj.row(t);
    edges_vj.row(2) = edge2_tj.row(t);

    // Make local stiffness ... strictly positive
    Gt_vv = (edges_vj * adjH_jj * edges_vj.transpose()) * exp(triangle_t1(t,0));

    // Make local stifness ... reverts to normal when triangle_t1 = 0, but not strictly positive
    //Gt_vv = (edges_vj * adjH_jj * edges_vj.transpose()) * (1.0 * triangle_t1(t,0));

    // Assemble
    Eigen::SparseMatrix<Type> Gt_ss( n_s, n_s );
    for(int v1 = 0; v1 < 3; v1++ ){
    for(int v2 = 0; v2 < 3; v2++ ){
      Gt_ss.coeffRef( s_tv(t,v1), s_tv(t,v2) ) = Gt_vv(v1,v2);
    }}
    G_ss += Gt_ss / ( 4.0 * area_t(t) );
  }
  return G_ss;
}

// Q_SAR( log_kappa, H, n_s, i_z, j_z, delta_z2 )
// Using a SAR (not standardized) to approximate spatial function
template<class Type>
Eigen::SparseMatrix<Type> Q_SAR(
    //Type rho,
    Type kappa,
    const matrix<Type> &H,
    int n_s,
    const vector<int> &i_z,
    const vector<int> &j_z,
    const matrix<Type> &delta_z2 ){

  //Type max_dist = delta_z2.cwiseAbs().maxCoeff();

  //d_z = sqrt((delta_z2 %*% H)^2 %*% matrix(1,nrow=2,ncol=1))[,1]
  matrix<Type> deltaprime_z2 = (delta_z2 * H); // / max_dist;
  vector<Type> dist_z( i_z.size() );
  for( int z = 0; z < i_z.size(); z++ ){
    dist_z(z) = pow( pow(deltaprime_z2(z,0),2.0) + pow(deltaprime_z2(z,1),2.0), 0.5 );
  }
  vector<Type> weight_z = exp( -kappa * dist_z );

  // Make weight matrix
  Eigen::SparseMatrix<Type> W_ss = vectorsToSparseMatrix(
    i_z, j_z, weight_z, n_s
  );

  // Row-standardize W_ss
  // Does *not* make sense when using anisotropy
  //for( int row = 0; row < n_s; row ++ ){
  //  Type rowsum = W_ss.row(row).sum();
  //  W_ss.row(row) /= rowsum;
  //}

  // Assemble
  Eigen::SparseMatrix<Type> I_ss( n_s, n_s );
  I_ss.setIdentity();
  Eigen::SparseMatrix<Type> IminusP_ss = I_ss - W_ss;
  Eigen::SparseMatrix<Type> Q = IminusP_ss.transpose() * IminusP_ss;
  return Q;
}

// New matrix-notation precision constructor for tail-down exponential stream network
template<class Type>
Eigen::SparseMatrix<Type> Q_network2(
    Type log_theta,
    int n_s,
    const vector<int> &parent_s,
    const vector<int> &child_s,
    const vector<Type> &dist_s ){

  // Compute vectors
  Type theta = exp(log_theta);
  vector<Type> v1_s = exp(-theta * dist_s) / (1.0 - exp(-2.0 * theta * dist_s));
  vector<Type> v2_s = exp(-2.0 * theta * dist_s) / (1.0 - exp(-2.0 * theta * dist_s));

  // Make components
  Eigen::SparseMatrix<Type> P_ss = vectorsToSparseMatrix( child_s,
                                                          parent_s,
                                                          v1_s,
                                                          n_s );
  Eigen::SparseMatrix<Type> D_ss = vectorsToSparseMatrix( child_s,
                                                          child_s,
                                                          v2_s,
                                                          n_s );
  Eigen::SparseMatrix<Type> invD_ss( n_s, n_s );
  for( int s=0; s<n_s; s++ ){
    D_ss.coeffRef(s,s) += 1.0;
    invD_ss.coeffRef(s,s) = 1.0 / D_ss.coeffRef(s,s);
  }

  // Assemble
  Eigen::SparseMatrix<Type> DminusP_ss = D_ss - P_ss;
  Eigen::SparseMatrix<Type> Q = DminusP_ss.transpose() * invD_ss * DminusP_ss;
  return Q;
}

// Functions to evaluate density for nearest-neighbors Gaussian process
template<class Type>
vector<Type> covariance_function(vector<Type> d, Type sigma2, Type inv_range) {
  return ((-d.array() * inv_range).exp() * sigma2).matrix();
}
template<class Type>
matrix<Type> covariance_function(matrix<Type> d, Type sigma2, Type inv_range) {
  return ((-d.array() * inv_range).exp() * sigma2).matrix();
}
/** \brief Object containing all data elements of an NNGP */
template<class Type>
struct nngp_data_t{
  vector<int> nn_index_flat;
  vector<int> nn_len;
  vector<Type> dist_to_nn_flat;
  vector<Type> dist_within_nn_flat;
  vector<int> gp_order;

  // Define output
  nngp_data_t(SEXP x){  /* x = List passed from R */
    nn_index_flat = asVector<int>(getListElement(x,"nn_index_flat"));
    nn_len = asVector<int>(getListElement(x,"nn_len"));
    dist_to_nn_flat = asVector<Type>(getListElement(x,"dist_to_nn_flat"));
    dist_within_nn_flat = asVector<Type>(getListElement(x,"dist_within_nn_flat"));
    gp_order = asVector<int>(getListElement(x,"gp_order"));
  }
};
// Density function .. sigma2 fixed at 1.0
template<class Type>
Type NNGP( Type sigma2,
           Type range,
           vector<Type> field_s,
           // Data
           nngp_data_t<Type> nngp_data ){

  // Using original order and ordered_structure, by calling gp_order(i) and gp_order(nn_ids), because:
  // 1.  using original version in new order is very slow! presumably the order is terrible for the inner Hessian
  // 2.  using re-ordered version of field is confusing for users
  Type nll = 0.0;
  Type inv_range = 1.0 / range; // Multiplication is faster than division
  vector<int> nn_index_flat = nngp_data.nn_index_flat;
  //vector<int> nn_start = nngp_data.nn_start;
  vector<int> nn_len = nngp_data.nn_len;
  vector<Type> dist_to_nn_flat = nngp_data.dist_to_nn_flat;
  vector<Type> dist_within_nn_flat = nngp_data.dist_within_nn_flat;
  vector<int> gp_order = nngp_data.gp_order;
  int n = field_s.size();

  int pos_mat = 0;
  int pos_vec = 0;
  for( int i = 0; i < n; i++ ){
    int k = nn_len(i);

    // No neighbors
    if(k == 0){
      nll -= dnorm( field_s(gp_order(i)), Type(0.0), sqrt(sigma2), true );
    }else{
      // Extract neighbor indices and distance to neighbors
      vector<int> nn_ids(k);
      vector<Type> dist_iN(k);
      for (int j = 0; j < k; j++) {
        nn_ids(j) = nn_index_flat(pos_vec);
        dist_iN(j) = dist_to_nn_flat(pos_vec);
        pos_vec++;
      }

      // Extract within-neighbor distance matrix (k x k) ...
      matrix<Type> dist_NN(k, k);
      for (int r = 0; r < k; r++) {
        for (int c = 0; c < k; c++) {
          dist_NN(r, c) = dist_within_nn_flat(pos_mat);
          pos_mat++;
        }
      }

      // Covariances
      matrix<Type> Sigma_NN = covariance_function(dist_NN, sigma2, inv_range);
      vector<Type> Sigma_iN = covariance_function(dist_iN, sigma2, inv_range);

      // Solve
      //vector<Type> a_i = solve(Sigma_NN, Sigma_iN);
      //vector<Type> a_i = Sigma_NN.ldlt().solve( Sigma_iN.matrix() );
      //Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > invSigma_NN;
      //invSigma_NN.compute(Sigma_NN);
      //vector<Type> a_i = invSigma_NN.solve(Sigma_iN);
      matrix<Type> Sigma_NN_inv = atomic::matinv(Sigma_NN);
      vector<Type> a_i = Sigma_NN_inv * Sigma_iN;

      // Residual variance
      Type resid_var = sigma2 - (a_i * Sigma_iN).sum(); // + 1e-12;

      // Conditional mean
      Type cond_mean = 0.0;
      for (int j = 0; j < k; j++) {
        cond_mean += a_i(j) * field_s(gp_order(nn_ids(j)));
      }

      // Likelihood
      nll -= dnorm( field_s(gp_order(i)), cond_mean, sqrt(resid_var), true );
    }
  }
  return nll;
}

// Function to calculate RAM matrices
template<class Type>
Eigen::SparseMatrix<Type> make_ram( matrix<int> ram,
                                    vector<Type> ram_start,
                                    vector<Type> beta_z,
                                    int n_c,
                                    int what ){

  Eigen::SparseMatrix<Type> out_cc(n_c, n_c);
  out_cc.setZero();
  Type tmp;

  for(int r=0; r<ram.rows(); r++){
    // Extract estimated or fixed value
    if(ram(r,3)>=1){
      tmp = beta_z(ram(r,3)-1);
    }else{
      tmp = ram_start(r);
    }
    // Rho_cc
    if( (ram(r,0)==1) && (what==0) ){
       out_cc.coeffRef( ram(r,1)-1, ram(r,2)-1 ) = tmp;
    }
    // Gammainv_cc
    if( (ram(r,0)==2) && (what==1) ){
      out_cc.coeffRef( ram(r,1)-1, ram(r,2)-1 ) = 1 / tmp;
    }
    // Gamma_cc
    if( (ram(r,0)==2) && (what==2) ){
      out_cc.coeffRef( ram(r,1)-1, ram(r,2)-1 ) = tmp;
    }
  }

  return out_cc;
}

// distribution/projection for gamma
// Based on:
//   https://github.com/meganferg/FergusonEtal_20250125_EBS_Beluga_DSM/blob/main/src/te_tw_DSM.cpp#L85C9-L85C40 (Megan Ferguson)
template<class Type>
Type gamma_distribution( vector<Type> gamma_k,
                         vector<int> Sdims,
                         vector<int> Sblock,
                         Eigen::SparseMatrix<Type> S_kk,
                         vector<Type> log_lambda ){

  using namespace density;
  Type nll = 0.0;
  int start_gamma = 0;   // Counter for gamma
  int start_k = 0;       // counter for columns of S_kk
  int index_lambda = 0;
  for(int z=0; z<Sdims.size(); z++){
    int ngamma_z = Sdims(z);
    // Recover gamma_segment
    vector<Type> gamma_segment = gamma_k.segment(start_gamma, ngamma_z);
    start_gamma += ngamma_z;

    // s() splines
    if( Sblock(z) == 1 ){
      // Recover S_block
      Eigen::SparseMatrix<Type> S_block = S_kk.block(start_k, start_k, ngamma_z, ngamma_z);  // Recover S_i

      /// Quadratic-form GRMF
      // GMRF(Q).Quadform(x) caused a problem calculating the log-determinant unnecessarily, resulting in valgrind error
      //nll -= Type(0.5)*Type(m_z)*log_lambda(z) - Type(0.5)*exp(log_lambda(z))*GMRF(S_block).Quadform(gamma_segment);
      nll -= Type(0.5)*Type(ngamma_z)*log_lambda(index_lambda) - Type(0.5) * exp(log_lambda(index_lambda)) *
                                                    (gamma_segment.matrix().transpose()*(S_block * gamma_segment.matrix())).sum();

      start_k += ngamma_z;
      index_lambda += 1;
    }

    // te() / ti() / t2() splines
    if( Sblock(z) >= 2 ){
      // Recover S_block ... Eigen::SparseMatrix<Type> appears to initialize as empty so no .setZero() equivalent needed
      Eigen::SparseMatrix<Type> S_block( ngamma_z, ngamma_z );
      for( int block=0; block < Sblock(z); block ++ ){
        Eigen::SparseMatrix<Type> S_tmp = S_kk.block(start_k, start_k, ngamma_z, ngamma_z);  // Recover S_i
        S_block += exp(log_lambda(index_lambda)) * S_tmp;
        start_k += ngamma_z;
        index_lambda++;
      }

      /// Quadratic-form GRMF
      nll -= 0.5*(atomic::logdet(matrix<Type>(S_block))) - 0.5 *
                                                    (gamma_segment.matrix().transpose()*(S_block * gamma_segment.matrix())).sum();
    }
  }
  return( nll );
}

// distribution/projection for omega
template<class Type>
tmbutils::array<Type> omega_distribution( tmbutils::array<Type> omega_sc,
                                 vector<int> model_options,
                                 Eigen::SparseMatrix<Type> Rho_cc,
                                 Eigen::SparseMatrix<Type> Gamma_cc,
                                 Eigen::SparseMatrix<Type> Gammainv_cc,
                                 Eigen::SparseMatrix<Type> Q_ss,
                                 Type range,
                                 nngp_data_t<Type> nngp_data,
                                 Type &nll ){

  if( omega_sc.size() > 0 ){
    int n_c = omega_sc.dim(1);
    using namespace density;
    Eigen::SparseMatrix<Type> I_cc( n_c, n_c );
    I_cc.setIdentity();
    if( omega_sc.size()>0 ){ // PARALLEL_REGION
      Eigen::SparseMatrix<Type> IminusRho_cc = I_cc - Rho_cc;
      if( model_options(1) == 0 ){
        // Separable precision ... Option-1
        //Eigen::SparseMatrix<Type> Linv_cc = Gammainv_cc * ( I_cc - Rho_cc );
        //Eigen::SparseMatrix<Type> Q_cc = Linv_cc.transpose() * Linv_cc;

        // Separable precision ... Option-2
        // Only compute Vinv_kk if Gamma_kk is full rank
        Eigen::SparseMatrix<Type> V_cc = Gamma_cc.transpose() * Gamma_cc;
        matrix<Type> Vinv_cc = invertSparseMatrix( V_cc );
        Eigen::SparseMatrix<Type> Vinv2_cc = asSparseMatrix( Vinv_cc );
        Eigen::SparseMatrix<Type> Q_cc = IminusRho_cc.transpose() * Vinv2_cc * IminusRho_cc;

        // GMRF for SEM:  separable variable-space
        nll += SEPARABLE( GMRF(Q_cc), GMRF(Q_ss) )( omega_sc );
        // Including this line with Makevars below seems to cause a crash:
        // PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS)
        // PKG_CXXFLAGS=$(SHLIB_OPENMP_CXXFLAGS)
      }else{
        // Rank-deficient (projection) method
        if( model_options(0) == 7 ){
          vector<Type> omega_s( omega_sc.rows() );
          for( int ci = 0; ci < omega_sc.cols(); ci++ ){
            omega_s = omega_sc.col(ci);
            nll += NNGP( Type(1.0), range, omega_s, nngp_data );
          }
        }else{
          nll += SEPARABLE( GMRF(I_cc), GMRF(Q_ss) )( omega_sc );
        }

        // Sparse inverse-product
        //Eigen::SparseMatrix<Type> IminusRho_cc = I_cc - Rho_cc;
        Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > inverseIminusRho_cc;
        inverseIminusRho_cc.compute(IminusRho_cc);

        // (I-Rho)^{-1} * Gamma * Epsilon
        matrix<Type> omega2_cs = Gamma_cc * omega_sc.matrix().transpose();
        matrix<Type> omega3_cs = inverseIminusRho_cc.solve(omega2_cs);
        omega_sc = omega3_cs.transpose();
      }
    }
  }

  return omega_sc;
}

// distribution/projection for xi
template<class Type>
Type xi_distribution( vector<int> model_options,
                      tmbutils::array<Type> xi_sl,
                      vector<Type> log_sigmaxi_l,
                      Eigen::SparseMatrix<Type> Q_ss,
                      Type range,
                      nngp_data_t<Type> nngp_data ){

  Type nll = 0;
  if( xi_sl.size() > 0 ){
    int n_l = xi_sl.cols();
    using namespace density;
    if( model_options(0) == 7 ){
      vector<Type> xi_s( xi_sl.rows() );
      for( int l=0; l<n_l; l++ ){
        xi_s = xi_sl.col(l);
        NNGP( exp(2.0*log_sigmaxi_l(l)), range, xi_s, nngp_data );
      }
    }else{
      for( int l=0; l<n_l; l++ ){
        nll += SCALE( GMRF(Q_ss), exp(log_sigmaxi_l(l)) )( xi_sl.col(l) );
      }
    }
  }

  return nll;
}

// distribution/projection for epsilon
// Passing const with pointer to simplify compilation (and avoid potential compiler errors during devtools checks)
template<class Type>
tmbutils::array<Type> epsilon_distribution(
    tmbutils::array<Type> epsilon_stc,
    const vector<int> &model_options,
    const Eigen::SparseMatrix<Type> &Rho_hh,
    const Eigen::SparseMatrix<Type> &Gamma_hh,
    const Eigen::SparseMatrix<Type> &Gammainv_hh,
    const Eigen::SparseMatrix<Type> &Q_ss,
    Type range,
    const nngp_data_t<Type> &nngp_data,
    Type &nll ){

  if( epsilon_stc.size() > 0 ){
    int n_s = epsilon_stc.dim(0);
    int n_t = epsilon_stc.dim(1);
    int n_c = epsilon_stc.dim(2);
    int n_h = n_t * n_c;
    using namespace density;
    Eigen::SparseMatrix<Type> I_hh( n_h, n_h );
    I_hh.setIdentity();
    int h;

    //if( epsilon_stc.size()>0 ){ // PARALLEL_REGION
    // Reshape for either model_options
    tmbutils::array<Type> epsilon_hs( n_h, n_s );
    for( int s=0; s<n_s; s++ ){
    for( int t=0; t<n_t; t++ ){
    for( int c=0; c<n_c; c++ ){
      h = c*n_t + t;
      epsilon_hs(h,s) = epsilon_stc(s,t,c);
    }}}
    Eigen::SparseMatrix<Type> IminusRho_hh = I_hh - Rho_hh;

    if( model_options(1) == 0 ){
      // Separable precision ... Option-1
      //Eigen::SparseMatrix<Type> Linv_hh = Gammainv_hh * ( I_hh - Rho_hh );
      //Eigen::SparseMatrix<Type> Q_hh = Linv_hh.transpose() * Linv_hh;

      // Separable precision ... Option-2
      // Only compute Vinv_kk if Gamma_kk is full rank
      Eigen::SparseMatrix<Type> V_hh = Gamma_hh.transpose() * Gamma_hh;
      matrix<Type> Vinv_hh = invertSparseMatrix( V_hh );
      Eigen::SparseMatrix<Type> Vinv2_hh = asSparseMatrix( Vinv_hh );
      Eigen::SparseMatrix<Type> Q_hh = IminusRho_hh.transpose() * Vinv2_hh * IminusRho_hh;

      // GMRF for DSEM:  non-separable time-variable, with separable space
      nll += SEPARABLE( GMRF(Q_ss), GMRF(Q_hh) )( epsilon_hs );
    }else{
      // Rank-deficient (projection) method
      if( model_options(0) == 7 ){
        vector<Type> epsilon_s( epsilon_hs.cols() );
        for( int hi = 0; hi < epsilon_hs.rows(); hi++ ){
          epsilon_s = epsilon_hs.matrix().row(hi);
          nll += NNGP( Type(1.0), range, epsilon_s, nngp_data );
        }
      }else{
        nll += SEPARABLE( GMRF(Q_ss), GMRF(I_hh) )( epsilon_hs );
      }

      // Sparse inverse-product
      //Eigen::SparseMatrix<Type> IminusRho_hh = I_hh - Rho_hh;
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
    }
    //}
  }

  return epsilon_stc;
}

// distribution/projection for epsilon
template<class Type>
tmbutils::array<Type> delta_distribution(
    tmbutils::array<Type> delta_tc,
    const vector<int> &model_options,
    const Eigen::SparseMatrix<Type> &Rho_hh,
    const Eigen::SparseMatrix<Type> &Gamma_hh,
    const Eigen::SparseMatrix<Type> &Gammainv_hh,
    Type &nll ){

  if( delta_tc.size() > 0 ){
    int n_t = delta_tc.dim(0);
    int n_c = delta_tc.dim(1);
    int n_h = n_t * n_c;
    using namespace density;
    Eigen::SparseMatrix<Type> I_hh( n_h, n_h );
    I_hh.setIdentity();
    int h;

    // Reshape for either model_options
    tmbutils::array<Type> delta_h1( n_h, 1 );
    for( int t=0; t<n_t; t++ ){
    for( int c=0; c<n_c; c++ ){
      h = c*n_t + t;
      delta_h1(h,0) = delta_tc(t,c);
    }}
    Eigen::SparseMatrix<Type> IminusRho_hh = I_hh - Rho_hh;

    if( model_options(1) == 0 ){
      Eigen::SparseMatrix<Type> V_hh = Gamma_hh.transpose() * Gamma_hh;
      matrix<Type> Vinv_hh = invertSparseMatrix( V_hh );
      Eigen::SparseMatrix<Type> Vinv2_hh = asSparseMatrix( Vinv_hh );
      Eigen::SparseMatrix<Type> Q_hh = IminusRho_hh.transpose() * Vinv2_hh * IminusRho_hh;

      // GMRF for DSEM:  non-separable time-variable, with separable space
      nll += GMRF(Q_hh)( delta_h1 );
    }else{
      // Rank-deficient (projection) method
      nll += GMRF(I_hh)( delta_h1 );

      // Sparse inverse-product
      //Eigen::SparseMatrix<Type> IminusRho_hh = I_hh - Rho_hh;
      Eigen::SparseLU< Eigen::SparseMatrix<Type>, Eigen::COLAMDOrdering<int> > inverseIminusRho_hh;
      inverseIminusRho_hh.compute(IminusRho_hh);

      // (I-Rho)^{-1} * Gamma * Epsilon
      matrix<Type> e2_h1 = Gamma_hh * delta_h1.matrix();
      matrix<Type> e3_h1 = inverseIminusRho_hh.solve(e2_h1);

      // Transformations
      for( int t=0; t<n_t; t++ ){
      for( int c=0; c<n_c; c++ ){
        h = c*n_t + t;
        delta_tc(t,c) = e3_h1(h,0);
      }}
    }
  }

  return delta_tc;
}

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Sparse array * matrix
// NAs in IVECTOR or IMATRIX get converted to -2147483648, so isNA doesn't work ... instead drop NAs from A prior to passing to TMB
template<class Type>
vector<Type> multiply_epsilon(
    const matrix<int> &A,
    const vector<Type> &weight,
    const tmbutils::array<Type> &x,
    int n_i ){

  vector<Type> out( n_i );
  out.setZero();
  if( x.size() > 0 ){
    for( int z=0; z<A.rows(); z++ ){
      out(A(z,0)) += weight(z) * x(A(z,1),A(z,2),A(z,3));
    }
  }
  return out;
}
template<class Type>
vector<Type> multiply_omega(
    const matrix<int> &A,
    const vector<Type> &weight,
    const tmbutils::array<Type> &x,
    int n_i ){

  vector<Type> out( n_i );
  out.setZero();
  if( x.size() > 0 ){
    for( int z=0; z<A.rows(); z++ ){
      out(A(z,0)) += weight(z) * x(A(z,1),A(z,2));
    }
  }
  return out;
}
template<class Type>
vector<Type> multiply_xi(
    const Eigen::SparseMatrix<Type> &A_is,
    const tmbutils::array<Type> &xi_sl,
    const matrix<Type> &W_il ){

  vector<Type> out( W_il.rows() );
  out.setZero();
  if( xi_sl.size() > 0 ){
    matrix<Type> xi_il = A_is * xi_sl.matrix();
    for( int i=0; i<xi_il.rows(); i++ ){
    for( int l=0; l<xi_il.cols(); l++ ){
      out(i) += xi_il(i,l) * W_il(i,l);
    }}
  }
  return out;
}
template<class Type>
vector<Type> multiply_delta(
    const tmbutils::array<Type> &delta_tc,
    const vector<int> &t_i,
    const vector<int> &c_i,
    int n_i ){

  vector<Type> delta_i( n_i );
  delta_i.setZero();
  if( delta_tc.size() > 0 ){
    for( int i=0; i<n_i; i++ ){
      delta_i(i) += delta_tc( t_i(i), c_i(i) );
    }
  }
  return delta_i;
}

// get sign of double, only for REPORT use
template<class Type>
Type sign(Type x){
  return x / pow(pow(x,2),0.5);
}

// dlnorm
template<class Type>
Type dlnorm(
    Type x,
    Type meanlog,
    Type sdlog,
    int give_log=0){

  //return 1/(sqrt(2*M_PI)*sd) * exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

// dlnorm
template<class Type>
Type dbinom_custom(
    Type x,
    Type log_prob,
    Type log_one_minus_prob,
    Type size,
    int give_log = 0 ){

  // size choose x
  Type logres = lgamma(size + 1.0) - lgamma(x + 1.0) - lgamma(size - x + 1.0);
  // p^x * (1-p)^(size-x)
  logres += x * log_prob  +  (size-x) * log_one_minus_prob;
  if(give_log) return logres; else return exp(logres);
}

// dlnorm
template<class Type>
Type devresid_binom(
    Type y,
    Type weight,  // weight = size
    Type mu ){

  Type y_weight = y * weight;
  Type mu_weight = mu * weight;
  Type p1, p2;
  if(y_weight == 0){
    p1 = 0.0;
  }else{
    p1 = y_weight * log(y_weight / mu_weight);
  }
  if( (weight - y_weight) == 0 ){
    p2 = 0.0;
  }else{
    p2 = (weight - y_weight) * log( (weight - y_weight) / (weight - mu_weight) );
  }
  Type devresid = sign(y - mu) * pow(2 * (p1 + p2), 0.5);
  return devresid;
}

// Deviance for the Tweedie
// https://en.wikipedia.org/wiki/Tweedie_distribution#Properties
template<class Type>
Type devresid_tweedie(
    Type y,
    Type mu,
    Type p ){

  Type c1 = pow( y, 2.0-p ) / (1.0-p) / (2.0-p);
  Type c2 = y * pow( mu, 1.0-p ) / (1.0-p);
  Type c3 = pow( mu, 2.0-p ) / (2.0-p);
  Type deviance = 2 * (c1 - c2 + c3 );
  Type devresid = sign( y - mu ) * pow( deviance, 0.5 );
  return devresid;
}

// Deviance for the student-T
// from chatGPT (experimental) ... converges on Gaussian as df -> Inf
template<class Type>
Type devresid_student(
    Type y,
    Type mu,
    Type sigma,
    Type df ){

  Type deviance = (df + 1.0) * log( 1.0 + pow(y-mu,2.0) / (df * sigma*sigma) );
  Type devresid = sign( y - mu ) * pow( deviance, 0.5 );
  return devresid;
}

// COPIED from sdmTMB.cpp in sdmTMB package on Nov. 14, 2025
template <class Type>
Type dstudent(
    Type x,
    Type mean,
    Type sigma,
    Type df,
    int give_log = 0) {

  // from metRology::dt.scaled()
  // dt((x - mean)/sd, df, ncp = ncp, log = TRUE) - log(sd)
  Type logres = dt((x - mean) / sigma, df, true) - log(sigma);
  if (give_log)
    return logres;
  else
    return exp(logres);
}

// Deviance for the Negative binomial
//template<class Type>
//Type devresid_nbinom1( Type y,
//                       Type logmu,
//                       Type logtheta ){
//
//  // var - mu = exp( log(mu) + log(theta) ) = theta * mu  -->  var = (theta+1) * mu
//  Type logp1 = dnbinom_robust( y, log(y + Type(1e-10)), log(y + Type(1e-10)) + logtheta, true );
//  Type logp2 = dnbinom_robust( y, logmu, logmu + logtheta, true );
//  Type deviance = 2 * (logp1 - logp2);
//  Type devresid = sign( y - exp(logmu) ) * pow( deviance, 0.5 );
//  return devresid;
//}
template<class Type>
Type devresid_nbinom2(
    Type y,
    Type logmu,
    Type logtheta ){

  // var - mu = exp( 2 * log(mu) - log(theta) ) = mu^2 / theta  -->  var = mu + mu^2 / theta
  Type logp1 = dnbinom_robust( y, log(y + Type(1e-10)), Type(2.0) * log(y + Type(1e-10)) - logtheta, true );
  Type logp2 = dnbinom_robust( y, logmu, Type(2.0) * logmu - logtheta, true );
  Type deviance = 2 * (logp1 - logp2);
  Type devresid = sign( y - exp(logmu) ) * pow( deviance, 0.5 );
  return devresid;
}


// distribution/projection for epsilon
// deviance = deviance1 + deviance2
template<class Type>
Type one_predictor_likelihood(
    Type &y,
    Type p,
    Type size,
    int link,
    int family,
    vector<Type> log_sigma_segment,
    Type &nll,
    Type &devresid,
    objective_function<Type>* of ){

  Type mu, logmu, log_one_minus_mu;
  //Type devresid = 0;

  switch( link ){
    case identity_link:
      mu = p;
      logmu = log( p );
      log_one_minus_mu = log( Type(1.0) - mu );
      break;
    case log_link:
      mu = exp(p);
      logmu = p;
      log_one_minus_mu = log( Type(1.0) - mu );
      break;
    case logit_link:
      mu = invlogit(p);
      //logmu = log(mu);
      logmu = p - logspace_add( Type(0.0), p );
      //log_one_minus_mu = log( invlogit(-1 * p) );
      log_one_minus_mu = 0.0 - logspace_add( Type(0.0), p );
      break;
    case cloglog_link:
      mu = Type(1.0) - exp( -1 * exp(p) );
      logmu = logspace_sub( Type(0.0), -1.0 * exp(p) );
      log_one_minus_mu = -1.0 * exp(p);
      break;
    default:
      error("Link not implemented.");
  }
  if( !R_IsNA(asDouble(y)) ){
    // Distribution
    switch( family ){
      case gaussian_family:
        nll = -1 * dnorm( y, mu, exp(log_sigma_segment(0)), true );
        devresid = y - mu;
        if(isDouble<Type>::value && of->do_simulate){
          y = rnorm( mu, exp(log_sigma_segment(0)) );
        }
        break;
      case tweedie_family:
        nll = -1 * dtweedie( y, mu, exp(log_sigma_segment(0)), 1.0 + invlogit(log_sigma_segment(1)), true );
        devresid = devresid_tweedie( y, mu, 1.0 + invlogit(log_sigma_segment(1)) );
        if(isDouble<Type>::value && of->do_simulate){
          y = rtweedie( mu, exp(log_sigma_segment(0)), 1.0 + invlogit(log_sigma_segment(1)) );
        }
        break;
      case lognormal_family:
        nll = -1 * dlnorm( y, logmu - 0.5*exp(2.0*log_sigma_segment(0)), exp(log_sigma_segment(0)), true );
        devresid = log(y) - ( logmu - 0.5*exp(2.0*log_sigma_segment(0)) );
        if(isDouble<Type>::value && of->do_simulate){
          y = exp(rnorm( logmu - 0.5*exp(2.0*log_sigma_segment(0)), exp(log_sigma_segment(0)) ));
        }
        break;
      case poisson_family:
        nll = -1 * dpois( y, mu, true );
        devresid = sign(y - mu) * pow(2*(y*log((Type(1e-10) + y)/mu) - (y-mu)), 0.5);
        if(isDouble<Type>::value && of->do_simulate){
          y = rpois( mu );
        }
        break;
      case binomial_family:
        //if(y==0){
        //  nll = -1 * log_one_minus_mu;
        //}else{
        //  nll = -1 * logmu;
        //}
        nll = -1 * dbinom_custom( y * size, logmu, log_one_minus_mu, size, true );
        if(isDouble<Type>::value && of->do_simulate){
          y = rbinom( size, mu );
        }
        // TODO:  Update deviance residual for Trials = size
        //devresid = sign(y - mu) * pow(-2*((1-y)*log(1.0-mu) + y*log(mu)), 0.5);
        devresid = devresid_binom( y, size, mu );
        break;
      case gamma_family: // shape = 1/CV^2;   scale = mean*CV^2
        nll = -1 * dgamma( y, exp(-2.0*log_sigma_segment(0)), mu*exp(2.0*log_sigma_segment(0)), true );
        devresid = sign(y - mu) * pow(2 * ( (y-mu)/mu - log(y/mu) ), 0.5);
        if(isDouble<Type>::value && of->do_simulate){
          y = rgamma( exp(-2.0*log_sigma_segment(0)), mu*exp(2.0*log_sigma_segment(0)) );
        }
        break;
      case nbinom1_family:   // dnbinom_robust( x, log(mu_i), log(var - mu) )
        // var - mu = exp( log(mu) + log(theta) ) = theta * mu  -->  var = (theta+1) * mu
        nll = -1 * dnbinom_robust( y, logmu, logmu + log_sigma_segment(0), true);
        devresid = devresid_nbinom2( y, logmu, logmu - log_sigma_segment(0) );    // theta = mu / phi
        if(isDouble<Type>::value && of->do_simulate){
          // rnbinom2( mu, var )
          y = rnbinom2( mu, mu * (Type(1.0) + exp(log_sigma_segment(0))) );
        }
        break;
      case nbinom2_family:  // dnbinom_robust( x, log(mu_i), log(var - mu) )
        // var - mu = exp( 2 * log(mu) - log(theta) ) = mu^2 / theta  -->  var = mu + mu^2 / theta
        nll = -1 * dnbinom_robust( y, logmu, Type(2.0) * logmu - log_sigma_segment(0), true);
        devresid = devresid_nbinom2( y, logmu, log_sigma_segment(0) );
        if(isDouble<Type>::value && of->do_simulate){
          // rnbinom2( mu, var )
          y = rnbinom2( mu, mu * (Type(1.0) + mu / exp(log_sigma_segment(0))) );
        }
        break;
      case student_family:  // dnbinom_robust( x, log(mu_i), log(var - mu) )
        // var - mu = exp( 2 * log(mu) - log(theta) ) = mu^2 / theta  -->  var = mu + mu^2 / theta
        nll = -1 * dstudent( y, mu, exp(log_sigma_segment(0)), 1.0 + exp(log_sigma_segment(1)), true);
        devresid = devresid_student( y, mu, exp(log_sigma_segment(0)), 1.0 + exp(log_sigma_segment(1)) );
        if(isDouble<Type>::value && of->do_simulate){
          y = mu + exp(log_sigma_segment(0)) * rt(1.0 + exp(log_sigma_segment(1)));
        }
        break;
      default:
        error("Distribution not implemented.");
    }
  }

  return mu;
}

// distribution/projection for epsilon
template<class Type>
Type two_predictor_likelihood(
    Type &y,
    Type p1,
    Type p2,
    Type size,
    vector<int> link,
    vector<int> family,
    vector<Type> log_sigma_segment,
    int poislink,
    Type &nll,
    Type &dev,
    objective_function<Type>* of ){

  Type mu1, logmu1, mu2, logmu2, log_one_minus_mu1;

  // First link
  if( poislink==0 ){
    mu1 = invlogit( p1 );
    logmu1 = log( mu1 );
    log_one_minus_mu1 = log( Type(1.0) - mu1 );
    // second link
    switch( link(1) ){
      case identity_link:
        mu2 = p2;
        logmu2 = log( p2 );
        break;
      case log_link:
        mu2 = exp(p2);
        logmu2 = p2;
        break;
      default:
        error("Link not implemented.");
    }
  }else{
    mu1 = Type(1.0) - exp( -1.0*exp(p1) );
    logmu1 = logspace_sub( Type(0.0), -1*exp(p1) );
    log_one_minus_mu1 = -1.0*exp(p1);
    mu2 = exp( p1 + p2 ) / mu1;
    logmu2 = p1 + p2 - logmu1;
  }
  if( !R_IsNA(asDouble(y)) ){
    // Distribution
    if(isDouble<Type>::value && of->do_simulate){
      y = rbinom( Type(1), mu1 );
    }
    if( y == 0 ){
      nll = -1 * log_one_minus_mu1;
      dev = -2 * log_one_minus_mu1;
    }
    if( y>0 ){  // Not if-else so y>0 triggered when simulating y>0
      nll = -1 * logmu1;
      dev = -2 * logmu1;
      //deviance1_i(i) = -2 * log_mu1(i);
      switch( family(1) ){
        case gaussian_family:
          nll -= dnorm( y, mu2, exp(log_sigma_segment(0)), true );
          dev += pow(y - mu2, 2.0);
          if(isDouble<Type>::value && of->do_simulate){
            y = rnorm( mu2, exp(log_sigma_segment(0)) );
          }
          break;
        //case tweedie_family:
        //  nll -= dtweedie( y, mu2, exp(log_sigma_segment(0)), 1.0 + invlogit(log_sigma_segment(1)), true );
        //  devresid += devresid_tweedie( y, mu, 1.0 + invlogit(log_sigma_segment(1)) );
        //  if(isDouble<Type>::value && of->do_simulate){
        //    y = rtweedie( mu2, exp(log_sigma_segment(0)), 1.0 + invlogit(log_sigma_segment(1)) );
        //  }
        //  break;
        case lognormal_family:
          nll -= dlnorm( y, logmu2 - 0.5*exp(2.0*log_sigma_segment(0)), exp(log_sigma_segment(0)), true );
          dev += pow( log(y) - (logmu2 - 0.5*exp(2.0*log_sigma_segment(0))), 2.0 );
          if(isDouble<Type>::value && of->do_simulate){
            y = exp(rnorm( logmu2 - 0.5*exp(2.0*log_sigma_segment(0)), exp(log_sigma_segment(0)) ));
          }
          break;
        //case poisson_family:
        //  nll -= dpois( y, mu2, true );
        //  devresid += sign(y - mu) * pow(2*(y*log((Type(1e-10) + y)/mu) - (y-mu)), 0.5);
        //  if(isDouble<Type>::value && of->do_simulate){
        //    y = rpois( mu2 );
        //  }
        //  break;
        // case 4:  // Bernoulli
        case gamma_family: // shape = 1/CV^2;   scale = mean*CV^2
          nll -= dgamma( y, exp(-2.0*log_sigma_segment(0)), mu2*exp(2.0*log_sigma_segment(0)), true );
          dev += 2 * ( (y-mu2)/mu2 - log(y/mu2) );
          if(isDouble<Type>::value && of->do_simulate){
            y = rgamma( exp(-2.0*log_sigma_segment(0)), mu2*exp(2.0*log_sigma_segment(0)) );
          }
          break;
        case student_family:  // dnbinom_robust( x, log(mu_i), log(var - mu) )
          // var - mu = exp( 2 * log(mu) - log(theta) ) = mu^2 / theta  -->  var = mu + mu^2 / theta
          nll -= dstudent( y, mu2, exp(log_sigma_segment(0)), 1.0 + exp(log_sigma_segment(1)), true);
          //dev = NAN;
          if(isDouble<Type>::value && of->do_simulate){
            y = mu2 + exp(log_sigma_segment(0)) * rt(1.0 + exp(log_sigma_segment(1)));
          }
          break;
        default:
          error("Distribution not implemented.");
      }
    }
  }

  return mu1 * mu2;
}

}  // namespace tinyVAST
