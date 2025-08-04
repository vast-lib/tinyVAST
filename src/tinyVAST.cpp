#define TMB_LIB_INIT R_init_tinyVAST
#include <TMB.hpp>   //Links in the TMB libraries

enum valid_family {
  gaussian_family  = 0,
  tweedie_family   = 1,
  lognormal_family = 2,
  poisson_family   = 3,
  bernoulli_family = 4,
  binomial_family = 4,
  gamma_family = 5,
  nbinom1_family = 6,
  nbinom2_family = 7
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
Eigen::SparseMatrix<Type> vectorsToSparseMatrix( vector<int> i,
                                                 vector<int> j,
                                                 vector<Type> x,
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

// make stiffness G1 from covariates in columns of E0, E1, E2
template<class Type>
Eigen::SparseMatrix<Type> G_spde_covariates( R_inla::spde_aniso_t<Type> spde,
                                              matrix<Type> H_jj,
                                              matrix<Type> V_zk,
                                              vector<Type> triangle_k ){

  // Extract objects
  vector<Type> area_t = spde.Tri_Area;
  matrix<Type> edge0_tj = spde.E0;
  matrix<Type> edge1_tj = spde.E1;
  matrix<Type> edge2_tj = spde.E2;
  matrix<int> s_tv = spde.TV;         // dim = n_t \times n_v
  matrix<Type> triangle_t1 = V_zk * triangle_k.matrix();
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

    // Make local stiffness
    Gt_vv = (edges_vj * adjH_jj * edges_vj.transpose()) * exp(triangle_t1(t,0));

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
// New matrix-notation precision constructor for tail-down exponential stream network
template<class Type>
Eigen::SparseMatrix<Type> Q_SAR( Type rho,
                                     matrix<Type> H,
                                     int n_s,
                                     vector<int> i_z,
                                     vector<int> j_z,
                                     matrix<Type> delta_z2 ){

  //d_z = sqrt((delta_z2 %*% H)^2 %*% matrix(1,nrow=2,ncol=1))[,1]
  matrix<Type> deltaprime_z2 = delta_z2 * H;
  vector<Type> dist2_z( i_z.size() );
  for( int z = 0; z < i_z.size(); z++ ){
    dist2_z(z) = pow(deltaprime_z2(z,0),2.0) + pow(deltaprime_z2(z,1),2.0);
  }
  vector<Type> weight_z = exp( -1.0 * dist2_z );

  // Make weight matrix
  Eigen::SparseMatrix<Type> W_ss = vectorsToSparseMatrix( i_z,
                                                          j_z,
                                                          weight_z,
                                                          n_s );

  // Row-normalize W_ss
  for( int row = 0; row < n_s; row ++ ){
    Type rowsum = W_ss.row(row).sum();
    W_ss.row(row) /= rowsum;
  }

  // Assemble
  Eigen::SparseMatrix<Type> I_ss( n_s, n_s );
  I_ss.setIdentity();
  Eigen::SparseMatrix<Type> IminusP_ss = I_ss - rho * W_ss;
  Eigen::SparseMatrix<Type> Q = IminusP_ss.transpose() * IminusP_ss;
  return Q;
}

// New matrix-notation precision constructor for tail-down exponential stream network
template<class Type>
Eigen::SparseMatrix<Type> Q_network2( Type log_theta,
                                     int n_s,
                                     vector<int> parent_s,
                                     vector<int> child_s,
                                     vector<Type> dist_s ){

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

//
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
array<Type> omega_distribution( array<Type> omega_sc,
                                 vector<int> model_options,
                                 Eigen::SparseMatrix<Type> Rho_cc,
                                 Eigen::SparseMatrix<Type> Gamma_cc,
                                 Eigen::SparseMatrix<Type> Gammainv_cc,
                                 Eigen::SparseMatrix<Type> Q_ss,
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
        nll += SEPARABLE( GMRF(I_cc), GMRF(Q_ss) )( omega_sc );

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

// distribution/projection for omega
template<class Type>
Type xi_distribution( array<Type> xi_sl,
                      vector<Type> log_sigmaxi_l,
                      Eigen::SparseMatrix<Type> Q_ss ){

  Type nll = 0;
  if( xi_sl.size() > 0 ){
    int n_l = xi_sl.dim(1);
    using namespace density;
    for( int l=0; l<n_l; l++ ){
      nll += SCALE( GMRF(Q_ss), exp(log_sigmaxi_l(l)) )( xi_sl.col(l) );
    }
  }

  return nll;
}

// distribution/projection for epsilon
template<class Type>
array<Type> epsilon_distribution( array<Type> epsilon_stc,
                                  vector<int> model_options,
                                  Eigen::SparseMatrix<Type> Rho_hh,
                                  Eigen::SparseMatrix<Type> Gamma_hh,
                                  Eigen::SparseMatrix<Type> Gammainv_hh,
                                  Eigen::SparseMatrix<Type> Q_ss,
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
    array<Type> epsilon_hs( n_h, n_s );
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
      nll += SEPARABLE( GMRF(Q_ss), GMRF(I_hh) )( epsilon_hs );

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
array<Type> delta_distribution( array<Type> delta_tc,
                                  vector<int> model_options,
                                  Eigen::SparseMatrix<Type> Rho_hh,
                                  Eigen::SparseMatrix<Type> Gamma_hh,
                                  Eigen::SparseMatrix<Type> Gammainv_hh,
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
    array<Type> delta_h1( n_h, 1 );
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
vector<Type> multiply_epsilon( matrix<int> A, vector<Type> weight, array<Type> x, int n_i ){
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
vector<Type> multiply_omega( matrix<int> A, vector<Type> weight, array<Type> x, int n_i ){
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
vector<Type> multiply_xi( Eigen::SparseMatrix<Type> A_is, array<Type> xi_sl, matrix<Type> W_il ){
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
vector<Type> multiply_delta( array<Type> delta_tc, vector<int> t_i, vector<int> c_i, int n_i ){
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
Type dlnorm( Type x,
             Type meanlog,
             Type sdlog,
             int give_log=0){

  //return 1/(sqrt(2*M_PI)*sd) * exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

// dlnorm
template<class Type>
Type dbinom_custom( Type x,
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
Type devresid_binom( Type y,
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
Type devresid_nbinom2( Type y,
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
Type one_predictor_likelihood( Type &y,
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
      default:
        error("Distribution not implemented.");
    }
  }

  return mu;
}

// distribution/projection for epsilon
template<class Type>
Type two_predictor_likelihood( Type y,
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
    if( y == 0 ){
      nll = -1 * log_one_minus_mu1;
      dev = -2 * log_one_minus_mu1;
      if(isDouble<Type>::value && of->do_simulate){
        y = rbinom( Type(1), mu1 );
      }
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
        default:
          error("Distribution not implemented.");
      }
    }
  }

  return mu1 * mu2;
}

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
    // E_gz.col(0) : Expansion Type ( 0=sum(D*Area);  1=weighted.mean(W_gz.col(1),w=D*Area))
    // E_gz.col(1) : Prior row for bivariate weighting

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
    H(i,i) = exp(ln_H_input(i));
  }

  REPORT( H );
  // Spatial settings
  Type log_tau = 0.0;
  Eigen::SparseMatrix<Type> Q_ss;
  // Using INLA with geometric anisotropy
  if( model_options(0)==1 ){
    DATA_STRUCT( spatial_list, R_inla::spde_aniso_t );
    // Build precision
    Q_ss = R_inla::Q_spde( spatial_list, exp(log_kappa), H );
    log_tau = log( 1.0 / (exp(log_kappa) * sqrt(4.0*M_PI)) );
    Type range = pow(8.0, 0.5) / exp( log_kappa );
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
    Q_ss = Q_SAR( rho, H, n_s, i_z, j_z, delta_z2 );
  }
  // Using SPDE with covariate-based anisotropy and geometric anisotropy
  if( model_options(0)==6 ){
    DATA_STRUCT( spatial_list, R_inla::spde_aniso_t );
    DATA_MATRIX( V_zk );
    PARAMETER_VECTOR( triangle_k );
    // Build precision
    Eigen::SparseMatrix<Type> G1 = G_spde_covariates( spatial_list, H, V_zk, triangle_k );
    REPORT( G1 );
    Q_ss = exp(4.0*log_kappa)*spatial_list.G0 + Type(2.0)*exp(2.0*log_kappa)*G1 + G1*spatial_list.G0_inv*G1;
    log_tau = log( 1.0 / (exp(log_kappa) * sqrt(4.0*M_PI)) );
    REPORT( Q_ss );
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
                                   Gamma_cc, Gammainv_cc, Q_ss, nll );
  omega2_sc = omega_distribution( omega2_sc, model_options, Rho2_cc,
                                   Gamma2_cc, Gammainv2_cc, Q_ss, nll );

  // spacetime_term
  epsilon_stc = epsilon_distribution( epsilon_stc, model_options, Rho_hh,
                                     Gamma_hh, Gammainv_hh, Q_ss, nll );
  epsilon2_stc = epsilon_distribution( epsilon2_stc, model_options, Rho2_hh,
                                       Gamma2_hh, Gammainv2_hh, Q_ss, nll );

  // time_term
  delta_tc = delta_distribution( delta_tc, model_options, Rho_time_hh,
                                     Gamma_time_hh, Gammainv_time_hh, nll );
  delta2_tc = delta_distribution( delta2_tc, model_options, Rho2_time_hh,
                                       Gamma2_time_hh, Gammainv2_time_hh, nll );


  // Distribution for spline components
  nll += gamma_distribution( gamma_k, Sdims, Sblock, S_kk, log_lambda );
  nll += gamma_distribution( gamma2_k, S2dims, S2block, S2_kk, log_lambda2 );

  // Distribution for SVC components
  nll += xi_distribution( xi_sl, log_sigmaxi_l, Q_ss );
  nll += xi_distribution( xi2_sl, log_sigmaxi2_l, Q_ss );

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
      vector<Type> phi_g( mu_g.size() );
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
