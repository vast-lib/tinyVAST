
#################
# Simulate bivariate VAR
#################

library(tinyVAST)
library(fmesher)

# Simulate settings
theta_xy = 0.2
n_x = n_y = 10
n_t = 20
B = rbind( c( 0.5, -0.25),
           c(-0.1,  0.50) )

# Simulate GMRFs
R = exp(-theta_xy * abs(outer(1:n_x, 1:n_y, FUN="-")) )
d1 = mvtnorm::rmvnorm(n_t, sigma=0.2*kronecker(R,R) )
d2 = mvtnorm::rmvnorm(n_t, sigma=0.2*kronecker(R,R) )
d = abind::abind( d1, d2, along=3 )

# Project through time and add mean
for( t in seq_len(n_t) ){
  if(t>1) d[t,,] = t(B%*%t(d[t-1,,])) + d[t,,]
}

# Shape into longform data-frame and add error
Data = data.frame( expand.grid(time=1:n_t, x=1:n_x, y=1:n_y, "var"=c("d1","d2")),
                   mu = exp(as.vector(d)))
Data$n = tweedie::rtweedie( n=nrow(Data), mu=Data$mu, phi=0.5, power=1.5 )

# make mesh
mesh = fm_mesh_2d( Data[,c('x','y')] )

# Define DSEM
dsem = "
  d1 -> d1, 1, b11
  d2 -> d2, 1, b22
  d2 -> d1, 1, b21
  d1 -> d2, 1, b12
  d1 <-> d1, 0, var1
  d2 <-> d2, 0, var1
"

# fit model
out = tinyVAST( spacetime_term = dsem,
           data = Data,
           formula = n ~ 0 + var,
           spatial_domain = mesh,
           family = tweedie() )
out

######################
# project
######################

#object = out
#extra_times = 21:25
newdata = Data[1:10,]
newdata$t = 25

project(
  out,
  extra_times = 21:25,
  newdata = newdata
)

#project <-
#function( object,
#          extra_times,
#          newdata,
#          what = "mu_g" ){
#
#
#  ##############
#  # Step 1: Generate uncertainty in historical period
#  ##############
#
#  # SKIP STEP
#
#  ##############
#  # Step 2: Augment objects
#  ##############
#
#  all_times = c( object$internal$times, extra_times )
#
#  ##############
#  # Step 3: Build object with padded bounds
#  ##############
#
#  new_control = object$internal$control
#  new_control$run_model = FALSE
#
#  newobj = tinyVAST(
#    formula = object$formula,
#    data = object$data,
#    time_term = object$internal$time_term,
#    space_term = object$internal$space_term,
#    spacetime_term = object$internal$spacetime_term,
#    family = object$internal$family,
#    space_columns = object$internal$space_columns,
#    spatial_domain = object$spatial_domain,
#    time_column = object$internal$time_column,
#    times = all_times,
#    variable_column = object$internal$variable_column,
#    variables = object$internal$variables,
#    distribution_column = object$internal$distribution_column,
#    delta_options = list( formula = object$internal$delta_formula,
#                          space_term = object$internal$delta_space_term,
#                          time_term = object$internal$delta_time_term,
#                          spacetime_term = object$internal$delta_spacetime_term,
#                          spatial_varying = object$internal$delta_spatial_varying ),
#    spatial_varying = object$internal$spatially_varying,
#    weights = object$internal$weights,
#    control = new_control
#  )
#
#  ##############
#  # Step 4: Merge ParList and ParList1
#  ##############
#
#  augment_epsilon <-
#  function( neweps_stc,
#            eps_stc,
#            beta_z,
#            model ){
#
#    if( length(beta_z) > 0 ){
#      #
#      mats = dsem::make_matrices(
#        beta_p = beta_z,
#        model = model,
#        variables = object$internal$variables,
#        times = all_times
#      )
#      Q_kk = Matrix::t(mats$IminusP_kk) %*% (Matrix::t(mats$G_kk) %*% mats$G_kk) %*% mats$IminusP_kk
#      Q_hh = Matrix::kronecker( Q_kk, Q_ss )
#
#      #
#      grid = expand.grid( s = seq_len(dim(neweps_stc)[1]),
#                          t = seq_len(dim(neweps_stc)[2]),
#                          c = seq_len(dim(neweps_stc)[3]) )
#      grid$num = seq_len(prod(dim(neweps_stc)))
#      observed_idx = subset( grid, t %in% object$internal$times )$num
#
#      #
#      simeps_stc = simulate_conditional_gmrf(
#        Q = Q_hh,
#        observed_idx = observed_idx,
#        x_obs = as.vector( eps_stc ),
#        n_sims = 1
#      )
#
#      # Compile
#      missing_indices = as.matrix(subset( grid, t %in% extra_times )[,1:3])
#      neweps_stc[missing_indices] = simeps_stc[,1]
#      observed_indices = as.matrix(subset( grid, t %in% object$internal$times )[,1:3])
#      neweps_stc[observed_indices] = eps_stc[observed_indices]
#    }
#    return(neweps_stc)
#  }
#  augment_delta <-
#  function( newdelta_tc,
#            delta_tc,
#            nu_z,
#            model ){
#
#    if( length(nu_z) > 0 ){
#      #
#      mats = dsem::make_matrices(
#        beta_p = nu_z,
#        model = model,
#        variables = object$internal$variables,
#        times = all_times
#      )
#      Q_kk = Matrix::t(mats$IminusP_kk) %*% (Matrix::t(mats$G_kk) %*% mats$G_kk) %*% mats$IminusP_kk
#
#      #
#      grid = expand.grid( t = seq_len(dim(newdelta_tc)[1]),
#                          c = seq_len(dim(newdelta_tc)[2]) )
#      grid$num = seq_len(prod(dim(newdelta_tc)))
#      observed_idx = subset( grid, t %in% object$internal$times )$num
#
#      #
#      simdelta_tc = simulate_conditional_gmrf(
#        Q = Q_kk,
#        observed_idx = observed_idx,
#        x_obs = as.vector( delta_tc ),
#        n_sims = 1
#      )
#
#      # Compile
#      missing_indices = as.matrix(subset( grid, t %in% extra_times )[,1:2])
#      newdelta_tc[missing_indices] = simdelta_stc[,1]
#      observed_indices = as.matrix(subset( grid, t %in% object$internal$times )[,1:2])
#      newdelta_tc[observed_indices] = delta_tc[observed_indices]
#    }
#    return(newdelta_tc)
#  }
#
#  #
#  parlist = object$internal$parlist
#  new_parlist = newobj$tmb_par
#  Q_ss = object$rep$Q_ss
#
#  # Replace epsilon
#  new_parlist$epsilon_stc = augment_epsilon(
#    beta_z = parlist$beta_z,
#    eps_stc = parlist$epsilon_stc,
#    neweps_stc = new_parlist$epsilon_stc,
#    model = object$internal$spacetime_term_ram$output$model
#  )
#  new_parlist$epsilon2_stc = augment_epsilon(
#    beta_z = parlist$beta2_z,
#    eps_stc = parlist$epsilon2_stc,
#    neweps_stc = new_parlist$epsilon2_stc,
#    model = object$internal$delta_spacetime_term_ram$output$model
#  )
#
#  # Replace delta
#  new_parlist$delta_tc = augment_delta(
#    nu_z = parlist$nu_z,
#    delta_tc = parlist$delta_tc,
#    newdelta_tc = new_parlist$delta_tc,
#    model = object$internal$time_term_ram$output$model
#  )
#  new_parlist$delta2_tc = augment_delta(
#    nu_z = parlist$nu2_z,
#    delta_tc = parlist$delta2_tc,
#    newdelta_tc = new_parlist$delta2_tc,
#    model = object$internal$delta_time_term_ram$output$model
#  )
#
#  # Replace other variables that are not changed
#  same_vars = setdiff( names(new_parlist), c("epsilon_stc","epsilon2_stc","delta_tc","delta2_tc") )
#  new_parlist[same_vars] = parlist[same_vars]
#
#  ##############
#  # Step 5: Re-build model
#  ##############
#
#  new_control$run_model = TRUE
#  new_control$tmb_par = new_parlist
#  new_control$nlminb_loops = 0
#  new_control$newton_loops = 0
#  new_control$getsd = FALSE
#  new_control$calculate_deviance_explained = FALSE
#
#  newobj = tinyVAST(
#    formula = object$formula,
#    data = object$data,
#    time_term = object$internal$time_term,
#    space_term = object$internal$space_term,
#    spacetime_term = object$internal$spacetime_term,
#    family = object$internal$family,
#    space_columns = object$internal$space_columns,
#    spatial_domain = object$spatial_domain,
#    time_column = object$internal$time_column,
#    times = all_times,
#    variable_column = object$internal$variable_column,
#    variables = object$internal$variables,
#    distribution_column = object$internal$distribution_column,
#    delta_options = list( formula = object$internal$delta_formula,
#                          space_term = object$internal$delta_space_term,
#                          time_term = object$internal$delta_time_term,
#                          spacetime_term = object$internal$delta_spacetime_term,
#                          spatial_varying = object$internal$delta_spatial_varying ),
#    spatial_varying = object$internal$spatially_varying,
#    weights = object$internal$weights,
#    control = new_control
#  )
#
#  ##############
#  # Step 6: simulate samples
#  ##############
#
#  pred = predict(
#    object = newobj,
#    newdata = newdata,
#    what = what
#  )
#  return(pred)
#}
