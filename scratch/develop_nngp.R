

# Test compile locally
if( FALSE ){
  setwd( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\src)' )
  TMB::compile('tinyVAST.cpp', framework = "TMBad")
}

######################
# Provide areal domain
######################

library(sf)

boundary = st_polygon(list(cbind(c(0,1,1,0,0),c(0,0,1,1,0))))
spatial_domain = st_make_grid( boundary, square = FALSE, n = c(4,4) )

######################
# Make NNGP neighborhoods
######################

library(Matrix)

#data = make_nngp_data(
#  coords = st_coordinates(st_centroid(spatial_domain)),
#  nn = 4
#)

st_adjacent <- function(m, ...) st_relate(m, m, pattern="F***1****", ...)
grid_A = st_adjacent(spatial_domain, sparse=TRUE)
A_ss = as(grid_A, "sparseMatrix")
A_ss = sweep(A_ss, MARGIN=1, FUN = "/", STAT = rowSums(A_ss))
IminusA_ss = Diagonal( n = length(spatial_domain) ) - 0.9 * A_ss
Q_ss = t(IminusA_ss) %*% IminusA_ss

omega_i = RTMB:::rgmrf0( n = 1, Q_ss )[,1]
plot( st_sf(spatial_domain,omega_i) )
p_i = 3 + omega_i - mean(omega_i)
y_i = rpois( n = length(spatial_domain), lambda = exp(p_i) )

##################
# Run locally
##################
library(checkmate)
library(TMB)
library(GpGp)
library(tinyVAST)
#source( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\R\utility.R)' )
source( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\R\internal.R)' )
source( R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\R\fit.R)' )

formula = response ~ 1
data = data.frame(
  response = y_i,
  setNames(data.frame(st_coordinates(st_centroid(spatial_domain))), c("x","y"))
)
#space_term = NULL
space_term = ""
family = poisson()

#formula
#data
time_term = NULL
#space_term = NULL
spacetime_term = NULL
#family = gaussian()
delta_options = list(formula = ~ 1)
spatial_varying = NULL
weights = NULL
#spatial_domain = NULL
development = list()
control = tinyVASTcontrol(
  nearest_neighbors = 4,
  gmrf_parameterization = "projection",
  extra_reporting = TRUE
)
# Indexing
space_columns = c("x","y")
time_column = "time"
times = NULL
variable_column = "var"
variables = NULL
distribution_column = "dist"

#######################
# COPY
#######################

  # https://roxygen2.r-lib.org/articles/rd-formatting.html#tables for roxygen formatting
  start_time = Sys.time()

  # General error checks
  assertClass(control, "tinyVASTcontrol")
  assertDataFrame(data)
  if(inherits(data,"tbl")) stop("`data` must be a data.frame and cannot be a tibble")

  # Add `development` stuff
  if( !is.null(development$vertex_formula) ){
    vertex_formula = development$vertex_formula
  }else{
    vertex_formula = ~ 0
  }
  if( !is.null(development$triangle_formula) ){
    triangle_formula = development$triangle_formula
  }else{
    triangle_formula = ~ 0
  }

  # Avoid name conflict
  if( "fake" %in% colnames(spatial_domain$triangle_covariates) ){
    stop("change colnames in `triangle_covariates` to avoid `fake`")
  }

    # Input conflicts
  matched_call = match.call()
  if( isTRUE(as.character(matched_call$family) == "family") ){
    stop("Naming argument `family` as `family` conflicts with function `cv::cv`, please use `family = Family` or other name", call. = FALSE)
  }
  if( !is(spatial_domain,"vertex_coords") ){
    if( vertex_formula != as.formula("~0") ){
      stop("specifying `vertex_formula` only makes sense when `spatial_domain` has class `vertex_coords`", call. = FALSE)
    }
    if( triangle_formula != as.formula("~0") ){
      stop("specifying `triangle_formula` only makes sense when `spatial_domain` has class `vertex_coords`", call. = FALSE)
    }
  }
  if( !is.null(control$nearest_neighbors) & (control$gmrf_parameterization != "projection") ){
    stop("`nearest neighbors` only works with `projection` gmrf_parameterization")
  }

  # Haven't tested for extra levels
  tmpdata = droplevels(data)
  if( !identical(tmpdata,data) ){
    stop("`data` has some factor with extra levels. Please retry after running `data = droplevels(data)` on the input `data`", call. = FALSE)
  }

  ##############
  # input telescoping
  ##############

  # Telescope family ... comes before adding `distribution` to `data`
  if( inherits(family,"family") ){
    family = list( "obs"=family )
  }

  # Defaults for missing columns of data
  # character so that make_dsem_ram arg covs has a character
  if( isFALSE(variable_column %in% colnames(data)) ){
    data = data.frame(data, matrix("response", nrow=nrow(data), ncol=1, dimnames=list(NULL,variable_column)))
  }
  if( isFALSE(time_column %in% colnames(data)) ){
    data = data.frame( data, matrix(1, nrow=nrow(data), ncol=1, dimnames=list(NULL,time_column)) )
  }
  if( isFALSE(distribution_column %in% colnames(data)) ){
    if( length(family)>1 ) stop("Must supply `dist` if using multiple `family` options")
    data = data.frame( data, matrix(names(family)[1], nrow=nrow(data), ncol=1, dimnames=list(NULL,distribution_column)) )
  }

  # Defaults for times
  if( is.null(spacetime_term) & is.null(delta_options$spacetime_term) & is.null(time_term) & is.null(delta_options$time_term) ){
    times = numeric(0)
  }
  if(is.null(times)) times = seq( min(data[,time_column]), max(data[,time_column]) )

  # Defaults for variables
  if( is.null(time_term) & is.null(spacetime_term) & is.null(space_term) & is.null(spatial_varying) &
      is.null(delta_options$spacetime_term) & is.null(delta_options$space_term) & is.null(delta_options$time_term) & is.null(delta_options$spatial_varying) ){
    variables = numeric(0)
    # Force spatial_domain = NULL so that estimate_kappa=FALSE
    spatial_domain = NULL
  }
  if(is.null(variables)) variables = unique( data[,variable_column] )

  # Turn of t_i and c_i when times and variables are missing, so that delta_k isn't built
  if( length(times) > 0 ){
    t_i = match( data[,time_column], times )
  }else{
    t_i = integer(0)
  }
  if( length(variables) > 0 ){
    c_i = match( data[,variable_column], variables )
    # variables can't have commas, because it conflicts with how `space_term` and `spacetime_term` are parsed
    if( any(grepl(pattern=",", x=variables, fixed=TRUE)) ){
      stop("`variables` cannot include any commas")
    }
  }else{
    c_i = integer(0)
  }

  # Deal with likelihood weights
  if( is.null(weights) ){
    weights_i = rep(1,nrow(data))
  }else{
    assertNumeric(weights, lower=0, finite=TRUE, len=nrow(data), any.missing=FALSE)
    weights_i = weights
  }

  ##############
  # spacetime_term_ram constructor
  ##############

  build_spacetime_term <-
  function( spacetime_term ){
    # (I-Rho)^-1 * Gamma * (I-Rho)^-1
    if( is.null(spacetime_term) ){
      output = list(
        ram = array( 0, dim=c(0,5), dimnames=list(NULL,c("heads","to","from","parameter","start")) ),
        model = array( 0, dim=c(0,8), dimnames=list(NULL,c("path","lag","name","start","parameter","first","second","direction")) )
      )
      class(output) = "dsem_ram"   # Define class for summary.tinyVAST
    }else if( isTRUE(is.character(spacetime_term)) ){
      output = make_dsem_ram( spacetime_term, times=times, variables=variables, quiet=isFALSE(control$verbose), covs=variables )
    }else if( is(spacetime_term,"spacetime_term_ram") | is(spacetime_term,"eof_ram") ){
      output = spacetime_term
    }else{
      stop("`spacetime_term` must be either `NULL` or a character-string", call. = FALSE)
    }

    # Identify arrow-type for each beta_j estimated in RAM
    which_nonzero = which(output$ram[,4]>0)
    param_type = tapply( output$ram[which_nonzero,1],
                         INDEX=output$ram[which_nonzero,4], FUN=max)

    # Error checks
    if( is(output, "dsem_ram") ){
      if( any((output$model[,'direction']==2) & (output$model[,2]!=0)) ){
        stop("All two-headed arrows should have lag=0")
      }
      if( isFALSE(all(c(output$model[,'first'],output$model[,'second']) %in% variables)) ){
        stop("Some variable in `spacetime_term` is not in `tsdata`")
      }
    }

    # Check for rank-deficient precision from RAM
    df_ram = data.frame(output$ram)
    ram_gamma = df_ram[df_ram$heads==2,,drop=FALSE]
    total_variance_h = tapply( as.numeric(ifelse(is.na(ram_gamma$start), 1, ram_gamma$start)),
            INDEX = ram_gamma$from, FUN=function(x)sum(abs(x)) )
    ram_rho = df_ram[df_ram$heads==1,,drop=FALSE]
    total_effect_h = tapply( as.numeric(ifelse(is.na(ram_rho$start), 1, ram_rho$start)),
            INDEX = ram_rho$from, FUN=function(x)sum(abs(x)) )
    if( any(total_variance_h==0) && control$gmrf_parameterization=="separable" ){
      stop("Must use gmrf_parameterization=`projection` for the spacetime_term RAM supplied")
    }

    # out
    out = list("output"=output, "param_type"=param_type)
    return(out)
  }
  spacetime_term_ram = build_spacetime_term(spacetime_term)
  delta_spacetime_term_ram = build_spacetime_term( delta_options$spacetime_term )

  ##############
  # time_term_ram constructor
  ##############

  time_term_ram = build_spacetime_term(time_term)
  delta_time_term_ram = build_spacetime_term( delta_options$time_term )

  ##############
  # space_term_ram constructor
  ##############

  build_space_term <-
  function( space_term ){
    if( is.null(space_term) ){
      output = list(
        ram = array( 0, dim=c(0,5), dimnames=list(NULL,c("heads","to","from","parameter","start")) ),
        model = array( 0, dim=c(0,8), dimnames=list(NULL,c("","","","","parameter","first","second","direction")) )
      )
      class(output) = "sem_ram" # Define class for summary.tinyVAST
    }else if( isTRUE(is.character(space_term)) ){
      output = make_sem_ram( space_term, variables=as.character(variables), quiet=isFALSE(control$verbose), covs=as.character(variables) )
    } else {
      stop("`space_term` must be either `NULL` or a character-string")
    }

    # Identify arrow-type for each beta_j estimated in RAM
    which_nonzero = which(output$ram[,4]>0)
    param_type = tapply( output$ram[which_nonzero,1],
                         INDEX=output$ram[which_nonzero,4], FUN=max)

    # Check for rank-deficient precision from RAM
    df_ram = data.frame(output$ram)
    ram2 = df_ram[df_ram$heads == 2,,drop=FALSE]
    total_variance_h = tapply( as.numeric(ifelse( is.na(ram2$start), 1, ram2$start)),
            INDEX = ram2$from, FUN=function(x)sum(abs(x)) )
    ram1 = df_ram[df_ram$heads == 1,,drop=FALSE]
    total_effect_h = tapply( as.numeric(ifelse(is.na(ram1$start), 1, ram1$start)),
            INDEX = ram1$from, FUN=function(x)sum(abs(x)) )
    if( any(total_variance_h==0) && control$gmrf_parameterization=="separable" ){
      stop("Must use options$gmrf_parameterization=`projection` for the space_term RAM supplied")
    }

    # out
    out = list("output"=output, "param_type"=param_type)
    return(out)
  }
  space_term_ram = build_space_term(space_term)
  delta_space_term_ram = build_space_term( delta_options$space_term )

  ##############
  # Spatial domain constructor
  ##############

  # n_Hpars is the number of params in H AND the dim of H
  n_Hpars = 2

  # Empty default
  nngp_data = make_nngp_data( what = "empty" )

  # Parse each `spatial_domain`
  if( is(spatial_domain,"vertex_coords") ){
    updated_vertex_formula = update.formula( old = vertex_formula, new = "~ . + 0" )
    Terms <- terms(updated_vertex_formula, data = spatial_domain$vertex_covariates)
    mf <- model.frame(
      updated_vertex_formula,
      data = spatial_domain$vertex_covariates,
      na.action = na.pass,
      drop.unused.levels = TRUE
    )
    # Parse covariate matrix
    R_sk = model.matrix(
      Terms,
      data = mf,
      na.action = na.pass
    )
    #V_zk = model.matrix(
    #  update.formula( old = ~ AK + GOA + BS, new = "~ . + 0" ),
    #  spatial_domain$triangle_covariates
    #)
    if( triangle_formula != as.formula("~0") ){
      gam_setup = mgcv::gam( update.formula( old = triangle_formula, new = "fake ~ . + 0" ),
                             data = cbind( "fake" = 0, spatial_domain$triangle_covariates ),
                             fit = FALSE )
      V_zk = cbind( offset = gam_setup$offset, gam_setup$X ) # First is always the offset
    }else{
      V_zk = matrix( 0, ncol = 1, nrow = nrow(spatial_domain$triangle_covariates) )
    }
    # covariate-based anisotropy
    n_s = spatial_domain$n
    spatial_method_code = 6
    spatial_list = make_anisotropy_spde( spatial_domain, covariates = R_sk )
    spatial_list$V_zk = V_zk
    A_is = fm_evaluator( spatial_domain, loc=as.matrix(data[,space_columns]) )$proj$A
    n_Hpars = ncol( spatial_list$E0 )
    estimate_kappa = TRUE
    kappa_startvalue = 1
    estimate_anisotropy = ifelse( isTRUE(control$use_anisotropy), TRUE, FALSE )
  }else if( is(spatial_domain,"fm_mesh_2d") ){
    # SPDE
    n_s = spatial_domain$n
    spatial_method_code = 1
    spatial_list = make_anisotropy_spde( spatial_domain )
    A_is = fm_evaluator( spatial_domain, loc=as.matrix(data[,space_columns]) )$proj$A
    n_Hpars = 2
    estimate_kappa = TRUE
    kappa_startvalue = sqrt(8) / mean(apply( data[,space_columns], MARGIN = 2, FUN = sd) )      # range = sqrt(8) / kappa
    estimate_anisotropy = ifelse( isTRUE(control$use_anisotropy), TRUE, FALSE )
  }else if( is(spatial_domain,"igraph") ) {
    # SAR ... not doing row-standardized, so currently can result in nonsense Q
    spatial_method_code = 2
    Adj = as_adjacency_matrix( spatial_domain, sparse=TRUE )
    n_s = nrow(Adj)
    Match = match( data[,space_columns], rownames(Adj) )
    if(any(is.na(Match))) stop("Check `spatial_domain` for SAR")
    A_is = sparseMatrix( i=1:nrow(data), j=Match, x=rep(1,nrow(data)) )
    # Turn off log_kappa if no edges (i.e., unconnected graph)
    estimate_kappa = ifelse( ecount(spatial_domain)>0, TRUE, FALSE )
    kappa_startvalue = 1
    estimate_anisotropy = FALSE
  }else if( is(spatial_domain,"sfnetwork_mesh") ){      # if( !is.null(space_term) )
    # stream network
    spatial_method_code = 4
    n_s = spatial_domain$n
    spatial_list = spatial_domain$table
    A_is = sfnetwork_evaluator( stream = spatial_domain$stream,
                                loc = as.matrix(data[,space_columns]) )
    estimate_kappa = TRUE
    kappa_startvalue = 1
    estimate_anisotropy = FALSE
  }else if( is_areal_sf(spatial_domain) & is.null(control$nearest_neighbors) ){
    # SAR with geometric anisotropy
    spatial_method_code = 5
    n_s = length(spatial_domain)
    if( control$sar_adjacency == "rook" ){
      st_adjacent <- function(m, ...) st_relate(m, m, pattern="F***1****", ...)
    }else{
      st_adjacent <- function(m, ...) st_relate(m, m, pattern = "F***T****", ...)
    }
    grid_A = st_adjacent(spatial_domain, sparse=TRUE)
    A_ss = as(grid_A, "sparseMatrix")
    grid_xy = st_coordinates(st_centroid(spatial_domain))
    tripA = mat2triplet(A_ss)
    spatial_list = list(
      i_z = tripA$i,
      j_z = tripA$j,
      delta_z2 = grid_xy[tripA$i,] - grid_xy[tripA$j,]
    )
    sf_coords = st_as_sf( data,
                          coords = space_columns,
                          crs = st_crs(spatial_domain) )
    s_i = as.integer(st_within( sf_coords, spatial_domain ))
    A_is = sparseMatrix( i = seq_along(s_i),
                         j = s_i,
                         x = 1,
                         dims = c(length(s_i),length(spatial_domain)) )
    estimate_kappa = TRUE
    kappa_startvalue = 1
    estimate_anisotropy = ifelse( isTRUE(control$use_anisotropy), TRUE, FALSE )
  }else if( is_areal_sf(spatial_domain) & !is.null(control$nearest_neighbors) ){
    spatial_method_code = 7
    n_s = length(spatial_domain)
    nngp_data = make_nngp_data(
      coords = st_coordinates(st_centroid(spatial_domain)),
      nn = control$nearest_neighbors
    )
    sf_coords = st_as_sf( data,
                          coords = space_columns,
                          crs = st_crs(spatial_domain) )
    s_i = as.integer(st_within( sf_coords, spatial_domain ))
    A_is = sparseMatrix( i = seq_along(s_i),
                         j = s_i,
                         x = 1,
                         dims = c(length(s_i),length(spatial_domain)) )
    estimate_kappa = TRUE
    kappa_startvalue = 1
    estimate_anisotropy = FALSE
  }else if( is.null(spatial_domain) ) {
    # Single-site
    spatial_method_code = 3
    n_s = 1
    A_is = matrix(1, nrow=nrow(data), ncol=1)    # dgCMatrix
    A_is = as(Matrix(A_is),"CsparseMatrix")
    spatial_list = list( "M0" = as(Matrix(1,nrow=1,ncol=1),"CsparseMatrix"),
                         "M1" = as(Matrix(0,nrow=1,ncol=1),"CsparseMatrix"),
                         "M2" = as(Matrix(0,nrow=1,ncol=1),"CsparseMatrix") )
    estimate_kappa = FALSE
    kappa_startvalue = 1
    estimate_anisotropy = FALSE
  }else{
    stop("`spatial_domain` is does not match options:  class fm_mesh_2d, igraph, sfnetwork_mesh, sfc_GEOMETRY, or NULL")
  }
  Atriplet = Matrix::mat2triplet(A_is)
  if( (n_s>1000) & isFALSE(control$suppress_user_warnings) ){
    warning("`spatial_domain` has over 1000 components, so the model may be extremely slow")
  }

  #
  Aepsilon_zz = cbind(Atriplet$i, Atriplet$j, t_i[Atriplet$i], c_i[Atriplet$i])
  which_Arows = which(apply( Aepsilon_zz, MARGIN=1, FUN=function(x) all(!is.na(x)) & any(x>0) ))
  which_Arows = which_Arows[ which(Atriplet$x[which_Arows] > 0) ]
  if( (nrow(spacetime_term_ram$output$ram)==0) & (nrow(delta_spacetime_term_ram$output$ram)==0) ){
    which_Arows = numeric(0)
  }
  Aepsilon_zz = Aepsilon_zz[which_Arows,,drop=FALSE]
  Aepsilon_z = Atriplet$x[which_Arows]

  #
  Aomega_zz = cbind(Atriplet$i, Atriplet$j, c_i[Atriplet$i])
  which_Arows = which(apply( Aomega_zz, MARGIN=1, FUN=function(x) all(!is.na(x)) ))
  which_Arows = which_Arows[ which(Atriplet$x[which_Arows] > 0) ]
  if( (nrow(space_term_ram$output$ram)==0) & (nrow(delta_space_term_ram$output$ram)==0) ){
    which_Arows = numeric(0)
  }
  Aomega_zz = Aomega_zz[which_Arows,,drop=FALSE]
  Aomega_z = Atriplet$x[which_Arows]

  ##############
  # Formula constructor
  ##############

  # Initial constructor of splines
  build_gam_basis <-
  function( formula ){
    gam_setup = mgcv::gam( formula, data = data, fit=FALSE ) # select doesn't do anything in this setup
    y_i = model.response(gam_setup$mf)  # OR USE: model.extract(gam_setup$mf, "response")
    offset_i = gam_setup$offset

    # Extract and combine penalization matrices
    #S_list = lapply( seq_along(gam_setup$smooth), function(x) gam_setup$smooth[[x]]$S[[1]] )
    #S_kk = Matrix::.bdiag(S_list)       # join S's in sparse matrix
    #Sdims = unlist(lapply(S_list,nrow)) # Find dimension of each S
    #if(is.null(Sdims)) Sdims = vector(length=0)

    # Warnings and errors
    not_allowed <- vapply( c("t2("), function(.x)
      length(grep(.x, x=gam_setup$term.names, fixed=TRUE)) > 0, FUN.VALUE = logical(1L)
    )
    if(any(not_allowed)) {
      stop("Found t2() smoothers. These are not yet implemented.", call. = FALSE)
    }
    which_experimental <- vapply( c("ti(", "te("), function(.x)
      length(grep(.x, x=gam_setup$term.names, fixed=TRUE)) > 0, FUN.VALUE = logical(1L)
    )
    if(any(which_experimental)) {
      if(isTRUE(control$verbose)) message("Found ti() or te() smoothers. These are experimental.")
    }

    # Extract and combine penalization matrices by block
    S_list = lapply( seq_along(gam_setup$smooth), function(x) gam_setup$smooth[[x]]$S )
    S_kk = Matrix::.bdiag( lapply(S_list, Matrix::.bdiag) )
    Sdims = unlist( lapply(S_list, FUN = function(List){nrow(List[[1]])}) )
    Sblock = unlist( lapply(S_list, length) )
    if( sum(Sdims * Sblock) != nrow(S_kk) ) stop("Check constructor for smooths, term `Sblock`")

    # Identify fixed vs. random effects
    search_names = c( "s(", "te(", "ti(", "t2(" )
    which_se = sapply( X = search_names,
                       FUN = function(x){ grep( pattern=x, x=gam_setup$term.names, fixed=TRUE ) } )
    which_se = sort( unlist(which_se), decreasing = FALSE )
    if( length(which_se) != sum(Sdims) ) stop("Check constructor for smooths, term `which_se`")

    # Extract design matrices for fixed and random effects and add names
    colnames(gam_setup$X) = gam_setup$term.names
    X_ij = gam_setup$X[,setdiff(seq_len(ncol(gam_setup$X)),which_se),drop=FALSE]
    Z_ik = gam_setup$X[,which_se,drop=FALSE]
    if(is.null(Sdims)) Sdims = numeric()
    if(is.null(Sblock)) Sblock = numeric()

    #
    out = list( "X_ij"=X_ij, "Z_ik"=Z_ik, "S_kk"=S_kk, "Sdims"=Sdims, "Sblock" = Sblock,
                "y_i"=y_i, "offset_i"=offset_i, "gam_setup"=gam_setup )
    return(out)
  }
  gam_basis = build_gam_basis( formula )
  delta_formula_with_response = update.formula( formula,
                                  paste(".", paste0(as.character(delta_options$formula),collapse="")) )
  delta_gam_basis = build_gam_basis( delta_formula_with_response )

  ##############
  # distribution/link
  ##############

  build_distributions <-
  function( family ){
    e_i = match( data[,distribution_column], names(family) )

    # Check for errors
    if( (any(is.na(e_i))) ){
      stop("`data[,distribution_column]` has values that don't match `names(family)`")
    }

    # Construct log_sigma based on family
    pad_length = function(x){if(length(x)==1) c(x,99L) else x}
    remove_last = function(x) x[-length(x)]
    #Nsigma_e = sapply( family, FUN=function(x){
    #                   switch( x$family[length(x$family)],
    #                     "gaussian" = 1,
    #                     "tweedie" = 2,
    #                     "lognormal" = 1,
    #                     "poisson" = 0,
    #                     "nbinom2" = 1,
    #                     "nbinom1" = 1,
    #                     "binomial" = 0,
    #                     "bernoulli" = 0,
    #                     "Gamma" = 1,
    #                     "student" = 2
    #                   )} )

    # Fixed values
    sigma_e = lapply( family, FUN=function(x){
                       switch( x$family[length(x$family)],
                         "gaussian" = NA,
                         "tweedie" = c(NA,NA),
                         "lognormal" = NA,
                         "poisson" = c(),
                         "nbinom2" = NA,
                         "nbinom1" = NA,
                         "binomial" = c(),
                         "bernoulli" = c(),
                         "Gamma" = c(NA),
                         "student" = c(NA, ifelse(is.null(x$df), NA, log(x$df-1) ) )
                       )} )
    Nsigma_e = sapply(sigma_e, length)
    Edims_ez = cbind( "start"=remove_last(cumsum(c(0,Nsigma_e))), "length"=Nsigma_e )

    #
    family_code = t(rbind(sapply( family, FUN=function(x){
                       pad_length(c("gaussian" = 0,
                         "tweedie" = 1,
                         "lognormal" = 2,
                         "poisson" = 3,
                         "binomial" = 4,
                         "bernoulli" = 4,
                         "Gamma" = 5,
                         "nbinom1" = 6,
                         "nbinom2" = 7,
                         "student" = 8)[x$family])
                       } )))
    link_code = t(rbind(sapply( family, FUN=function(x){
                       pad_length(c("identity" = 0,
                         "log" = 1,
                         "logit" = 2,
                         "cloglog" = 3 )[x$link])
                       } )))
    components = apply( family_code, MARGIN=1,
                         FUN=function(x)sum(x != 99L) )
    poisson_link_delta = sapply( family, FUN=function(x){
                         as.integer(isTRUE(x$type == "poisson_link_delta"))} )
    out = list( "family_code" = cbind(family_code),
                "link_code" = cbind(link_code),
                "components" = components,
                "poisson_link_delta" = poisson_link_delta,
                "e_i" = e_i,
                "Nsigma_e" = Nsigma_e,
                "Edims_ez" = Edims_ez,
                "sigma_e" = sigma_e )
    return(out)
  }
  distributions = build_distributions( family )

  # Convert weights to size for binomial as special case
  is_binomial = (distributions$family_code[distributions$e_i,1] == 4) & (distributions$family_code[distributions$e_i,2] == 99)
  size_i = ifelse( is_binomial, weights_i, 1 )
  weights_i = ifelse( is_binomial, 1, weights_i )

  ##############
  # Build spatially varying
  ##############

  build_spatial_varying <-
  function( spatial_varying ){
    if( is.null(spatial_varying) ){
      spatial_varying = ~ 0
    }
    W_il = model.matrix( spatial_varying, data = data )

    out = list( "W_il" = W_il )
    return(out)
  }
  SVC = build_spatial_varying( spatial_varying )
  delta_SVC = build_spatial_varying( delta_options$spatial_varying )

  ##############
  # Build inputs
  # All interactions among features should come here
  ##############

  # Make options
  model_options = c(
    spatial_method_code,
    ifelse( control$gmrf_parameterization=="separable", 0, 1),
    ifelse( isFALSE(control$get_rsr), 0, 1),
    ifelse( isFALSE(control$extra_reporting), 0, 1),
    0
  )

  # make dat
  tmb_data = list(
    model_options = model_options,
    y_i = gam_basis$y_i,
    X_ij = gam_basis$X_ij,
    Z_ik = gam_basis$Z_ik,
    W_il = SVC$W_il,
    X2_ij = delta_gam_basis$X_ij,
    Z2_ik = delta_gam_basis$Z_ik,
    W2_il = delta_SVC$W_il,
    t_i = ivector_minus_one(t_i, "t_i"), # -1 to convert to CPP index, and keep as integer-vector
    c_i = ivector_minus_one(c_i, "c_i"), # -1 to convert to CPP index, and keep as integer-vector
    offset_i = gam_basis$offset_i,
    weights_i = weights_i,
    size_i = size_i,
    family_ez = distributions$family_code,
    link_ez = distributions$link_code,
    components_e = distributions$components,
    poislink_e = distributions$poisson_link_delta,
    e_i = ivector_minus_one(distributions$e_i, "e_i"), # -1 to convert to CPP index, and keep as integer-vector
    Edims_ez = distributions$Edims_ez,
    S_kk = gam_basis$S_kk,
    Sdims = gam_basis$Sdims,
    Sblock = gam_basis$Sblock,
    S2_kk = delta_gam_basis$S_kk,
    S2dims = delta_gam_basis$Sdims,
    S2block = delta_gam_basis$Sblock,
    Aepsilon_zz = Aepsilon_zz - 1,     # Index form, i, s, t
    Aepsilon_z = Aepsilon_z,
    Aomega_zz = Aomega_zz - 1,     # Index form, i, s, t
    Aomega_z = Aomega_z,
    A_is = A_is,
    nngp_data = nngp_data,
    ram_space_term = as.matrix(na.omit(space_term_ram$output$ram[,1:4])),
    ram_space_term_start = as.numeric(space_term_ram$output$ram[,5]),
    ram_time_term = as.matrix(na.omit(time_term_ram$output$ram[,1:4])),
    ram_time_term_start = as.numeric(time_term_ram$output$ram[,5]),
    ram_spacetime_term = as.matrix(na.omit(spacetime_term_ram$output$ram[,1:4])),
    ram_spacetime_term_start = as.numeric(spacetime_term_ram$output$ram[,5]),
    ram2_space_term = as.matrix(na.omit(delta_space_term_ram$output$ram[,1:4])),
    ram2_space_term_start = as.numeric(delta_space_term_ram$output$ram[,5]),
    ram2_time_term = as.matrix(na.omit(delta_time_term_ram$output$ram[,1:4])),
    ram2_time_term_start = as.numeric(delta_time_term_ram$output$ram[,5]),
    ram2_spacetime_term = as.matrix(na.omit(delta_spacetime_term_ram$output$ram[,1:4])),
    ram2_spacetime_term_start = as.numeric(delta_spacetime_term_ram$output$ram[,5]),
    X2_gj = matrix(0,ncol=ncol(delta_gam_basis$X_ij),nrow=0),
    Z2_gk = matrix(0,ncol=ncol(delta_gam_basis$Z_ik),nrow=0),
    W2_gl = matrix(0,ncol=ncol(delta_SVC$W_il),nrow=0),
    X_gj = matrix(0,ncol=ncol(gam_basis$X_ij),nrow=0),
    Z_gk = matrix(0,ncol=ncol(gam_basis$Z_ik),nrow=0),
    W_gl = matrix(0,ncol=ncol(SVC$W_il),nrow=0),
    AepsilonG_zz = matrix(0,nrow=0,ncol=4),
    AepsilonG_z = numeric(0),
    AomegaG_zz = matrix(0,nrow=0,ncol=4),
    AomegaG_z = numeric(0),
    A_gs = sparseMatrix(i=integer(), j=integer(), x=numeric(), dims=c(0,ncol(A_is))),
    t_g = integer(0),
    c_g = integer(0),
    offset_g = integer(0),
    e_g = integer(0),
    W_gz = matrix(0,nrow=0,ncol=2),
    V_gz = matrix(0,nrow=0,ncol=2)
  )
  if( spatial_method_code %in% c(1,3) ){
    tmb_data$spatial_list = spatial_list
  }else if( spatial_method_code %in% 6 ){
    tmb_data$spatial_list = spatial_list
    #tmb_data$V_zk = V_zk
  }else if( spatial_method_code %in% 7 ){
    tmb_data$nngp_data = nngp_data
    #tmb_data$V_zk = V_zk
  }else if( spatial_method_code %in% 2 ){
    tmb_data$Adj = Adj
  }else if( spatial_method_code %in% 4 ){
    tmb_data$graph_sz = as.matrix(spatial_list[,c('from','to')]) - 1
    tmb_data$dist_s = spatial_list[,c('dist')]
  }else if( spatial_method_code %in% 5 ){
    tmb_data$i_z = spatial_list$i_z - 1
    tmb_data$j_z = spatial_list$j_z - 1
    tmb_data$delta_z2 = spatial_list$delta_z2
  }

  # make params
  # as.numeric(.) ensures class-numeric even for length=0 (when it would be class-logical), so matches output from obj$env$parList()
  tmb_par = list(
    alpha_j = as.numeric(rep(0,ncol(tmb_data$X_ij))),  # Spline coefficients
    gamma_k = as.numeric(rep(0,sum(tmb_data$Sdims))),  # Spline coefficients ... length is sum() <= ncol(tmb_data$Z_ik) ... not equal when using te/ti/t2
    beta_z = as.numeric(ifelse(spacetime_term_ram$param_type==1, 0.01, 1)),  # as.numeric(.) ensures class-numeric even for length=0 (when it would be class-logical), so matches output from obj$env$parList()
    theta_z = as.numeric(ifelse(space_term_ram$param_type==1, 0.01, 1)),
    nu_z = as.numeric(ifelse(time_term_ram$param_type==1, 0.01, 1)),
    log_lambda = as.numeric(rep(0,sum(tmb_data$Sblock))), #Log spline penalization coefficients
    log_sigmaxi_l = as.numeric(rep(0,ncol(tmb_data$W_il))),

    alpha2_j = as.numeric(rep(0,ncol(tmb_data$X2_ij))),  # Spline coefficients
    gamma2_k = as.numeric(rep(0,sum(tmb_data$S2dims))),  # Spline coefficients ... length is sum() <= ncol(tmb_data$Z_ik) ... not equal when using te/ti/t2
    beta2_z = as.numeric(ifelse(delta_spacetime_term_ram$param_type==1, 0.01, 1)),
    theta2_z = as.numeric(ifelse(delta_space_term_ram$param_type==1, 0.01, 1)),
    nu2_z = as.numeric(ifelse(delta_time_term_ram$param_type==1, 0.01, 1)),
    log_lambda2 = as.numeric(rep(0,sum(tmb_data$S2block))), #Log spline penalization coefficients
    log_sigmaxi2_l = as.numeric(rep(0,ncol(tmb_data$W2_il))),

    #log_sigma = rep( 0, sum(distributions$Nsigma_e) ),
    log_sigma = as.numeric(ifelse( is.na(unlist(distributions$sigma_e)), 0, unlist(distributions$sigma_e) )),
    epsilon_stc = array(0, dim=c(n_s, length(times), length(variables))),
    omega_sc = array(0, dim=c(n_s, length(variables))),
    delta_tc = array(0, dim=c(length(times), length(variables))),
    xi_sl = array(0, dim=c(n_s, ncol(tmb_data$W_il))),

    epsilon2_stc = array(0, dim=c(n_s, length(times), length(variables))),
    omega2_sc = array(0, dim=c(n_s, length(variables))),
    delta2_tc = array(0, dim=c(length(times), length(variables))),
    xi2_sl = array(0, dim=c(n_s, ncol(tmb_data$W2_il))),

    eps = numeric(0),
    log_kappa = log(kappa_startvalue),
    ln_H_input = as.numeric(rep(0, n_Hpars))
  )

  # Telescoping
  if( nrow(spacetime_term_ram$output$ram)==0 ){
    tmb_par$epsilon_stc = tmb_par$epsilon_stc[,numeric(0),,drop=FALSE]   # Keep c original length so n_c is detected correctly
  }
  if( nrow(space_term_ram$output$ram)==0 ){
    tmb_par$omega_sc = tmb_par$omega_sc[,numeric(0),drop=FALSE]
  }
  if( nrow(time_term_ram$output$ram)==0 ){
    tmb_par$delta_tc = tmb_par$delta_tc[,numeric(0),drop=FALSE]
  }
  if( nrow(delta_spacetime_term_ram$output$ram)==0 ){
    tmb_par$epsilon2_stc = tmb_par$epsilon2_stc[,numeric(0),,drop=FALSE]   # Keep c original length so n_c is detected correctly
  }
  if( nrow(delta_space_term_ram$output$ram)==0 ){
    tmb_par$omega2_sc = tmb_par$omega2_sc[,numeric(0),drop=FALSE]
  }
  if( nrow(delta_time_term_ram$output$ram)==0 ){
    tmb_par$delta2_tc = tmb_par$delta2_tc[,numeric(0),drop=FALSE]
  }
  if( spatial_method_code %in% 6 ){
    # First value corresponds to `barrier` from `triangle_formula = offset(barrier)`
    tmb_par$triangle_k = c( log(control$barrier_stiffness), rep(0, ncol(V_zk)-1) )
  }

  # Turn off initial conditions ... cutting from model
  #if( control$estimate_delta0==FALSE ){
  #  tmb_par$delta0_c = numeric(0)
  #}

  # Turn of log_kappa when not needed
  #  Now mapping it off, because inclusion for igraph depends on number of edges
  #if( isFALSE(estimate_kappa) ){
  #  tmb_par = tmb_par[-match("log_kappa",names(tmb_par))]
  #}

  # User-supplied parameters
  # Check for obvious issues ... no NAs except in RAMstart
  name_set = c( "ram_space_term_start", "ram_time_term_start", "ram_spacetime_term_start",
                "ram2_space_term_start", "ram2_time_term_start", "ram2_spacetime_term_start")
  if( any(is.na(tmb_data[-match(name_set,names(tmb_data))])) ){
    stop("Check `tmb_data` for NAs")
  }

  # Map off alpha2_j if no delta-models
  tmb_map = list()
  if( all(distributions$components==1) ){
    tmb_map$alpha2_j = factor( rep(NA,length(tmb_par$alpha2_j)) )
  }

  # Map off log_kappa as needed
  if( isFALSE(estimate_kappa) ){
    tmb_map$log_kappa = factor(NA)
  }

  # Map off ln_H_input
  if( isFALSE(estimate_anisotropy) ){
    tmb_map$ln_H_input = factor( c(NA, NA, seq_along(tmb_par$ln_H_input[-c(1:2)])) )
  }

  # Map off offset for triangle_k
  if( spatial_method_code %in% 6 ){
    tmb_map$triangle_k = factor( c(NA, seq_len(ncol(tmb_data$spatial_list$V_zk)-1)) )
  }

  # Map off fixed log_sigma when unlist(distributions$sigma_e) has a value
  tmb_map$log_sigma = factor(ifelse( is.na(unlist(distributions$sigma_e)), seq_along(tmb_par$log_sigma), NA ))

  #
  tmb_random = c("gamma_k","epsilon_stc","omega_sc","delta_tc","xi_sl",
                 "gamma2_k","epsilon2_stc","omega2_sc","delta2_tc","xi2_sl")
  if( isTRUE(control$reml) ){
    tmb_random = union( tmb_random, c("alpha_j","alpha2_j") )
  }
  if( !is.null(control$tmb_random) ){
    tmb_random = control$tmb_random
  }

  # User-specified tmb_par
  if( !is.null(control$tmb_par) ){
    # Check shape but not numeric values, and give informative error
    attr(tmb_par,"check.passed") = attr(control$tmb_par,"check.passed")
    # Compare dimensions by multiplying both by zero
    if( isFALSE(control$suppress_user_warnings) ){
      warning("Supplying `control$tmb_par`:  use carefully as it may crash your terminal")
    }
    if( isTRUE(all.equal(control$tmb_par, tmb_par, tolerance=Inf)) ){
      tmb_par = control$tmb_par
    }else{
      stop("Not using `control$tmb_par` because it has some difference from `tmb_par` built internally")
    }
  }

  # user-specified tmb_map
  if( !is.null(control$tmb_map) ){
    if( isTRUE(all(names(control$tmb_map) %in% names(tmb_par))) ){
      tmb_map = control$tmb_map
    }else{
      stop("Not using `control$tmb_map` because it has some element not present in from `tmb_par`")
    }
  }

  ##############
  # Fit model
  ##############
  if( isFALSE(control$run_model) ){
    out = list(
      tmb_data = tmb_data,
      tmb_par = tmb_par,
      tmb_map = tmb_map,
      tmb_random = tmb_random
    )
    return(out)
  }

  # Debugging speedup ... if walking through manually, save and reload inputs in case of R crash
  if( FALSE ){
    setwd( R'(C:\Users\James.Thorson\Desktop\Temporary (can be deleted!))' )
    save.image( "debugging.RData" )
    #load.image( "debugging.RData" )
  }

  # Optional compression
  if( !is.null(development$tmbad.sparse_hessian_compress) ){
    config(
      tmbad.sparse_hessian_compress = development$tmbad.sparse_hessian_compress,
      DLL = "tinyVAST"
    )
  }

  if( FALSE ){
    setwd(R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\src)')
    dyn.unload(dynlib("tinyVAST"))
    compile("tinyVAST.cpp" , framework = "TMBad" )
    dyn.load(dynlib("tinyVAST"))
  }
  #browser()
  obj = MakeADFun( data = tmb_data,
                   parameters = tmb_par,
                   map = tmb_map,
                   random = tmb_random,
                   DLL = "tinyVAST",
                   profile = control$profile,
                   silent = control$silent )  #
