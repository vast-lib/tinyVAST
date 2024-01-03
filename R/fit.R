#' @title Fit vector autoregressive spatio-temporal model
#'
#' @description Fits a vector autoregressive spatio-temporal model using
#'  a minimal feature-set and a widely used interface.
#'
#' @inheritParams dsem::make_dsem_ram
#'
#' @param sem Specification for structural equation model structure for
#'        constructing a space-variable interaction.
#'        \code{sem=NULL} disables the space-variable interaction, and
#'        see [make_sem_ram()] for more description.
#' @param dsem Specification for time-series structural equation model structure
#'        including lagged or simultaneous effects for
#'        constructing a space-variable interaction.
#'        \code{dsem=NULL} disables the space-variable interaction, and see
#'        [make_dsem_ram()]  or [make_eof_ram()]
#'        for more description.
#' @param data Data-frame of predictor, response, and offset variables.  Also includes
#'        variables that specify space, time, variables, and the distribution for samples,
#'        as identified by argument \code{data_colnames}.
#' @param data_colnames A list that indicates what columns of `data` are used
#'        to indicate different space, time, variables, and distributions.
#'        Space and variable are then used to interpret argument `sem`,
#'        space, time, and variables are used to interpret argument `dsem`, and
#'        distribution is used to interpret argument `family`.
#' @param formula Formula with response on left-hand-side and predictors on right-hand-side,
#'        parsed by `mgcv` and hence allowing `s(.)` for splines or `offset(.)` for
#'        an offset.
#' @param family a function returning a class \code{family}, including [gaussian()],
#'        [lognormal()], or [tweedie()].  Alternatively, it can be a named list of
#'        these functions, with names that match levels of
#'        \code{data$data_colnames$distribution}.  Finally, the user can specify
#'        [independent_delta()] to specify a delta model.
#' @param delta_options a named list with slots for \code{delta_formula},
#'        \code{delta_sem}, and \code{delta_dsem}.  These follow the same format as
#'        \code{family}, \code{sem}, and \code{dsem}, but specify options for the
#'        second linear predictor of a delta model, and are only used (or estimable)
#'        when [independent_delta()] for some samples.
#' @param spatial_graph Object that represents spatial relationships, either using
#'        _fmesher_ [fm_mesh_2d()] to apply the SPDE method,
#'        _igraph_ [make_empty_graph()] for independent time-series,
#'        _igraph_ [make_graph()] to apply a simultaneous autoregressive (SAR)
#'        process, [sfnetwork_mesh()] for stream networks,
#'        or `NULL` to specify a single site.
#' @param control Output from [tinyVASTcontrol()], used to define user
#'        settings, and see documentation for that function for details.
#' @param times A integer vector listing the set of times in order.
#'        If \code{times=NULL}, then it is filled in as the vector of integers
#'        from the minimum to maximum value of \code{data$data_colnames$time}
#' @param variables A character vector listing the set of variables.
#'        if \code{variables=NULL}, then it is filled in as the unique values
#'        from \code{data$data_colnames$variables}
#'
#' @details
#' `tinyVAST` includes four basic inputs that specify the model structure:
#' * `formula` specifies covariates and splines in a Generalized Additive Model;
#' * `dsem` specifies interactions among variables and over time, constructing
#'   the space-time-variable interaction.
#' * `sem` specifies interactions among variables and over time, constructing
#'   the space-variable interaction.
#' * `spatial_graph` specifies spatial correlations
#'
#' the default `dsem=NULL` turns off all multivariate and temporal indexing, such
#' that `spatial_graph` is then ignored, and the model collapses
#' to a standard model using \code{\link[mgcv]{gam}}.  To specify a univeriate spatial model,
#' the user must specify both `spatial_graph` and `dsem=""`, where the latter
#' is then parsed to include a single exogenous variance for the single variable
#'
#' | \strong{Model type} | \strong{How to specify} |
#' | --- | --- |
#' | Generalized additive model | specify `spatial_graph=NULL` and `dsem=""`, and then use `formula` to specify splines and covariates |
#' | Dynamic structural equation model (including vector autoregressive, dynamic factor analysis, ARIMA, and structural equation models) | specify `spatial_graph=NULL` and use `dsem` to specify interactions among variables and over time |
#' | Univeriate spatial model | specify `spatial_graph` and `dsem=""`, where the latter is then parsed to include a single exogenous variance for the single variable |
#' | Multivariate spatial model | specify `spatial_graph` and use `dsem` (without any lagged effects) to specify spatial interactions |
#' | Vector autoregressive spatio-temporal model | specify `spatial_graph` and use `dsem=""` to specify interactions among variables and over time, where spatio-temporal variables are constructed via the separable interaction of `dsem` and `spatial_graph` |
#'
#' @importFrom igraph as_adjacency_matrix ecount
#' @importFrom sem specifyModel specifyEquations
#' @importFrom corpcor pseudoinverse
#' @importFrom methods is as
#' @importFrom fmesher fm_evaluator fm_mesh_2d fm_fem
#' @importFrom stats .getXlevels gaussian lm model.frame model.matrix
#'   model.offset model.response na.omit nlminb optimHess pnorm rnorm terms
#'   update.formula
#'
#' @examples
#' # Simulate a 2D AR1 spatial process with a cyclic confounder w
#' n_x = n_y = 25
#' n_w = 10
#' R_xx = exp(-0.4 * abs(outer(1:n_x, 1:n_x, FUN="-")) )
#' R_yy = exp(-0.4 * abs(outer(1:n_y, 1:n_y, FUN="-")) )
#' z = mvtnorm::rmvnorm(1, sigma=kronecker(R_xx,R_yy) )
#'
#' # Simulate nuissance parameter z from oscillatory (day-night) process
#' w = sample(1:n_w, replace=TRUE, size=length(z))
#' Data = data.frame( expand.grid(x=1:n_x, y=1:n_y), w=w, z=as.vector(z) + cos(w/n_w*2*pi))
#' Data$n = Data$z + rnorm(nrow(Data), sd=1)
#'
#' # Add columns for multivariate and temporal dimensions
#' Data$var = "n"
#'
#' # make mesh
#' mesh = fmesher::fm_mesh_2d( Data[,c('x','y')], n=100 )
#'
#' # fit model
#' out = fit( data = Data,
#'            formula = n ~ s(w),
#'            spatial_graph = mesh,
#'            control = tinyVASTcontrol(quiet=TRUE, trace=0),
#'            sem = "n <-> n, sd_n" )
#'
#' @useDynLib tinyVAST, .registration = TRUE
#' @export
fit <-
function( data,
          formula,
          sem = NULL,
          dsem = NULL,
          family = gaussian(),
          delta_options = list(delta_formula = ~ 1),
          data_colnames = list("space"=c("x","y"), "variable"="var", "time"="time", "distribution"="dist"),
          times = NULL,
          variables = NULL,
          spatial_graph = NULL,
          control = tinyVASTcontrol(),
          ... ){

  # https://roxygen2.r-lib.org/articles/rd-formatting.html#tables for roxygen formatting
  start_time = Sys.time()

  # General error checks
  if( isFALSE(is(control, "tinyVASTcontrol")) ) stop("`control` must be made by `tinyVASTcontrol()`")
  if( !is.data.frame(data) ) stop("`data` must be a data frame")

  ##############
  # input telescoping
  ##############

  # Telescope family ... comes before adding `data_colnames$distribution` to `data`
  if( inherits(family,"family") ){
    family = list( "obs"=family )
  }

  # Defaults for missing columns of data
  if( !(data_colnames$variable %in% colnames(data)) ){
    data = data.frame(data, matrix(1, nrow=nrow(data), ncol=1, dimnames=list(NULL,data_colnames$variable)))
  }
  if( !(data_colnames$time %in% colnames(data)) ){
    data = data.frame( data, matrix(1, nrow=nrow(data), ncol=1, dimnames=list(NULL,data_colnames$time)) )
  }
  if( !(data_colnames$distribution %in% colnames(data)) ){
    if( length(family)>1 ) stop("Must supply `dist` if using multiple `family` options")
    data = data.frame( data, matrix(names(family)[1], nrow=nrow(data), ncol=1, dimnames=list(NULL,data_colnames$distribution)) )
  }

  # Defaults for times
  if( is.null(dsem) & is.null(delta_options$delta_dsem) ){
    times = numeric(0)
  }
  if(is.null(times)) times = seq( min(data[,data_colnames$time]), max(data[,data_colnames$time]) )

  # Defaults for variables
  if( is.null(dsem) & is.null(sem) & is.null(delta_options$delta_dsem) & is.null(delta_options$delta_sem) ){
    variables = numeric(0)
  }
  if(is.null(variables)) variables = unique( data[,data_colnames$variable] )


  # Turn of t_i and c_i when times and variables are missing, so that delta_k isn't built
  if( length(times) > 0 ){
    t_i = match( data[,data_colnames$time], times )
  }else{ t_i = integer(0) }
  if( length(variables) > 0 ){
    c_i = match( data[,data_colnames$var], variables )
  }else{ c_i = integer(0) }

  ##############
  # DSEM RAM constructor
  ##############

  build_dsem <-
  function( dsem ){
    # (I-Rho)^-1 * Gamma * (I-Rho)^-1
    if( is.null(dsem) ){
      output = list(
        ram = array( 0, dim=c(0,5), dimnames=list(NULL,c("heads","to","from","parameter","start")) ),
        model = array( 0, dim=c(0,8), dimnames=list(NULL,c("","","","","parameter","first","second","direction")) )
      )
    }else if( isTRUE(is.character(dsem)) ){
      output = make_dsem_ram( dsem, times=times, variables=variables, quiet=control$quiet, covs=variables )
    }else if( is(dsem,"dsem_ram") | is(dsem,"eof_ram") ){
      output = dsem
    }else{
      stop("`dsem` must be either `NULL` or a character-string")
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
      if( !all(c(output$model[,'first'],output$model[,'second']) %in% variables) ){
        stop("Some variable in `dsem` is not in `tsdata`")
      }
    }

    # Check for rank-deficient precision from RAM
    ram_gamma = subset( data.frame(output$ram), heads==2 )
    total_variance_h = tapply( as.numeric(ifelse(is.na(ram_gamma$start), 1, ram_gamma$start)),
            INDEX = ram_gamma$from, FUN=\(x)sum(abs(x)) )
    ram_rho = subset( data.frame(output$ram), heads==1 )
    total_effect_h = tapply( as.numeric(ifelse(is.na(ram_rho$start), 1, ram_rho$start)),
            INDEX = ram_rho$from, FUN=\(x)sum(abs(x)) )
    if( any(total_variance_h==0) & control$gmrf_parameterization=="separable" ){
      stop("Must use gmrf_parameterization=`projection` for the dsem RAM supplied")
    }

    # out
    out = list("output"=output, "param_type"=param_type)
    return(out)
  }
  dsem_ram = build_dsem(dsem)
  delta_dsem_ram = build_dsem( delta_options$delta_dsem )

  ##############
  # SEM RAM constructor
  ##############

  build_sem <-
  function( sem ){
    if( is.null(sem) ){
      output = list(
        ram = array( 0, dim=c(0,5), dimnames=list(NULL,c("heads","to","from","parameter","start")) ),
        model = array( 0, dim=c(0,8), dimnames=list(NULL,c("","","","","parameter","first","second","direction")) )
      )
    }else if( isTRUE(is.character(sem)) ){
      output = make_sem_ram( sem, variables=as.character(variables), quiet=control$quiet, covs=as.character(variables) )
    } else {
      stop("`sem` must be either `NULL` or a character-string")
    }

    # Identify arrow-type for each beta_j estimated in RAM
    which_nonzero = which(output$ram[,4]>0)
    param_type = tapply( output$ram[which_nonzero,1],
                         INDEX=output$ram[which_nonzero,4], FUN=max)

    # Check for rank-deficient precision from RAM
    ram2 = subset( data.frame(output$ram), heads==2 )
    total_variance_h = tapply( as.numeric(ifelse( is.na(ram2$start), 1, ram2$start)),
            INDEX = ram2$from, FUN=\(x)sum(abs(x)) )
    ram1 = subset( data.frame(output$ram), heads==1 )
    total_effect_h = tapply( as.numeric(ifelse(is.na(ram1$start), 1, ram1$start)),
            INDEX = ram1$from, FUN=\(x)sum(abs(x)) )
    if( any(total_variance_h==0) & control$gmrf_parameterization=="separable" ){
      stop("Must use options$gmrf_parameterization=`projection` for the sem RAM supplied")
    }

    # out
    out = list("output"=output, "param_type"=param_type)
    return(out)
  }
  sem_ram = build_sem(sem)
  delta_sem_ram = build_sem( delta_options$delta_sem )

  ##############
  # Spatial domain constructor
  ##############

  if( is(spatial_graph,"fm_mesh_2d") ){
    # SPDE
    n_s = spatial_graph$n
    spatial_method_code = 1
    spatial_list = fm_fem( spatial_graph )
    spatial_list = list("M0"=spatial_list$c0, "M1"=spatial_list$g1, "M2"=spatial_list$g2)
    A_is = fm_evaluator( spatial_graph, loc=as.matrix(data[,data_colnames$space]) )$proj$A
    estimate_kappa = TRUE
  }else if( is(spatial_graph,"igraph") ) {
    # SAR
    spatial_method_code = 2
    Adj = as_adjacency_matrix( spatial_graph, sparse=TRUE )
    n_s = nrow(Adj)
    Match = match( data[,data_colnames$space], rownames(Adj) )
    if(any(is.na(Match))) stop("Check `spatial_graph` for SAR")
    A_is = sparseMatrix( i=1:nrow(data), j=Match, x=rep(1,nrow(data)) )
    # Turn off log_kappa if no edges (i.e., unconnected graph)
    estimate_kappa = ifelse( ecount(spatial_graph)>0, TRUE, FALSE )
  }else if( is(spatial_graph,"sfnetwork_mesh") ){      # if( !is.null(sem) )
    # stream network
    spatial_method_code = 4
    n_s = spatial_graph$n
    spatial_list = spatial_graph$table
    A_is = sfnetwork_evaluator( stream = spatial_graph$stream,
                                loc = as.matrix(data[,data_colnames$space]) )
    estimate_kappa = TRUE
  }else{
    # Single-site
    spatial_method_code = 3
    n_s = 1
    A_is = matrix(1, nrow=nrow(data), ncol=1)    # dgCMatrix
    A_is = as(Matrix(A_is),"CsparseMatrix")
    spatial_list = list( "M0" = as(Matrix(1,nrow=1,ncol=1),"CsparseMatrix"),
                         "M1" = as(Matrix(0,nrow=1,ncol=1),"CsparseMatrix"),
                         "M2" = as(Matrix(0,nrow=1,ncol=1),"CsparseMatrix") )
    estimate_kappa = FALSE
  }
  Atriplet = Matrix::mat2triplet(A_is)

  #
  Aepsilon_zz = cbind(Atriplet$i, Atriplet$j, t_i[Atriplet$i], c_i[Atriplet$i])
  which_Arows = which(apply( Aepsilon_zz, MARGIN=1, FUN=\(x) all(!is.na(x)) & any(x>0) ))
  which_Arows = which_Arows[ which(Atriplet$x[which_Arows] > 0) ]
  if( (nrow(dsem_ram$output$ram)==0) & (nrow(delta_dsem_ram$output$ram)==0) ){
    which_Arows = numeric(0)
  }
  Aepsilon_zz = Aepsilon_zz[which_Arows,,drop=FALSE]
  Aepsilon_z = Atriplet$x[which_Arows]

  #
  Aomega_zz = cbind(Atriplet$i, Atriplet$j, c_i[Atriplet$i])
  which_Arows = which(apply( Aomega_zz, MARGIN=1, FUN=\(x) all(!is.na(x)) ))
  which_Arows = which_Arows[ which(Atriplet$x[which_Arows] > 0) ]
  if( (nrow(sem_ram$output$ram)==0) & (nrow(delta_sem_ram$output$ram)==0) ){
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
    S_list = lapply( seq_along(gam_setup$smooth), \(x) gam_setup$smooth[[x]]$S[[1]] )
    S_kk = .bdiag(S_list)         # join S's in sparse matrix
    Sdims = unlist(lapply(S_list,nrow)) # Find dimension of each S
    if(is.null(Sdims)) Sdims = vector(length=0)

    # Get covariates
    which_se = grep( pattern="s(", x=gam_setup$term.names, fixed=TRUE )
    X_ij = gam_setup$X[,setdiff(seq_len(ncol(gam_setup$X)),which_se),drop=FALSE]
    Z_ik = gam_setup$X[,which_se,drop=FALSE]

    #
    out = list( "X_ij"=X_ij, "Z_ik"=Z_ik, "S_kk"=S_kk, "Sdims"=Sdims, "y_i"=y_i,
                "offset_i"=offset_i, "gam_setup"=gam_setup )
    return(out)
  }
  gam_basis = build_gam_basis( formula )
  delta_formula_with_response = update.formula( formula,
                                  paste(".", paste0(as.character(delta_options$delta_formula),collapse="")) )
  delta_gam_basis = build_gam_basis( delta_formula_with_response )

  ##############
  # distribution/link
  ##############

  build_distributions <-
  function( family ){
    e_i = match( data[,data_colnames$distribution], names(family) )

    # Check for errors
    if( (any(is.na(e_i))) ){
      stop("`data[,data_colnames$distribution]` has values that don't match `names(family)`")
    }

    # Construct log_sigma based on family
    pad_length = function(x){if(length(x)==1) c(x,NA) else x}
    remove_last = \(x) x[-length(x)]
    Nsigma_e = sapply( family, FUN=\(x){
                       switch( x$family[length(x$family)],
                         "gaussian" = 1,
                         "tweedie" = 2,
                         "lognormal" = 1,
                         "poisson" = 0
                       )} )
    Edims_ez = cbind( "start"=remove_last(cumsum(c(0,Nsigma_e))), "length"=Nsigma_e )

    #
    # family = list("obs"=gaussian(),"y"=poisson())
    # family = list("obs"=independent_delta(),"y"=independent_delta())
    # family = list("obs"=gaussian(),"y"=independent_delta())
    family_code = t(rbind(sapply( family, FUN=\(x){
                       pad_length(c("gaussian" = 0,
                         "tweedie" = 1,
                         "lognormal" = 2,
                         "poisson" = 3,
                         "bernoulli" = 4 )[x$family])
                       } )))
    link_code = t(rbind(sapply( family, FUN=\(x){
                       pad_length(c("identity" = 0,
                         "log" = 1,
                         "logit" = 2 )[x$link])
                       } )))
    components = apply( family_code, MARGIN=1,
                                  FUN=\(x)sum(!is.na(x)) )
    out = list( "family_code" = cbind(family_code),
                "link_code" = cbind(link_code),
                "components" = components,
                "e_i" = e_i,
                "Nsigma_e" = Nsigma_e,
                "Edims_ez" = Edims_ez )
    return(out)
  }
  distributions = build_distributions( family )

  ##############
  # Build inputs
  # All interactions among features should come here
  ##############

  # make dat
  tmb_data = list(
    spatial_options = c(spatial_method_code, ifelse(control$gmrf_parameterization=="separable",0,1) ),
    y_i = gam_basis$y_i,
    X_ij = gam_basis$X_ij,
    Z_ik = gam_basis$Z_ik,

    X2_ij = delta_gam_basis$X_ij,
    Z2_ik = delta_gam_basis$Z_ik,

    t_i = t_i - 1, # -1 to convert to CPP index
    c_i = c_i - 1, # -1 to convert to CPP index
    offset_i = gam_basis$offset_i,
    family_ez = distributions$family_code,
    link_ez = distributions$link_code,
    components_e = distributions$components,
    e_i = distributions$e_i - 1, # -1 to convert to CPP index
    Edims_ez = distributions$Edims_ez,
    S_kk = gam_basis$S_kk,
    Sdims = gam_basis$Sdims,

    S2_kk = delta_gam_basis$S_kk,
    S2dims = delta_gam_basis$Sdims,

    Aepsilon_zz = Aepsilon_zz - 1,     # Index form, i, s, t
    Aepsilon_z = Aepsilon_z,
    Aomega_zz = Aomega_zz - 1,     # Index form, i, s, t
    Aomega_z = Aomega_z,
    ram_sem = as.matrix(na.omit(sem_ram$output$ram[,1:4])),
    ram_sem_start = as.numeric(sem_ram$output$ram[,5]),
    ram_dsem = as.matrix(na.omit(dsem_ram$output$ram[,1:4])),
    ram_dsem_start = as.numeric(dsem_ram$output$ram[,5]),

    ram2_sem = as.matrix(na.omit(delta_sem_ram$output$ram[,1:4])),
    ram2_sem_start = as.numeric(delta_sem_ram$output$ram[,5]),
    ram2_dsem = as.matrix(na.omit(delta_dsem_ram$output$ram[,1:4])),
    ram2_dsem_start = as.numeric(delta_dsem_ram$output$ram[,5]),
    X2_gj = matrix(0,ncol=ncol(delta_gam_basis$X_ij),nrow=0),
    Z2_gk = matrix(0,ncol=ncol(delta_gam_basis$Z_ik),nrow=0),

    X_gj = matrix(0,ncol=ncol(gam_basis$X_ij),nrow=0),
    Z_gk = matrix(0,ncol=ncol(gam_basis$Z_ik),nrow=0),
    AepsilonG_zz = matrix(0,nrow=0,ncol=4),
    AepsilonG_z = numeric(0),
    AomegaG_zz = matrix(0,nrow=0,ncol=4),
    AomegaG_z = numeric(0),
    t_g = integer(0),
    c_g = integer(0),
    offset_g = integer(0),
    e_g = integer(0),
    W_gz = matrix(0,nrow=0,ncol=2),
    V_gz = matrix(0,nrow=0,ncol=2)
  )
  if( spatial_method_code %in% c(1,3) ){
    tmb_data$spatial_list = spatial_list
  }else if( spatial_method_code %in% 2 ){
    tmb_data$Adj = Adj
  }else if( spatial_method_code %in% 4 ){
    tmb_data$graph_sz = as.matrix(spatial_list[,c('from','to')]) - 1
    tmb_data$dist_s = spatial_list[,c('dist')]
  }

  # make params
  tmb_par = list(
    alpha_j = rep(0,ncol(tmb_data$X_ij)),  # Spline coefficients
    gamma_k = rep(0,ncol(tmb_data$Z_ik)),  # Spline coefficients
    beta_z = as.numeric(ifelse(dsem_ram$param_type==1, 0.01, 1)),  # as.numeric(.) ensures class-numeric even for length=0 (when it would be class-logical), so matches output from obj$env$parList()
    theta_z = as.numeric(ifelse(sem_ram$param_type==1, 0.01, 1)),
    log_lambda = rep(0,length(tmb_data$Sdims)), #Log spline penalization coefficients

    alpha2_j = rep(0,ncol(tmb_data$X2_ij)),  # Spline coefficients
    gamma2_k = rep(0,ncol(tmb_data$Z2_ik)),  # Spline coefficients
    beta2_z = as.numeric(ifelse(delta_dsem_ram$param_type==1, 0.01, 1)),  # as.numeric(.) ensures class-numeric even for length=0 (when it would be class-logical), so matches output from obj$env$parList()
    theta2_z = as.numeric(ifelse(delta_sem_ram$param_type==1, 0.01, 1)),
    log_lambda2 = rep(0,length(tmb_data$S2dims)), #Log spline penalization coefficients

    log_sigma = rep( 0, sum(distributions$Nsigma_e) ),
    delta0_c = rep(0, length(variables)),
    epsilon_stc = array(0, dim=c(n_s, length(times), length(variables))),
    omega_sc = array(0, dim=c(n_s, length(variables))),

    epsilon2_stc = array(0, dim=c(n_s, length(times), length(variables))),
    omega2_sc = array(0, dim=c(n_s, length(variables))),

    eps = numeric(0),
    log_kappa = log(1)
  )

  # Telescoping
  if( nrow(dsem_ram$output$ram)==0 ){
    tmb_par$epsilon_stc = tmb_par$epsilon_stc[,numeric(0),,drop=FALSE]   # Keep c original length so n_c is detected correctly
  }
  if( nrow(sem_ram$output$ram)==0 ){
    tmb_par$omega_sc = tmb_par$omega_sc[,numeric(0),drop=FALSE]
  }
  if( nrow(delta_dsem_ram$output$ram)==0 ){
    tmb_par$epsilon2_stc = tmb_par$epsilon2_stc[,numeric(0),,drop=FALSE]   # Keep c original length so n_c is detected correctly
  }
  if( nrow(delta_sem_ram$output$ram)==0 ){
    tmb_par$omega2_sc = tmb_par$omega2_sc[,numeric(0),drop=FALSE]
  }

  # Turn off initial conditions
  if( control$estimate_delta0==FALSE ){
    tmb_par$delta0_c = numeric(0)
  }

  # Turn of log_kappa when not needed
  #  Now mapping it off, because inclusion for igraph depends on number of edges
  #if( isFALSE(estimate_kappa) ){
  #  tmb_par = tmb_par[-match("log_kappa",names(tmb_par))]
  #}

  # User-supplied parameters
  if( !is.null(control$tmb_par) ){
    # Check shape but not numeric values, and give informative error
    attr(tmb_par,"check.passed") = attr(control$tmb_par,"check.passed")
    if( isTRUE(all.equal(control$tmb_par, tmb_par, tolerance=Inf)) ){
      tmb_par = control$tmb_par
    }else{
      stop("Not using `control$tmb_par` because it has some difference from `tmb_par` built internally")
    }
  }

  # Check for obvious issues ... no NAs except in RAMstart
  if( any(is.na(tmb_data[-match(c("ram_sem_start","ram_dsem_start"),names(tmb_data))])) ){
    stop("Check `tmb_data` for NAs")
  }

  # Map off alpha2_j if no delta-models
  tmb_map = list()
  if( all(distributions$components==1) ){
    tmb_map$alpha2_j = factor( rep(NA,length(tmb_par$alpha2_j)) )
  }
  if( isFALSE(estimate_kappa) ){
    tmb_map$log_kappa = factor(NA)
  }

  ##############
  # Fit model
  ##############

  if( FALSE ){
    setwd(R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\src)')
    dyn.unload(dynlib("tinyVAST"))
    compile("tinyVAST.cpp" , framework = "TMBad" )
    dyn.load(dynlib("tinyVAST"))
  }
  obj = MakeADFun( data = tmb_data,
                   parameters = tmb_par,
                   map = tmb_map,
                   random = c("gamma_k","epsilon_stc","omega_sc","gamma2_k","epsilon2_stc","omega2_sc"),
                   DLL = "tinyVAST",
                   profile = control$profile )  #
  #openmp( ... , DLL="tinyVAST" )
  obj$env$beSilent()
  # L = rep$IminusRho_hh %*% rep$Gamma_hh

  # Optimize
  #start_time = Sys.time()
  opt = list( "par"=obj$par )
  for( i in seq_len(max(0,control$nlminb_loops)) ){
    if( isFALSE(control$quiet) ) message("Running nlminb_loop #", i)
    opt = nlminb( start = opt$par,
                  obj = obj$fn,
                  gr = obj$gr,
                  control = list( eval.max = control$eval.max,
                                  iter.max = control$iter.max,
                                  trace = control$trace ) )
  }
  #Sys.time() - start_time
  #opt

  # Newtonsteps
  for( i in seq_len(max(0,control$newton_loops)) ){
    if( isFALSE(control$quiet) ) message("Running newton_loop #", i)
    g = as.numeric( obj$gr(opt$par) )
    h = optimHess(opt$par, fn=obj$fn, gr=obj$gr)
    opt$par = opt$par - solve(h, g)
    opt$objective = obj$fn(opt$par)
  }

  # Run sdreport
  if( isTRUE(control$getsd) ){
    if( isFALSE(control$quiet) ) message("Running sdreport")
    Hess_fixed = optimHess( par=opt$par, fn=obj$fn, gr=obj$gr )
    sdrep = sdreport( obj, hessian.fixed=Hess_fixed )
  }else{
    Hess_fixed = sdrep = NULL
  }

  # bundle and return output
  internal = list(
    dsem_ram = dsem_ram,                                  # for `add_predictions`
    sem_ram = sem_ram,                                    # for `add_predictions`
    dsem = dsem,
    sem = sem,
    gam_setup = gam_basis$gam_setup,
    delta_dsem_ram = delta_dsem_ram,                                  # for `add_predictions`
    delta_sem_ram = delta_sem_ram,                                    # for `add_predictions`
    delta_dsem = delta_options$delta_dsem,
    delta_sem = delta_options$delta_sem,
    data_colnames = data_colnames,
    delta_gam_setup = delta_gam_basis$gam_setup,
    times = times,
    variables = variables,
    parlist = obj$env$parList(par=obj$env$last.par.best),
    Hess_fixed = Hess_fixed,
    control = control,
    family = family                                       # for `add_predictions`
  )
  out = structure( list(
    formula = formula,
    data = data,
    obj = obj,
    opt = opt,
    rep = obj$report(obj$env$last.par.best),
    sdrep = sdrep,
    tmb_inputs = list(tmb_data=tmb_data, tmb_par=tmb_par, tmb_map=tmb_map),
    call = match.call(),
    spatial_graph = spatial_graph,
    data_colnames = data_colnames,
    run_time = Sys.time() - start_time,
    internal = internal
  ), class="tinyVAST" )
  return(out)
}

#' @title Control parameters for tinyVAST
#'
#' @inheritParams stats::nlminb
#' @inheritParams TMB::MakeADFun
#'
#' @param getsd Boolean indicating whether to call [TMB::sdreport()]
#' @param newton_loops Integer number of newton steps to do after running `nlminb`
#' @param tmb_par list of parameters for starting values, with shape identical to
#'        `fit(...)$internal$parlist`
#'
#' @export
tinyVASTcontrol <-
function( nlminb_loops = 1,
          newton_loops = 0,
          eval.max = 1000,
          iter.max = 1000,
          getsd = TRUE,
          quiet = FALSE,
          trace = 1,
          profile = c(),
          tmb_par = NULL,
          gmrf_parameterization = c("separable","projection"),
          estimate_delta0 = FALSE ){

  gmrf_parameterization = match.arg(gmrf_parameterization)

  # Return
  structure( list(
    nlminb_loops = nlminb_loops,
    newton_loops = newton_loops,
    eval.max = eval.max,
    iter.max = iter.max,
    getsd = getsd,
    quiet = quiet,
    trace = trace,
    profile = profile,
    tmb_par = tmb_par,
    gmrf_parameterization = gmrf_parameterization,
    estimate_delta0 = estimate_delta0
  ), class = "tinyVASTcontrol" )
}

#' @title Print fitted tinyVAST object
#'
#' @description Prints output from fitted tinyVAST model
#'
#' @method print tinyVAST
#' @export
print.tinyVAST <-
function( x,
          ... ){
  print(x[c('call','opt','sdrep','run_time')])
}


#' @title summarize tinyVAST
#'
#' @description summarize parameters from a fitted tinyVAST
#'
#' @details
#' tinyVAST includes an "arrow and lag" notation, which specifies the set of
#' path coefficients and exogenous variance parameters to be estimated. Function \code{fit}
#' then estimates the maximum likelihood value for those coefficients and parameters
#' by maximizing the log-marginal likelihood.
#'
#' However, many users will want to associate individual parameters and standard errors
#' with the path coefficients that were specified using the "arrow and lag" notation.
#' This task is complicated in
#' models where some path coefficients or variance parameters are specified to share a single value a priori,
#' or were assigned a name of NA and hence assumed to have a fixed value a priori (such that
#' these coefficients or parameters have an assigned value but no standard error).
#' The \code{summary} function therefore compiles the MLE for coefficients (including duplicating
#' values for any path coefficients that assigned the same value) and standard error
#' estimates, and outputs those in a table that associates them with the user-supplied path and parameter names.
#' It also outputs the z-score and a p-value arising from a two-sided Wald test (i.e.
#' comparing the estimate divided by standard error against a standard normal distribution).
#'
#' @param object Output from \code{\link{fit}}
#' @param what What component to summarize
#' @param ... Not used
#'
#' @method summary tinyVAST
#' @export
summary.tinyVAST <-
function( object,
          what = c("sem","dsem"),
          ... ){

  #
  what = match.arg(what)

  # SEM component
  if( what=="sem" ){
    model = object$internal$sem_ram$output$ram
    if(nrow(model)>0){
      model$to = as.character(object$internal$variables)[model$to]
      model$from = as.character(object$internal$variables)[model$from]
    }
    ParHat = object$obj$env$parList()

    #
    coefs = data.frame( model, "Estimate"=c(NA,ParHat$theta_z)[ as.numeric(model[,'parameter'])+1 ] ) # parameter=0 outputs NA
    coefs$Estimate = ifelse( is.na(coefs$Estimate), as.numeric(model[,'start']), coefs$Estimate )
    if( "sdrep" %in% names(object) ){
      SE = as.list( object$sdrep, report=FALSE, what="Std. Error")
      coefs = data.frame( coefs, "Std_Error"=c(NA,SE$theta_z)[ as.numeric(model[,'parameter'])+1 ] ) # parameter=0 outputs NA
      coefs = data.frame( coefs, "z_value"=coefs[,'Estimate']/coefs[,'Std_Error'] )
      coefs = data.frame( coefs, "p_value"=pnorm(-abs(coefs[,'z_value'])) * 2 )
    }
  }

  # DSEM component
  if( what=="dsem" ){
    if( is(object$internal$dsem_ram$output,"dsem_ram") ){
      model = object$internal$dsem_ram$output$model
      model = data.frame( heads = model[,'direction'],
                          to = model[,'second'],
                          from = model[,'first'],
                          parameter = model[,'parameter'],
                          start = model[,'start'],
                          lag = model[,'lag'] )
    }else if( is(object$internal$dsem_ram$output,"eof_ram") ){
      model = object$internal$dsem_ram$output$model
      vars = object$internal$dsem_ram$output$variances
      model = data.frame( heads = c( rep(1,nrow(model)), rep(2,nrow(vars)) ),
                          to = c( as.character(model$to), as.character(vars$to) ),
                          from = c( as.character(model$from), as.character(vars$from) ),
                          parameter = c( model$parameter, vars$parameter ),
                          start = NA,
                          lag = NA )
    }else{
      stop("Class for what=`dsem` not implemented for `summary(.)`")
    }
    ParHat = object$obj$env$parList()

    #
    coefs = data.frame( model, "Estimate"=c(NA,ParHat$beta_z)[ as.numeric(model[,'parameter'])+1 ] ) # parameter=0 outputs NA
    coefs$Estimate = ifelse( is.na(coefs$Estimate), as.numeric(model[,'start']), coefs$Estimate )
    if( "sdrep" %in% names(object) ){
      SE = as.list( object$sdrep, report=FALSE, what="Std. Error")
      coefs = data.frame( coefs, "Std_Error"=c(NA,SE$beta_z)[ as.numeric(model[,'parameter'])+1 ] ) # parameter=0 outputs NA
      coefs = data.frame( coefs, "z_value"=coefs[,'Estimate']/coefs[,'Std_Error'] )
      coefs = data.frame( coefs, "p_value"=pnorm(-abs(coefs[,'z_value'])) * 2 )
    }
  }

  return(coefs)
}

#' Calculate residuals
#'
#' @title Calculate deviance or response residuals for tinyVAST
#'
#' @param object Output from [fit()]
#' @param type which type of residuals to compute (only option is `"deviance"` or `"response"` for now)
#' @param ... Note used
#'
#' @method residuals tinyVAST
#' @export
residuals.tinyVAST <-
function( object,
          type = c("deviance","response"),
          ... ){

  # https://stats.stackexchange.com/questions/1432/what-do-the-residuals-in-a-logistic-regression-mean
  # Normal deviance residuals
  if( FALSE ){
    x = rnorm(10)
    y = x + rnorm(10)
    Glm = glm( y ~ 1 + x, family="gaussian")
    mu = predict(Glm,type="response")
    r1 = y - mu
    r1 - resid(Glm)
  }
  # Poisson deviance residuals
  if( FALSE ){
    x = rnorm(10)
    y = rpois(10, exp(x))
    Glm = glm( y ~ 1 + x, family="poisson")
    mu = predict(Glm,type="response")
    # https://stats.stackexchange.com/questions/398098/formula-for-deviance-residuals-for-poisson-model-with-identity-link-function
    r1 = sign(y - mu) * sqrt(2*(y*log((y+1e-10)/mu) - (y-mu)))
    r1 - resid(Glm)
  }
  # Binomial deviance residuals
  if( FALSE ){
    p = 0.5
    y = rbinom(10, prob=p, size=1)
    Glm = glm( y ~ 1, family="binomial")
    mu = predict(Glm, type="response")
    r1 = sign(y - mu) * sqrt(-2*(((1-y)*log(1-mu) + y*log(mu))))
    r1 - resid(Glm)
  }
  # Gamma deviance residuals
  if( FALSE ){
    mu = 1
    cv = 0.8
    y = rgamma( n=10, shape=1/cv^2, scale=mu*cv^2 )
    Glm = glm( y ~ 1, family=Gamma(link='log'))
    mu = predict(Glm, type="response")
    r1 = sign(y - mu) * sqrt(2 * ( (y-mu)/mu - log(y/mu) ))
    r1 - resid(Glm)
  }

  # Poisson: sign(y - mu) * sqrt(2*(ifelse(y==0, 0, y*log(y/mu)) - (y-mu)))
  # Binomial:  -2 * ((1-y)*log(1-mu) + y*log(mu))
  # Gamma: 2 * ( (y-mu)/mu - log(y/mu) )

  # Easy of use
  mu = object$rep$mu
  Y = object$tmb_inputs$tmb_data$Y
  #familycode_j = object$tmb_inputs$data$familycode_j
  report = object$rep

  #
  type = match.arg(type)
  if( type == "deviance" ){
    resid = report$devresid
  }
  if( type == "response" ){
    resid = Y - mu
  }

  return(resid)
}

# Extract the (marginal) log-likelihood of a tinyVAST model
#
# @return object of class \code{logLik} with attributes
#   \item{val}{log-likelihood}
#   \item{df}{number of parameters}
#' @importFrom stats logLik
#' @export
logLik.tinyVAST <- function(object, ...) {
  val = -1 * object$opt$objective
  # Get df including profiled parameters
  df = length( object$opt$par ) +
       sum(names(object$obj$env$last.par) %in% object$internal$control$profile)
  # S3 object "logLik"
  out = structure( val,
             df = df,
             class = "logLik")
  return(out)
}

#' Extract Variance-Covariance Matrix
#'
#' extract the covariance of fixed effects, or both fixed and random effects.
#'
#' @param object output from `fit`
#' @param which whether to extract the covariance among fixed effects, random effects, or both
#' @param ... ignored, for method compatibility
#' @importFrom stats vcov
#' @method vcov tinyVAST
#' @export
vcov.tinyVAST <-
function( object,
          which = c("fixed", "random", "both"),
          ...) {

  which = match.arg(which)

  if( which=="fixed" ){
    V = object$sdrep$cov.fixed
    if(is.null(V)){
      warning("Please re-run `tinyVAST` with `getsd=TRUE`, or confirm that the model is converged")
    }
  }
  if( which=="random" ){
    V = solve(object$obj$env$spHess(random=TRUE))
  }
  if( which=="both" ){
    H = object$sdrep$jointPrecision
    if(is.null(H)){
      warning("Please re-run `tinyVAST` with `getsd=TRUE` and `getJointPrecision=TRUE`, or confirm that the model is converged")
      V = NULL
    }else{
      V = solve(H)
    }
  }

  return( V )
}

