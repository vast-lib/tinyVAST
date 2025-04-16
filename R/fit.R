#' @title Fit vector autoregressive spatio-temporal model
#'
#' @description Fits a vector autoregressive spatio-temporal (VAST) model using
#'  a minimal feature-set and a widely used interface.
#'
#' @param formula Formula with response on left-hand-side and predictors on right-hand-side,
#'        parsed by `mgcv` and hence allowing `s(.)` for splines or `offset(.)` for
#'        an offset.
#' @param data Data-frame of predictor, response, and offset variables.  Also includes
#'        variables that specify space, time, variables, and the distribution for samples,
#'        as identified by arguments `variable_column`, `time_column`, `space_columns`,
#'        and `distribution_column`.
#' @param time_term Specification for time-series structural equation model structure for
#'        constructing a time-variable interaction that defines a time-varying intercept
#'        for each variable (i.e., applies uniformly across space).
#'        \code{time_term=NULL} disables the space-variable interaction;
#'        see [make_dsem_ram()] for notation.
#' @param space_term Specification for structural equation model structure for
#'        constructing a space-variable interaction.
#'        \code{space_term=NULL} disables the space-variable interaction;
#'        see [make_sem_ram()] for notation.
#' @param spacetime_term Specification for time-series structural equation model structure
#'        including lagged or simultaneous effects for
#'        constructing a time-variable interaction, which is then combined in
#'        a separable process with the spatial correlation to form a
#'        space-time-variable interaction (i.e., the interaction occurs locally at each site).
#'        `spacetime_term=NULL` disables the space-variable interaction; see
#'        [make_dsem_ram()]  or [make_eof_ram()].
#' @param family A function returning a class \code{family}, including [gaussian()],
#'        [lognormal()], [tweedie()],  [binomial()],  [Gamma()], [poisson()],
#'        [nbinom1()], or [nbinom2()].
#'        Alternatively, can be a named list of
#'        these functions, with names that match levels of
#'        \code{data$distribution_column} to allow different
#'        families by row of data. Delta model families are possible, and see
#'        \code{\link[tinyVAST:families]{Families}} for delta-model options,
#' @param space_columns A string or character vector that indicates
#'        the column(s) of `data` indicating the location of each sample.
#'        When `spatial_domain` is an `igraph` object, `space_columns` is a string with
#'        with levels matching the names of vertices of that object.
#'        When `spatial_domain` is an `fmesher` or `sfnetwork` object,
#'        space_columns is a character vector indicating columns of `data` with
#'        coordinates for each sample.
#' @param spatial_domain Object that represents spatial relationships, either using
#'        [fmesher::fm_mesh_2d()] to apply the SPDE method,
#'        [igraph::make_empty_graph()] for independent time-series,
#'        [igraph::make_graph()] to apply a simultaneous autoregressive (SAR)
#'        process, [sfnetwork_mesh()] for stream networks,
#'        or `NULL` to specify a single site.  If using `igraph` then the
#'        graph must have vertex names \code{V(graph)$name} that match
#'        levels of \code{data[,'space_columns']}
#' @param time_column A character string indicating the column of `data`
#'        listing the time-interval for each sample, from the set of times
#'        in argument `times`.
#' @param times A integer vector listing the set of times in order.
#'        If `times=NULL`, then it is filled in as the vector of integers
#'        from the minimum to maximum value of `data$time`.  Alternatively,
#'        it could be the minimum value of `data$time` through future years,
#'        such that the model can forecast those future years.
#' @param variable_column A character string indicating the column of `data`
#'        listing the variable for each sample, from the set of times
#'        in argument `variables`.
#' @param variables A character vector listing the set of variables.
#'        if `variables=NULL`, then it is filled in as the unique values
#'        from `data$variable_columns`.
#' @param distribution_column A character string indicating the column of `data`
#'        listing the distribution for each sample, from the set of names
#'        in argument `family`.
#'        if `variables=NULL`, then it is filled in as the unique values
#'        from `data$variables`.
#' @param delta_options a named list with slots for \code{formula},
#'        \code{space_term}, and \code{spacetime_term}. These specify options for the
#'        second linear predictor of a delta model, and are only used (or estimable)
#'        when a \code{\link[tinyVAST:families]{delta family}} is used for some samples.
#' @param spatial_varying a formula specifying spatially varying coefficients.
#' @param weights A numeric vector representing optional likelihood weights for the
#'        data likelihood. Weights do not have to sum to one and are not internally modified.
#'        Thee weights argument needs to be a vector and not a name of the variable in the data frame.
#' @param control Output from [tinyVASTcontrol()], used to define user
#'        settings.
#' @param ... Not used.
#'
#' @details
#' `tinyVAST` includes four basic inputs that specify the model structure:
#' * `formula` specifies covariates and splines in a Generalized Additive Model;
#' * `space_term` specifies interactions among variables and over time, constructing
#'   the space-variable interaction.
#' * `spacetime_term` specifies interactions among variables and over time, constructing
#'   the space-time-variable interaction.
#' * `spatial_domain` specifies spatial correlations
#'
#' the default `spacetime_term=NULL` and `space_term=NULL` turns off all multivariate
#' and temporal indexing, such that `spatial_domain` is then ignored, and the model collapses
#' to a generalized additive model using \code{\link[mgcv]{gam}}.  To specify a univariate spatial model,
#' the user must specify `spatial_domain` and either `space_term=""` or `spacetime_term=""`, where the latter
#' two are then parsed to include a single exogenous variance for the single variable
#'
#' | \strong{Model type} | \strong{How to specify} |
#' | --- | --- |
#' | Generalized additive model | specify `spatial_domain=NULL` `space_term=""` and `spacetime_term=""`, and then use `formula` to specify splines and covariates |
#' | Dynamic structural equation model (including vector autoregressive, dynamic factor analysis, ARIMA, and structural equation models) | specify `spatial_domain=NULL` and use `spacetime_term` to specify interactions among variables and over time |
#' | Univariate spatio-temporal model, or multiple independence spatio-temporal variables | specify `spatial_domain` and `spacetime_term=""`, where the latter is then parsed to include a single exogenous variance for the single variable |
#' | Multivariate spatial model including interactions | specify `spatial_domain` and use `space_term` to specify spatial interactions |
#' | Vector autoregressive spatio-temporal model (i.e., lag-1 interactions among variables) | specify `spatial_domain` and use `spacetime_term=""` to specify interactions among variables and over time, where spatio-temporal variables are constructed via the separable interaction of `spacetime_term` and `spatial_domain` |
#'
#' @importFrom igraph as_adjacency_matrix ecount
#' @importFrom sem specifyModel specifyEquations
#' @importFrom corpcor pseudoinverse
#' @importFrom methods is as
#' @importFrom fmesher fm_evaluator fm_mesh_2d fm_fem
#' @importFrom stats .getXlevels gaussian lm model.frame model.matrix update
#'   model.offset model.response na.omit nlminb optimHess pnorm rnorm terms
#'   update.formula binomial poisson predict
#' @importFrom TMB MakeADFun sdreport
#' @importFrom checkmate assertClass assertDataFrame checkInteger checkNumeric assertNumeric
#' @importFrom Matrix Cholesky solve Matrix diag t
#' @importFrom abind abind
#' @importFrom insight get_response get_data
#' @importFrom cv GetResponse cv
#'
#' @seealso Details section of [make_dsem_ram()] for a summary of the math involved with constructing the DSEM, and \doi{10.1111/2041-210X.14289} for more background on math and inference
#' @seealso \doi{10.48550/arXiv.2401.10193} for more details on how GAM, SEM, and DSEM components are combined from a statistical and software-user perspective
#' @seealso [summary.tinyVAST()] to visualize parameter estimates related to SEM and DSEM model components
#'
#' @return
#' An object (list) of class `tinyVAST`. Elements include:
#' \describe{
#' \item{data}{Data-frame supplied during model fitting}
#' \item{spatial_domain}{the spatial domain supplied during fitting}
#' \item{formula}{the formula specified during model fitting}
#' \item{obj}{The TMB object from \code{\link[TMB]{MakeADFun}}}
#' \item{opt}{The output from \code{\link[stats]{nlminb}}}
#' \item{opt}{The report from \code{obj$report()}}
#' \item{sdrep}{The output from \code{\link[TMB]{sdreport}}}
#' \item{tmb_inputs}{The list of inputs passed to \code{\link[TMB]{MakeADFun}}}
#' \item{call}{A record of the function call}
#' \item{run_time}{Total time to run model}
#' \item{interal}{Objects useful for package function, i.e., all arguments
#'                passed during the call}
#' \item{deviance_explained}{output from \code{\link{deviance_explained}}}
#' }
#'
#' @examples
#' # Simulate a seperable two-dimensional AR1 spatial process
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
#' # Add columns for multivariate and/or temporal dimensions
#' Data$var = "n"
#'
#' # make SPDE mesh for spatial term
#' mesh = fmesher::fm_mesh_2d( Data[,c('x','y')], n=100 )
#'
#' # fit model with cyclic confounder as GAM term
#' out = tinyVAST( data = Data,
#'                 formula = n ~ s(w),
#'                 spatial_domain = mesh,
#'                 space_term = "n <-> n, sd_n" )
#'
#' # Run crossvalidation
#' CV = cv::cv( out )
#' CV
#'
#' @useDynLib tinyVAST, .registration = TRUE
#' @export
tinyVAST <-
function( formula,
          data,
          time_term = NULL,
          space_term = NULL,
          spacetime_term = NULL,
          family = gaussian(),
          space_columns = c("x","y"),
          spatial_domain = NULL,
          time_column = "time",
          times = NULL,
          variable_column = "var",
          variables = NULL,
          distribution_column = "dist",
          delta_options = list(formula = ~ 1),
          spatial_varying = NULL,
          weights = NULL,
          control = tinyVASTcontrol(),
          ... ){

  # https://roxygen2.r-lib.org/articles/rd-formatting.html#tables for roxygen formatting
  start_time = Sys.time()

  # General error checks
  #if( isFALSE(is(control, "tinyVASTcontrol")) ) stop("`control` must be made by `tinyVASTcontrol()`", call. = FALSE)
  #if( !is.data.frame(data) ) stop("`data` must be a data frame", call. = FALSE)
  assertClass(control, "tinyVASTcontrol")
  assertDataFrame(data)
  if(inherits(data,"tbl")) stop("`data` must be a data.frame and cannot be a tibble")

  # Input conflicts
  matched_call = match.call()
  if( isTRUE(as.character(matched_call$family) == "family") ){
    stop("Naming argument `family` as `family` conflicts with function `cv::cv`, please use `family = Family` or other name")
  }

  # Haven't tested for extra levels
  tmpdata = droplevels(data)
  if( !identical(tmpdata,data) ){
    stop("`data` has some factor with extra levels. Please retry after running `data = droplevels(data)` on the input `data`")
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
            INDEX = ram_gamma$from, FUN=\(x)sum(abs(x)) )
    ram_rho = df_ram[df_ram$heads==1,,drop=FALSE]
    total_effect_h = tapply( as.numeric(ifelse(is.na(ram_rho$start), 1, ram_rho$start)),
            INDEX = ram_rho$from, FUN=\(x)sum(abs(x)) )
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
            INDEX = ram2$from, FUN=\(x)sum(abs(x)) )
    ram1 = df_ram[df_ram$heads == 1,,drop=FALSE]
    total_effect_h = tapply( as.numeric(ifelse(is.na(ram1$start), 1, ram1$start)),
            INDEX = ram1$from, FUN=\(x)sum(abs(x)) )
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

  if( is(spatial_domain,"fm_mesh_2d") ){
    # SPDE
    n_s = spatial_domain$n
    spatial_method_code = 1
    spatial_list = fm_fem( spatial_domain )
    spatial_list = list("M0"=spatial_list$c0, "M1"=spatial_list$g1, "M2"=spatial_list$g2)
    A_is = fm_evaluator( spatial_domain, loc=as.matrix(data[,space_columns]) )$proj$A
    estimate_kappa = TRUE
  }else if( is(spatial_domain,"igraph") ) {
    # SAR
    spatial_method_code = 2
    Adj = as_adjacency_matrix( spatial_domain, sparse=TRUE )
    n_s = nrow(Adj)
    Match = match( data[,space_columns], rownames(Adj) )
    if(any(is.na(Match))) stop("Check `spatial_domain` for SAR")
    A_is = sparseMatrix( i=1:nrow(data), j=Match, x=rep(1,nrow(data)) )
    # Turn off log_kappa if no edges (i.e., unconnected graph)
    estimate_kappa = ifelse( ecount(spatial_domain)>0, TRUE, FALSE )
  }else if( is(spatial_domain,"sfnetwork_mesh") ){      # if( !is.null(space_term) )
    # stream network
    spatial_method_code = 4
    n_s = spatial_domain$n
    spatial_list = spatial_domain$table
    A_is = sfnetwork_evaluator( stream = spatial_domain$stream,
                                loc = as.matrix(data[,space_columns]) )
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
  if(n_s>1000) warning("`spatial_domain` has over 1000 components, so the model may be extremely slow")

  #
  Aepsilon_zz = cbind(Atriplet$i, Atriplet$j, t_i[Atriplet$i], c_i[Atriplet$i])
  which_Arows = which(apply( Aepsilon_zz, MARGIN=1, FUN=\(x) all(!is.na(x)) & any(x>0) ))
  which_Arows = which_Arows[ which(Atriplet$x[which_Arows] > 0) ]
  if( (nrow(spacetime_term_ram$output$ram)==0) & (nrow(delta_spacetime_term_ram$output$ram)==0) ){
    which_Arows = numeric(0)
  }
  Aepsilon_zz = Aepsilon_zz[which_Arows,,drop=FALSE]
  Aepsilon_z = Atriplet$x[which_Arows]

  #
  Aomega_zz = cbind(Atriplet$i, Atriplet$j, c_i[Atriplet$i])
  which_Arows = which(apply( Aomega_zz, MARGIN=1, FUN=\(x) all(!is.na(x)) ))
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
    S_list = lapply( seq_along(gam_setup$smooth), \(x) gam_setup$smooth[[x]]$S[[1]] )
    S_kk = Matrix::.bdiag(S_list)       # join S's in sparse matrix
    Sdims = unlist(lapply(S_list,nrow)) # Find dimension of each S
    if(is.null(Sdims)) Sdims = vector(length=0)

    # Get covariates
    not_allowed <- vapply(c("t2(", "te("), \(.x)
      length(grep(.x, x=gam_setup$term.names, fixed=TRUE)) > 0, FUN.VALUE = logical(1L)
    )
    if (any(not_allowed)) {
      stop("Found t2() or te() smoothers. These are not yet implemented.", call. = FALSE)
    }
    which_se = grep( pattern="s(", x=gam_setup$term.names, fixed=TRUE )

    # Extract and add names
    colnames(gam_setup$X) = gam_setup$term.names
    X_ij = gam_setup$X[,setdiff(seq_len(ncol(gam_setup$X)),which_se),drop=FALSE]
    Z_ik = gam_setup$X[,which_se,drop=FALSE]

    #
    out = list( "X_ij"=X_ij, "Z_ik"=Z_ik, "S_kk"=S_kk, "Sdims"=Sdims, "y_i"=y_i,
                "offset_i"=offset_i, "gam_setup"=gam_setup )
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
    remove_last = \(x) x[-length(x)]
    Nsigma_e = sapply( family, FUN=\(x){
                       switch( x$family[length(x$family)],
                         "gaussian" = 1,
                         "tweedie" = 2,
                         "lognormal" = 1,
                         "poisson" = 0,
                         "nbinom2" = 1,
                         "nbinom1" = 1,
                         "binomial" = 0,
                         "bernoulli" = 0,
                         "Gamma" = 1
                       )} )
    Edims_ez = cbind( "start"=remove_last(cumsum(c(0,Nsigma_e))), "length"=Nsigma_e )

    #
    family_code = t(rbind(sapply( family, FUN=\(x){
                       pad_length(c("gaussian" = 0,
                         "tweedie" = 1,
                         "lognormal" = 2,
                         "poisson" = 3,
                         "binomial" = 4,
                         "bernoulli" = 4,
                         "Gamma" = 5,
                         "nbinom1" = 6,
                         "nbinom2" = 7)[x$family])
                       } )))
    link_code = t(rbind(sapply( family, FUN=\(x){
                       pad_length(c("identity" = 0,
                         "log" = 1,
                         "logit" = 2,
                         "cloglog" = 3 )[x$link])
                       } )))
    components = apply( family_code, MARGIN=1,
                         FUN=\(x)sum(x != 99L) )
    poisson_link_delta = sapply( family, FUN=\(x){
                         as.integer(isTRUE(x$type == "poisson_link_delta"))} )
    out = list( "family_code" = cbind(family_code),
                "link_code" = cbind(link_code),
                "components" = components,
                "poisson_link_delta" = poisson_link_delta,
                "e_i" = e_i,
                "Nsigma_e" = Nsigma_e,
                "Edims_ez" = Edims_ez )
    return(out)
  }
  distributions = build_distributions( family )

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
  spatial_options = c(
    spatial_method_code,
    ifelse( control$gmrf_parameterization=="separable", 0, 1),
    ifelse( isFALSE(control$get_rsr), 0, 1)
  )

  # make dat
  tmb_data = list(
    spatial_options = spatial_options,
    y_i = gam_basis$y_i,
    X_ij = gam_basis$X_ij,
    Z_ik = gam_basis$Z_ik,
    W_il = SVC$W_il,
    X2_ij = delta_gam_basis$X_ij,
    Z2_ik = delta_gam_basis$Z_ik,
    W2_il = delta_SVC$W_il,
    t_i = ivector_minus_one(t_i), # -1 to convert to CPP index, and keep as integer-vector
    c_i = ivector_minus_one(c_i), # -1 to convert to CPP index, and keep as integer-vector
    offset_i = gam_basis$offset_i,
    weights_i = weights_i,
    family_ez = distributions$family_code,
    link_ez = distributions$link_code,
    components_e = distributions$components,
    poislink_e = distributions$poisson_link_delta,
    e_i = ivector_minus_one(distributions$e_i), # -1 to convert to CPP index, and keep as integer-vector
    Edims_ez = distributions$Edims_ez,
    S_kk = gam_basis$S_kk,
    Sdims = gam_basis$Sdims,
    S2_kk = delta_gam_basis$S_kk,
    S2dims = delta_gam_basis$Sdims,
    Aepsilon_zz = Aepsilon_zz - 1,     # Index form, i, s, t
    Aepsilon_z = Aepsilon_z,
    Aomega_zz = Aomega_zz - 1,     # Index form, i, s, t
    Aomega_z = Aomega_z,
    A_is = A_is,
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
    beta_z = as.numeric(ifelse(spacetime_term_ram$param_type==1, 0.01, 1)),  # as.numeric(.) ensures class-numeric even for length=0 (when it would be class-logical), so matches output from obj$env$parList()
    theta_z = as.numeric(ifelse(space_term_ram$param_type==1, 0.01, 1)),
    nu_z = as.numeric(ifelse(time_term_ram$param_type==1, 0.01, 1)),
    log_lambda = rep(0,length(tmb_data$Sdims)), #Log spline penalization coefficients
    log_sigmaxi_l = rep(0,ncol(tmb_data$W_il)),

    alpha2_j = rep(0,ncol(tmb_data$X2_ij)),  # Spline coefficients
    gamma2_k = rep(0,ncol(tmb_data$Z2_ik)),  # Spline coefficients
    beta2_z = as.numeric(ifelse(delta_spacetime_term_ram$param_type==1, 0.01, 1)),  # as.numeric(.) ensures class-numeric even for length=0 (when it would be class-logical), so matches output from obj$env$parList()
    theta2_z = as.numeric(ifelse(delta_space_term_ram$param_type==1, 0.01, 1)),
    nu2_z = as.numeric(ifelse(delta_time_term_ram$param_type==1, 0.01, 1)),
    log_lambda2 = rep(0,length(tmb_data$S2dims)), #Log spline penalization coefficients
    log_sigmaxi2_l = rep(0,ncol(tmb_data$W2_il)),

    log_sigma = rep( 0, sum(distributions$Nsigma_e) ),
    epsilon_stc = array(0, dim=c(n_s, length(times), length(variables))),
    omega_sc = array(0, dim=c(n_s, length(variables))),
    delta_tc = array(0, dim=c(length(times), length(variables))),
    xi_sl = array(0, dim=c(n_s, ncol(tmb_data$W_il))),

    epsilon2_stc = array(0, dim=c(n_s, length(times), length(variables))),
    omega2_sc = array(0, dim=c(n_s, length(variables))),
    delta2_tc = array(0, dim=c(length(times), length(variables))),
    xi2_sl = array(0, dim=c(n_s, ncol(tmb_data$W2_il))),

    eps = numeric(0),
    log_kappa = log(1)
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
  if( isFALSE(estimate_kappa) ){
    tmb_map$log_kappa = factor(NA)
  }

  # 
  tmb_random = c("gamma_k","epsilon_stc","omega_sc","delta_tc","xi_sl",
                 "gamma2_k","epsilon2_stc","omega2_sc","delta2_tc","xi2_sl")
  if( isTRUE(control$reml) ){
    tmb_random = union( tmb_random, c("alpha_j","alpha2_j") )
  }

  # User-specified tmb_par
  if( !is.null(control$tmb_par) ){
    # Check shape but not numeric values, and give informative error
    attr(tmb_par,"check.passed") = attr(control$tmb_par,"check.passed")
    # Compare dimensions by multiplying both by zero
    if( isTRUE(all.equal(control$tmb_par*0, tmb_par*0, tolerance=Inf)) ){
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

  if( FALSE ){
    setwd(R'(C:\Users\James.Thorson\Desktop\Git\tinyVAST\src)')
    dyn.unload(dynlib("tinyVAST"))
    compile("tinyVAST.cpp" , framework = "TMBad" )
    dyn.load(dynlib("tinyVAST"))
  }
  obj = MakeADFun( data = tmb_data,
                   parameters = tmb_par,
                   map = tmb_map,
                   random = tmb_random,
                   DLL = "tinyVAST",
                   profile = control$profile,
                   silent = control$silent )  #
  #openmp( ... , DLL="tinyVAST" )
  #obj$env$beSilent()
  # L = rep$IminusRho_hh %*% rep$Gamma_hh

  # Optimize
  #start_time = Sys.time()
  if( isTRUE(control$suppress_nlminb_warnings) ){
    do_nlminb = function( ... ) suppressWarnings(nlminb( ... ))
  }else{
    do_nlminb = function( ... ) nlminb( ... )
  }
  opt = list( "par"=obj$par )
  for( i in seq_len(max(0,control$nlminb_loops)) ){
    if( isTRUE(control$verbose) ) message("Running nlminb_loop #", i)
    opt = do_nlminb( start = opt$par,
                  objective = obj$fn,
                  gradient = obj$gr,
                  control = list( eval.max = control$eval.max,
                                  iter.max = control$iter.max,
                                  trace = control$trace ) )
  }
  #Sys.time() - start_time
  #opt

  # Newtonsteps
  for( i in seq_len(max(0,control$newton_loops)) ){
    if( isTRUE(control$verbose) ) message("Running newton_loop #", i)
    g = as.numeric( obj$gr(opt$par) )
    h = optimHess(opt$par, fn=obj$fn, gr=obj$gr)
    opt$par = opt$par - solve(h, g)
    opt$objective = obj$fn(opt$par)
  }

  # Run sdreport
  if( isTRUE(control$getsd) ){
    if( isTRUE(control$verbose) ) message("Running sdreport")
    Hess_fixed = optimHess( par=opt$par, fn=obj$fn, gr=obj$gr )
    sdrep = sdreport( obj,
                      hessian.fixed = Hess_fixed,
                      getJointPrecision = control$getJointPrecision )
  }else{
    Hess_fixed = sdrep = NULL
  }

  # bundle and return output
  internal = list(
    spacetime_term_ram = spacetime_term_ram,                                  # for `add_predictions`
    space_term_ram = space_term_ram,                                    # for `add_predictions`
    time_term_ram = time_term_ram,                                    # for `add_predictions`
    spacetime_term = spacetime_term,
    space_term = space_term,
    time_term = time_term,
    gam_setup = gam_basis$gam_setup,
    spatial_varying = spatial_varying,
    SVC = SVC,
    delta_spacetime_term_ram = delta_spacetime_term_ram,                                  # for `add_predictions`
    delta_space_term_ram = delta_space_term_ram,                                    # for `add_predictions`
    delta_time_term_ram = delta_time_term_ram,                                    # for `add_predictions`
    delta_spacetime_term = delta_options$spacetime_term,
    delta_space_term = delta_options$space_term,
    delta_time_term = delta_options$time_term,
    space_columns = space_columns,
    time_column = time_column,
    variable_column = variable_column,
    distribution_column = distribution_column,
    delta_gam_setup = delta_gam_basis$gam_setup,
    delta_spatial_varying = delta_options$spatial_varying,
    delta_SVC = delta_SVC,
    times = times,
    variables = variables,
    parlist = obj$env$parList(par=obj$env$last.par.best),
    Hess_fixed = Hess_fixed,
    control = control,
    family = family                                       # for `add_predictions`
  )
  out = list(
    data = data,
    spatial_domain = spatial_domain,
    formula = formula,
    obj = obj,
    opt = opt,
    rep = obj$report(obj$env$last.par.best),
    sdrep = sdrep,
    tmb_inputs = list(tmb_data=tmb_data, tmb_par=tmb_par, tmb_map=tmb_map, tmb_random=tmb_random),
    call = matched_call,
    run_time = Sys.time() - start_time,
    internal = internal
  )

  # Run deviance_explained()
  if( isTRUE(control$calculate_deviance_explained) ){
    out$deviance_explained = deviance_explained( out )
  }else{
    out$deviance_explained = "Not run; use `deviance_explained()` to calculate"
  }

  class(out) = "tinyVAST"
  return(out)
}

#' @title Control parameters for tinyVAST
#'
#' @inheritParams stats::nlminb
#' @inheritParams TMB::MakeADFun
#'
#' @param nlminb_loops Integer number of times to call [stats::nlminb()].
#' @param newton_loops Integer number of Newton steps to do after running
#'   [stats::nlminb()].
#' @param getsd Boolean indicating whether to call [TMB::sdreport()]
#' @param tmb_par list of parameters for starting values, with shape identical
#'   to `tinyVAST(...)$internal$parlist`
#' @param tmb_map input passed to [TMB::MakeADFun] as argument `map`, over-writing
#'   the version `tinyVAST(...)$tmb_inputs$tmb_map` and allowing detailed control
#'   over estimated parameters (advanced feature)
#' @param eval.max Maximum number of evaluations of the objective function
#'   allowed. Passed to `control` in [stats::nlminb()].
#' @param iter.max Maximum number of iterations allowed. Passed to `control` in
#'   [stats::nlminb()].
#' @param verbose Output additional messages about model steps during fitting?
#' @param silent Disable terminal output for inner optimizer?
#' @param trace Parameter values are printed every `trace` iteration
#'   for the outer optimizer. Passed to
#'   `control` in [stats::nlminb()].
#' @param gmrf_parameterization Parameterization to use for the Gaussian Markov 
#'        random field, where the default `separable` constructs a full-rank and
#'        separable precision matrix, and the alternative `projection` constructs
#'        a full-rank and IID precision for variables over time, and then projects
#'        this using the inverse-cholesky of the precision, where this projection
#'        allows for rank-deficient covariance.
#' @param reml Logical: use REML (restricted maximum likelihood) estimation rather than
#'        maximum likelihood? Internally, this adds the fixed effects to the 
#'        list of random effects to integrate over.
#' @param getJointPrecision whether to get the joint precision matrix.  Passed
#'        to \code{\link[TMB]{sdreport}}.
#' @param calculate_deviance_explained whether to calculate proportion of deviance
#'        explained.  See [deviance_explained()]
#' @param run_model whether to run the model of export TMB objects prior to compilation
#'        (useful for debugging)
#' @param suppress_nlminb_warnings whether to suppress uniformative warnings
#'        from \code{nlminb} arising when a function evaluation is NA, which
#'        are then replaced with Inf and avoided during estimation
#' @param get_rsr Experimental option, whether to report restricted spatial
#'        regression (RSR) adjusted estimator for covariate responses
#'
#' @return
#' An object (list) of class `tinyVASTcontrol`, containing either default or
#' updated values supplied by the user for model settings
#'
#' @export
tinyVASTcontrol <-
function( nlminb_loops = 1,
          newton_loops = 0,
          eval.max = 1000,
          iter.max = 1000,
          getsd = TRUE,
          silent = getOption("tinyVAST.silent", TRUE),
          trace = getOption("tinyVAST.trace", 0),
          verbose = getOption("tinyVAST.verbose", FALSE),
          profile = c(),
          tmb_par = NULL,
          tmb_map = NULL,
          gmrf_parameterization = c("separable","projection"),
          #estimate_delta0 = FALSE,
          reml = FALSE,
          getJointPrecision = FALSE,
          calculate_deviance_explained = TRUE,
          run_model = TRUE,
          suppress_nlminb_warnings = TRUE,
          get_rsr = FALSE ){

  gmrf_parameterization = match.arg(gmrf_parameterization)

  # Return
  structure( list(
    nlminb_loops = nlminb_loops,
    newton_loops = newton_loops,
    eval.max = eval.max,
    iter.max = iter.max,
    getsd = getsd,
    silent = silent,
    trace = trace,
    verbose = verbose,
    profile = profile,
    tmb_par = tmb_par,
    tmb_map = tmb_map,
    gmrf_parameterization = gmrf_parameterization,
    #estimate_delta0 = FALSE,
    reml = reml,
    getJointPrecision = getJointPrecision,
    calculate_deviance_explained = calculate_deviance_explained,
    run_model = run_model,
    suppress_nlminb_warnings = suppress_nlminb_warnings,
    get_rsr = get_rsr
  ), class = "tinyVASTcontrol" )
}



