

#' @title Extract covariance
#'
#' @description Extract the covariance resulting from a specified path structure
#' and estimated parameters
#'
#' @param object Output from \code{\link{tinyVAST}}
#' @param what Which SEM or DSEM term to extract
#' @param pred Extract the term \code{what} for which linear predictor
#' @param n_times The number of times to include when calculating covariance for a DSEM
#'        component, i.e., \code{time_term} or \code{spacetime_term}.  If missing, the default
#'        is to use the one more than the maximum specified lag (e.g., \code{n_times=2}
#'        by default when the maximum \code{lag=1})
#'
#' @details
#' tinyVAST constructs the covariance from specified path structure and estimated parameters
#'
#' @return
#' The covariance matrix among variables
#'
#' @examples
#' # Simulate settings
#' set.seed(101)
#' theta_xy = 0.4
#' n_x = n_y = 10
#' n_c = 5
#' rho = 0.8
#' resid_sd = 0.5
#'
#' # Simulate GMRFs
#' R_s = exp(-theta_xy * abs(outer(1:n_x, 1:n_y, FUN="-")) )
#' R_ss = kronecker(X=R_s, Y=R_s)
#' delta_fs = mvtnorm::rmvnorm(n_c, sigma=R_ss )
#'
#' # Simulate loadings for two factors
#' L_cf = matrix( rnorm(n_c^2), nrow=n_c )
#' L_cf[,3:5] = 0
#' L_cf = L_cf + resid_sd * diag(n_c)
#'
#' # Simulate correlated densities
#' d_cs = L_cf %*% delta_fs
#'
#' # Shape into longform data-frame and add error
#' Data = data.frame( expand.grid(species=1:n_c, x=1:n_x, y=1:n_y),
#'                    "var"="logn", "z"=exp(as.vector(d_cs)) )
#' Data$n = tweedie::rtweedie( n=nrow(Data), mu=Data$z, phi=0.5, power=1.5 )
#'
#' # make mesh
#' mesh = fmesher::fm_mesh_2d( Data[,c('x','y')] )
#'
#' # Specify factor model with two factors and additional independent variance with shared SD
#' sem = "
#'   f1 -> 1, l1
#'   f1 -> 2, l2
#'   f1 -> 3, l3
#'   f1 -> 4, l4
#'   f1 -> 5, l5
#'   f2 -> 2, l6
#'   f2 -> 3, l7
#'   f2 -> 4, l8
#'   f2 -> 5, l9
#'   f1 <-> f1, NA, 1
#'   f2 <-> f2, NA, 1
#'   1 <-> 1, sd, 1
#'   2 <-> 2, sd, 1
#'   3 <-> 3, sd, 1
#'   4 <-> 4, sd, 1
#'   5 <-> 5, sd, 1
#' "
#'
#' # fit model
#' out = tinyVAST( space_term = sem,
#'            data = Data,
#'            formula = n ~ 0 + factor(species),
#'            spatial_domain = mesh,
#'            family = tweedie(),
#'            variables = c( "f1", "f2", 1:n_c ),
#'            space_columns = c("x","y"),
#'            variable_column = "species",
#'            time_column = "time",
#'            distribution_column = "dist" )
#'
#' Extract covariance among species and factors, where
#' estimated covariance is obtained by ignoring factors
#' V = extract_covariance( out, what = "space_term", pred = "one" )
#'
#' @export
extract_covariance <-
function( object,
          what = c("space_term", "time_term", "spacetime_term"),
          pred = c("one", "two"),
          n_times = NULL ){

  pred = match.arg(pred)
  what = match.arg(what)

  # Extract stuff
  variables = object$internal$variables
  if(what=="space_term"){
    if( pred == "one" ){
      parname = "theta_z"
      ram = object$internal$space_term_ram$output$ram
    }else{
      parname = "theta2_z"
      ram = object$internal$delta_space_term_ram$output$ram
    }
    model = data.frame( "direction" = ram$heads, "lag" = 0, "first" = variables[ram$from], "second" = variables[ram$to],
                        "start" = ram$start, "parameter" = ram$parameter)

    if(is.null(n_times)){
      times = 1
    }else{
      times = seq_len(n_times)
    }
    var_names = variables
  }
  if(what=="time_term"){
    if( pred == "one" ){
      parname = "nu_z"
      model = as.data.frame(object$internal$time_term_ram$output$model)
    }else{
      parname = "nu2_z"
      model = as.data.frame(object$internal$delta_time_term_ram$output$model)
    }
    if(is.null(n_times)){
      times = seq_len( max(as.numeric(model$lag)+1) )
    }else{
      times = seq_len(n_times)
    }
    var_names = apply( expand.grid(paste0("_lag",times-1),variables),
                       FUN=function(vec){paste0(vec[2],vec[1])}, MARGIN=1 )
  }
  if(what=="spacetime_term"){
    if( pred == "one" ){
      parname = "beta_z"
      model = as.data.frame(object$internal$spacetime_term_ram$output$model)
    }else{
      parname = "beta2_z"
      model = as.data.frame(object$internal$delta_spacetime_term_ram$output$model)
    }
    if(is.null(n_times)){
      times = seq_len( max(as.numeric(model$lag)+1) )
    }else{
      times = seq_len(n_times)
    }
    var_names = apply( expand.grid(paste0("_lag",times-1),variables),
                       FUN=function(vec){paste0(vec[2],vec[1])}, MARGIN=1 )
  }

  if( nrow(model) > 0 ){
    # Unpack stuff
    if(is.null(object$internal$parlist)){
      object$internal$parlist = object$obj$env$parList()
    }
    # Extract path matrix
    matrices = dsem:::make_matrices(
      beta_p = object$internal$parlist[[parname]],
      model = model,
      times = times,
      variables = variables
    )
    P_kk = matrices$P_kk
    G_kk = matrices$G_kk

    # Calcualte covariance (not precision, in case G is not invertible)
    IminusP_kk = Matrix::Diagonal(n = nrow(P_kk)) - P_kk
    #Q_kk = t(IminusP_kk) %*% solve(t(G_kk) %*% G_kk) %*% IminusP_kk
    V_kk = Matrix::solve(IminusP_kk) %*% Matrix::t(G_kk) %*% G_kk %*% Matrix::solve(Matrix::t(IminusP_kk))
    dimnames(V_kk) = list( var_names, var_names )
  }else{
    V_kk = "selected term not present in model"
  }

  # Return
  return(V_kk)
}
