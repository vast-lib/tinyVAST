#' @title Make a RAM (Reticular Action Model) from a SEM (structural equation model)
#'
#' @description \code{make_sem_ram} converts SEM arrow notation to \code{ram} describing SEM parameters
#'
#' @param sem structural equation model structure, passed to either \code{\link[sem]{specifyModel}}
#'        or \code{\link[sem]{specifyEquations}} and then parsed to control
#'        the set of path coefficients and variance-covariance parameters
#' @param variables A character vector listing the set of variables
#'
#' @inheritParams sem::specifyModel
#'
#' @return An S3-class \code{"sem_ram"} containing:
#' \describe{
#'   \item{\code{model}}{Output from \code{\link[sem]{specifyEquations}} or \code{\link[sem]{specifyModel}}
#'                       that defines paths and parameters}
#'   \item{\code{ram}}{reticular action module (RAM) describing dependencies}
#' }
#'
#' @export
make_sem_ram <-
function( sem,
          variables,
          quiet = FALSE,
          covs = variables ){

  #
  model = tryCatch(
    specifyModel( text=sem, exog.variances=TRUE, endog.variances=TRUE, covs=as.character(covs), quiet=quiet ),
    error = function(e) e
  )
  if( isFALSE("semmod" %in% class(model)) ){
    model = tryCatch(
      specifyEquations( text=sem, exog.variances=TRUE, endog.variances=TRUE, covs=as.character(covs) ),
      error = function(e) e
    )
  }
  if( isFALSE("semmod" %in% class(model)) ){
    stop("Must supply either input for `sem::specifyModel` or `sem::specifyEquations`")
  }

  #vars = sapply( vars, FUN=function(char){gsub("-", "", gsub(" ", "", char))} )
  n.paths = nrow(model)
  par.names = model[, 2]
  startvalues = model[,3]

  # EXCERPT FROM `getAnywhere("sem.semmod")`
  heads = from = to = rep(0, n.paths)
  for (p in 1:n.paths) {
    #path = sem:::parse.path(model[p, 1])
    path = parse_path(model[p, 1])
    heads[p] = abs(path$direction)
    to[p] = path$second
    from[p] = path$first
    if (path$direction == -1) {
      to[p] = path$first
      from[p] = path$second
    }
  }
  missing_vars = setdiff( c(from,to), variables )
  if( length(missing_vars) > 0 ) stop( "Check `build_ram`:", paste0(missing_vars,sep=", ") )

  ram = data.frame(matrix(0, nrow=p, ncol=5))
  pars = na.omit(unique(par.names))
  ram[, 1] = heads
  ram[, 2] = apply(outer(variables, to, "=="), 2, which)
  ram[, 3] = apply(outer(variables, from, "=="), 2, which)
  par.nos = apply(outer(pars, par.names, "=="), 2, which)
  if(length(par.nos) > 0){
    ram[, 4] = unlist(lapply(par.nos, function(x) if (length(x)==0){0}else{x}))
  }
  ram[, 5] = startvalues
  colnames(ram) = c("heads", "to", "from", "parameter", "start")

  #
  out = list( "ram"=ram, "model"=model )
  class(out) = "sem_ram"
  return(out)
}
