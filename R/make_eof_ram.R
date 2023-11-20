#' @title Make a RAM (Reticular Action Model)
#'
#' @description \code{make_eof_ram} converts SEM arrow notation to \code{ram} describing SEM parameters
#'
#' @param times A character vector listing the set of times in order
#' @param variables A character vector listing the set of variables
#' @param n_eof Number of EOF modes of variability to estimate
#' @param remove_na Boolean indicating whether to remove NA values from RAM (default) or not.
#'            \code{remove_NA=FALSE} might be useful for exploration and diagnostics for
#'            advanced users
#'
#' @return A reticular action module (RAM) describing dependencies
#'
#' @examples
#' # Two EOFs for two variables
#' make_eof_ram( times = 2010:2020, variables = c("pollock","cod"), n_eof=2 )
#'
#' @export
make_eof_ram <-
function( times,
          variables,
          n_eof,
          remove_na = TRUE,
          standard_deviations = "unequal" ){
  # Docs : https://roxygen2.r-lib.org/articles/formatting.html

  ####### Error checks

  if( isFALSE((standard_deviations %in% c("equal","unequal")) | is.numeric(standard_deviations)) ){
    stop("Check `standard_deviations` in `make_eof_ram`")
  }

  ######## Step 1 -- Make model
  EOF_names = paste0("EOF_",seq_len(n_eof))
  L_tz = matrix(NA, nrow=length(times), ncol=n_eof, dimnames=list(times,EOF_names))
  L_tz[lower.tri(L_tz, diag=TRUE)] = 1:sum(lower.tri(L_tz, diag=TRUE))

  #
  model = data.frame( expand.grid("to"=rownames(L_tz),"from"=colnames(L_tz)), "parameter"=as.vector(L_tz) )
  variances = data.frame( "to" = c(EOF_names,as.character(variables)),
                          "from" = c(EOF_names,as.character(variables)),
                          "parameter" = c(rep(0,n_eof),max(model[,3],na.rm=TRUE)+1:length(variables)) )
  if( standard_deviations == "equal" ){
    variances$parameter = ifelse( variances$parameter==0, 0, min(ifelse(variances$parameter==0,NA,variances$parameter),na.rm=TRUE) )
  }else if( is.numeric(standard_deviations) ){
    variances$parameter = ifelse( variances$parameter==0, 0, NA )
  }
  #model = rbind( model, variances )

  ####### Step 2 -- Make RAM

  # Global stuff
  Q_names = expand.grid( "times"=c(EOF_names,times), "variables"=variables )
  ram = NULL  # heads, to, from, parameter, start

  # Loop through paths
  for( from in seq_len(nrow(Q_names)) ){
  for( to in seq_len(nrow(Q_names)) ){
    if( (as.character(Q_names[to,1]) %in% rownames(L_tz)) & (as.character(Q_names[from,1]) %in% colnames(L_tz)) ){
      if( !is.na(L_tz[as.character(Q_names[to,1]),as.character(Q_names[from,1])]) ){
        ram_new = c( 1, to, from, L_tz[as.character(Q_names[to,1]),as.character(Q_names[from,1])], 0.01 )
        ram = rbind( ram, ram_new )
      }
    }
  }}

  # Loop through variances
  for( from in seq_len(nrow(Q_names)) ){
    varnum = ifelse( Q_names[from,'times'] %in% times,
                     variances[ match(Q_names[from,'variables'], variances[,'to']), 'parameter'],
                     variances[ match(Q_names[from,'times'], variances[,'to']), 'parameter'] )
    ram_new = c( 2, from, from, varnum, ifelse(varnum==0, 1, NA) )
    if( is.na(ram_new[4]) ){
      ram_new[4:5] = c(0, standard_deviations)
    }
    ram = rbind( ram, ram_new )
  }
  dimnames(ram) = list(NULL, c('heads','to','from','parameter','start'))

  if( isTRUE(remove_na) ){
    which_keep = which(apply( ram[,1:4], MARGIN=1, FUN=\(x)!any(is.na(x)) ))
    ram = ram[ which_keep, ]
  }

  #
  out = list( "model" = model,
              "ram" = ram,
              "variances" = variances,
              "standard_deviations" = standard_deviations )
  class(out) = "eof_ram"
  return(out)
}
