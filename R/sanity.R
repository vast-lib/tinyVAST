
#' Sanity check of a tinyVAST model
#'
#' @param object Fitted model from [tinyVAST()].
#' @param big_sd_log10 Value to check size of standard errors against. A value
#'   of 2 would indicate that standard errors greater than `10^2` (i.e., 100)
#'   should be flagged.
#' @param gradient_thresh Gradient threshold to issue warning.
#' @param silent Logical: suppress messages? Useful to set to `TRUE` if running
#'   large numbers of models and just interested in returning sanity list
#'   objects.
#'
#' @return An invisible named list of checks.
#' @export
sanity <- 
function( object, 
          big_sd_log10 = 2, 
          gradient_thresh = 0.001, 
          silent = FALSE) {

  assertClass(object, "tinyVAST")
  hessian_ok <- eigen_values_ok <- gradients_ok <- se_magnitude_ok <- FALSE
  nlminb_ok <- FALSE

  simplify_msg <- "Try simplifying the model, adjusting the mesh, or adding priors"

  if (identical(object$opt$convergence, 0L)) {
    msg <- "Non-linear minimizer suggests successful convergence"
    if (!silent) cli_alert_success(msg)
    nlminb_ok <- TRUE
  } else {
    msg <- "Non-linear minimizer did not converge: do not trust this model!"
    if (!silent) cli::cli_alert_danger(msg)
    if (!silent) cli::cli_alert_info(simplify_msg)
    cat("\n")
  }

  if( is.null(object$sdrep$pdHess) ){
    msg <- "Non-positive-definite Hessian matrix: model may not have converged"
  }else if (isFALSE(object$sdrep$pdHess)) {
    msg <- "Hessian matrix not detected: consider re-running to check if positive definite"
    if (!silent) cli::cli_alert_danger(msg)
    if (!silent) cli::cli_alert_info(simplify_msg)
    cat("\n")
  } else {
    msg <- "Hessian matrix is positive definite"
    if (!silent) cli::cli_alert_success(msg)
    hessian_ok <- TRUE
  }

  if (isTRUE(object$convergence_diagnostics$bad_eig)) {
    msg <- "Extreme or very small eigenvalues detected: model may not have converged"
    if (!silent) cli::cli_alert_danger(msg)
    if (!silent) cli::cli_alert_info(simplify_msg)
    cat("\n")
  } else {
    msg <- "No extreme or very small eigenvalues detected"
    if (!silent) cli::cli_alert_success(msg)
    eigen_values_ok <- TRUE
  }

  g <- object$internal$convergence_diagnostics$final_grads
  np <- names(object$obj$par)
  for (i in seq_along(g)) {
    if (abs(g[i]) > gradient_thresh) {
      if (!silent) {
        cli::cli_alert_danger(c(
          "`", np[i],
          paste0("` gradient > ", gradient_thresh)
        ))
      }
      msg <- "standardize covariates, and/or simplify the model"
      if (!silent) {
        cli::cli_alert_info(msg)
        cat("\n")
      }
    }
  }

  if (all(abs(g) <= gradient_thresh)) {
    msg <- "No gradients with respect to fixed effects are >= "
    if (!silent) cli::cli_alert_success(paste0(msg, gradient_thresh))
    gradients_ok <- TRUE
  }

  obj <- object$tmb_obj
  random <- unique(names(obj$env$par[obj$env$random]))

  pl <- as.list(object$sdrep, "Estimate")
  ple <- as.list(object$sdrep, "Std. Error")
  pars <- names(object$obj$par)
  pl <- pl[pars]
  ple <- ple[pars]
  np <- names(ple)
  se_na_ok <- TRUE
  for (i in seq_along(ple)) {
    if (any(is.na(ple[i]))) {
      if (!silent) cli::cli_alert_danger(c("`", np[i], paste0("` standard error is NA")))
      par_message(np[i])
      if (!silent) {
        cli::cli_alert_info(simplify_msg)
        cat("\n")
      }
      se_na_ok <- FALSE
    }
  }
  if (se_na_ok) {
    msg <- "No fixed-effect standard errors are NA"
    if (!silent) cli::cli_alert_success(msg)
  }

  est <- as.list(object$sd_report, "Estimate")
  se <- as.list(object$sd_report, "Std. Error")
  fixed <- !(names(est) %in% random)
  est <- est[fixed]
  se <- se[fixed]

  too_big <- function(se) {
    if (any(!is.na(se))) {
      se_max <- max(se, na.rm = TRUE)
      if (any(log10(abs(se_max)) > big_sd_log10)) {
        return(TRUE)
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  }

  se_big <- lapply(se, too_big)

  for (i in seq_along(se_big)) {
    if (isTRUE(se_big[[i]])) {
      msg <- paste0("` standard error may be large")
      if (!silent) cli::cli_alert_danger(c("`", names(se_big)[i], msg))
      par_message(names(se_big)[i])
      if (!silent) {
        cli::cli_alert_info(simplify_msg)
        cat("\n")
      }
    }
  }

  if (all(unlist(lapply(se_big, is.null)))) {
    msg <- "No standard errors look unreasonably large"
    if (!silent) cli::cli_alert_success(msg)
    se_magnitude_ok <- TRUE
  }

#  # tidying objects with different names -- fixed won't have model or group_name
#  b <- tidy(object, conf.int = TRUE, silent = TRUE)
#  br <- tidy(object, "ran_pars", conf.int = TRUE, silent = TRUE)
#  b[, names(br)[!names(br) %in% names(b)]] <- NA
#  b <- rbind(b, br)
#
#  if (isTRUE(object$family$delta)) {
#    b2 <- tidy(object, conf.int = TRUE, model = 2, silent = TRUE)
#    br2 <- tidy(object, "ran_pars", conf.int = TRUE, model = 2, silent = TRUE)
#    b2[, names(br2)[!names(br2) %in% names(b2)]] <- NA
#    b2 <- rbind(b2, br2)
#
#    # also check for missing names -- will occur when one model has no random effects
#    b[, names(b2)[!names(b2) %in% names(b)]] <- NA
#    b2[, names(b)[!names(b) %in% names(b2)]] <- NA
#    b <- rbind(b, b2)
#  }
#  s <- grep("sigma", b$term)
#  sigmas_ok <- TRUE
#  if (length(s)) {
#    for (i in s) {
#      if (b$estimate[i] < 0.01) {
#        msg <- "` is smaller than 0.01"
#        if (!silent) cli::cli_alert_danger(c("`", b$term[i], msg))
#        par_message(b$term[i])
#        msg <- "Consider omitting this part of the model"
#        if (!silent) {
#          cli::cli_alert_info(msg)
#          cat("\n")
#        }
#        sigmas_ok <- FALSE
#      }
#    }
#  }
#  if (sigmas_ok) {
#    msg <- "No sigma parameters are < 0.01"
#    if (!silent) cli::cli_alert_success(msg)
#  }
#
#  if (length(s)) {
#    for (i in s) {
#      if (b$estimate[i] > 100) {
#        msg <- "` is larger than 100"
#        if (!silent) cli::cli_alert_danger(c("`", b$term[i], msg))
#        par_message(b$term[i])
#        msg <- "Consider simplifying the model or adding priors"
#        if (!silent) {
#          cli::cli_alert_info(msg)
#          cat("\n")
#        }
#        sigmas_ok <- FALSE
#      }
#    }
#  }
#  if (sigmas_ok) {
#    msg <- "No sigma parameters are > 100"
#    if (!silent) cli::cli_alert_success(msg)
#  }

#  r1 <- diff(range(object$data[[object$spde$xy_cols[1]]]))
#  r2 <- diff(range(object$data[[object$spde$xy_cols[2]]]))
#  r <- max(r1, r2)
#  range_ok <- TRUE
#  if ("range" %in% b$term) {
#    if (max(b$estimate[b$term == "range"]) > r * 1.5) {
#      msg <- "A `range` parameter looks fairly large (> 1.5 the greatest distance in data)"
#      if (!silent) {
#        cli::cli_alert_danger(msg)
#        cli::cli_alert_info(simplify_msg)
#        cli::cli_alert_info("Also make sure your spatial coordinates are not too big or small (e.g., work in UTM km instead of UTM m)", wrap = TRUE)
#        cat("\n")
#      }
#      range_ok <- FALSE
#    } else {
#      nr <- length(grep("range", b$term))
#      if (nr == 1L) msg <- "Range parameter doesn't look unreasonably large"
#      if (nr > 1L) msg <- "Range parameters don't look unreasonably large"
#      if (!silent) cli::cli_alert_success(msg)
#    }
#  }

  ret <- named_list(
    hessian_ok, eigen_values_ok, nlminb_ok, 
    gradients_ok, se_magnitude_ok, se_na_ok 
  )
  all_ok <- all(unlist(ret))
  ret <- c(ret, all_ok = all_ok)
  invisible(ret)
}
