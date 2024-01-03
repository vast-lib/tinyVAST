# Copied from sdmTMB, modified from glmmTMB
# extra stuff for Effects package, class, etc.
add_to_family <- function(x) {
  # x <- c(x, list(link = link), make.link(link))
  # Effect.default/glm.fit
  if (is.null(x$aic)) {
    x <- c(x, list(aic = function(...) NA_real_))
  }
  if (is.null(x$initialize)) {
    x <- c(x, list(initialize = expression({
      mustart <- y + 0.1
    })))
  }
  if (is.null(x$dev.resids)) {
    # can't return NA, glm.fit is unhappy
    x <- c(x, list(dev.resids = function(y, mu, wt) {
      rep(0, length(y))
    }))
  }
  class(x) <- "family"
  x
}

#' Additional families
#'
#' Additional families compatible with \pkg{tinyVAST}.
#'
#' @param link Link.
#' @export
#' @rdname families
#' @name Families
#'
#' @details Copied from sdmTMB
#'
#' @return
#' A list with elements common to standard R family objects including `family`,
#' `link`, `linkfun`, and `linkinv`.
#'
#' @examples
#' lognormal(link = "log")
lognormal <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("identity", "log", "inverse")
  if (linktemp %in% okLinks)
    stats <- stats::make.link(linktemp)
  else if (is.character(link))
    stats <- stats::make.link(link)
  x <- c(list(family = "lognormal", link = linktemp), stats)
  add_to_family(x)
}

# Copied from sdmTMB
#' @export
#' @examples
#' tweedie(link = "log")
#' @rdname families
tweedie <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("log", "identity")
  if (linktemp %in% okLinks)
    stats <- stats::make.link(linktemp)
  else if (is.character(link))
    stats <- stats::make.link(link)

  x <- c(list(family = "tweedie", link = linktemp), stats)
  add_to_family(x)
}

# Delta model with independent terms
#' @param family2 Delta model second linear predictor family
#' @param link1 First linear predictor link
#' @param link2 Second linear predictor link
#' @export
#' @rdname families
independent_delta <-
function( family2 = c("lognormal"),
          link1 = c("logit"),
          link2 = c("log") ) {

  family2 = match.arg(family2)
  link1 = match.arg(link1)
  link2 = match.arg(link2)

  out = list( family = c("bernoulli",family2),
        link = c(link1,link2) )
  class(out) = "family"
  return(out)
}
