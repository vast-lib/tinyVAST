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

#bernoulli <-
#function( link = "log" ){
#  linktemp <- substitute(link)
#  if (!is.character(linktemp))
#    linktemp <- deparse(linktemp)
#  okLinks <- c("log", "identity", "cloglog")
#  if (linktemp %in% okLinks)
#    stats <- stats::make.link(linktemp)
#  else if (is.character(link))
#    stats <- stats::make.link(link)
#
#  x <- c(list(family = "bernoulli", link = linktemp), stats)
#  add_to_family(x)
#}

# Copied from sdmTMB
#' @param link1 First linear predictor link
#' @param link2 Second linear predictor link
#' @export
#' @importFrom stats binomial
#' @examples
#' delta_lognormal()
#' @rdname families
delta_lognormal <- function(link1 = "logit", link2 = "log") {
  link1 <- match.arg(link1)
  link2 <- match.arg(link2)
  f1 <- binomial(link = "logit")
  f2 <- lognormal(link = "log")
  structure(
    list(f1, f2, delta = TRUE,
      link = c("logit", "log"),
      family = c("binomial", "lognormal"),
      clean_name = "delta_lognormal(link1 = 'logit', link2 = 'log')"),
    class = "family")
}

# Copied from sdmTMB
#' @param link1 Link for first part of delta/hurdle model.
#' @param link2 Link for second part of delta/hurdle model.
#' @export
#' @importFrom stats Gamma binomial
#' @examples
#' delta_gamma()
#' @rdname families
delta_gamma <- function(link1 = "logit", link2 = "log") {
  link1 <- match.arg(link1)
  link2 <- match.arg(link2)
  f1 <- binomial(link = "logit")
  f2 <- Gamma(link = "log")
  structure(
    list(f1, f2, delta = TRUE, 
      link = c("logit", "log"),
      family = c("binomial", "Gamma"),
      clean_name = "delta_gamma(link1 = 'logit', link2 = 'log')"),
    class = "family")
}

# Copied from sdmTMB
#' @examples
#' delta_poisson_link_gamma()
#' @rdname families
#' @details `delta_poisson_link_gamma()` is the Poisson-link (complementary
#'   log-log) delta model (Thorson 2018).
#' @references
#' *Poisson-link families*:
#'
#' Thorson, J.T. 2018. Three problems with the conventional delta-model for
#' biomass sampling data, and a computationally efficient alternative. Canadian
#' Journal of Fisheries and Aquatic Sciences, 75(9), 1369-1382.
#' \doi{10.1139/cjfas-2017-0266}
#' @export
delta_poisson_link_gamma <- function(link1 = "log", link2 = "log") {
  link1 <- match.arg(link1)
  link2 <- match.arg(link2)
  f1 <- binomial(link = "log")
  f2 <- Gamma(link = "log")
  structure(
    list(f1, f2, delta = TRUE,
      link = c("log", "log"),
      family = c("binomial", "Gamma"), type = "poisson_link_delta",
    clean_name = "delta_poisson_link_gamma(link1 = 'log', link2 = 'log')"),
  class = "family")
}

