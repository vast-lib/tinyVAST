# modified from glmmTMB
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

#' @export
#' @importFrom sdmTMB lognormal
sdmTMB::lognormal

#' @export
#' @importFrom sdmTMB nbinom2
sdmTMB::nbinom2

#' @export
#' @importFrom sdmTMB nbinom1
sdmTMB::nbinom1

#' @export
#' @importFrom sdmTMB tweedie
sdmTMB::tweedie

#' Additional families
#'
#' Additional families compatible with [tinyVAST()].
#'
#' @param link Link.
#' @export
#' @rdname families
#' @name Families
#'
#' @return
#' A list with elements common to standard R family objects including `family`,
#' `link`, `linkfun`, and `linkinv`. Delta/hurdle model families also have
#' elements `delta` (logical) and `type` (standard vs. Poisson-link).
#' @examples
#' delta_lognormal()
delta_lognormal <- function(link1, link2 = "log", type = c("standard", "poisson-link")) {
  type <- match.arg(type)
  if( missing(link1) ){
    if(type=="standard") link1 = "logit"
    if(type=="poisson-link") link1 = "log"
  }
  l1 <- substitute(link1)
  if (!is.character(l1)) l1 <- deparse(l1)
  l2 <- substitute(link2)
  if (!is.character(l2)) l2 <- deparse(l2)
  f1 <- binomial(link = l1)
  f2 <- lognormal(link = l2)
  if (type == "poisson-link") {
    .type <- "poisson_link_delta"
    clean_name <- paste0("delta_lognormal(link1 = '", l1, "', link2 = '", l2, "', type = 'poisson-link')")
  } else {
    .type <- "standard"
    clean_name <- paste0("delta_lognormal(link1 = '", l1, "', link2 = '", l2, "')")
  }
  structure(list(f1, f2, delta = TRUE, link = c("log", "log"),
    family = c("binomial", "lognormal"), type = .type,
    clean_name = clean_name), class = "family")
}

#' @param link1 Link for first part of delta/hurdle model.
#' @param link2 Link for second part of delta/hurdle model.
#' @param type Delta/hurdle family type. `"standard"` for a classic hurdle
#'   model. `"poisson-link"` for a Poisson-link delta model (Thorson 2018).
#' @export
#' @importFrom stats Gamma binomial
#' @examples
#' delta_gamma()
#' @rdname families
#' @references
#' *Poisson-link delta families*:
#'
#' Thorson, J.T. 2018. Three problems with the conventional delta-model for
#' biomass sampling data, and a computationally efficient alternative. Canadian
#' Journal of Fisheries and Aquatic Sciences, 75(9), 1369-1382.
#' \doi{10.1139/cjfas-2017-0266}
delta_gamma <- function(link1, link2 = "log", type = c("standard", "poisson-link")) {
  type <- match.arg(type)
  if( missing(link1) ){
    if(type=="standard") link1 = "logit"
    if(type=="poisson-link") link1 = "log"
  }
  l1 <- substitute(link1)
  if (!is.character(l1)) l1 <- deparse(l1)
  l2 <- substitute(link2)
  if (!is.character(l2)) l2 <- deparse(l2)
  f1 <- binomial(link = l1)
  f2 <- Gamma(link = l2)
  if (type == "poisson-link") {
    .type <- "poisson_link_delta"
    clean_name <- paste0("delta_gamma(link1 = '", l1, "', link2 = '", l2, "', type = 'poisson-link')")
  } else {
    .type <- "standard"
    clean_name <- paste0("delta_gamma(link1 = '", l1, "', link2 = '", l2, "')")
  }
  structure(list(f1, f2, delta = TRUE, link = c(l1, l2),
    type = .type, family = c("binomial", "Gamma"),
    clean_name = clean_name), class = "family")
}

