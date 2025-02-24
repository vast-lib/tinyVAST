.onAttach <- function(libname, pkgname) {
}

rm_wsp <- function (x) {
  # from brms:::rm_wsp()
  # VIA sdmTMB smoothers.R
  out <- gsub("[ \t\r\n]+", "", x, perl = TRUE)
  dim(out) <- dim(x)
  out
}

all_terms <- function (x) {
  # from brms:::all_terms()
  # VIA sdmTMB smoothers.R
  if (!length(x)) {
    return(character(0))
  }
  if (!inherits(x, "terms")) {
    x <- terms(stats::as.formula(x))
  }
  rm_wsp(attr(x, "term.labels"))
}

get_smooth_terms <- function(terms) {
  # from brms:::all_terms()
  # VIA sdmTMB smoothers.R
  x1 <- grep("s\\(", terms)
  x2 <- grep("t2\\(", terms)
  c(x1, x2)
}

remove_s_and_t2 <- function(formula) {
  # FROM sdmTMB smoothers.R
  terms <- stats::terms(formula)
  terms_labels <- attr(terms, "term.labels")
  drop <- grep("s\\(", terms_labels)
  dropt2 <- grep("t2\\(", terms_labels)
  tdrop <- terms_labels[union(drop, dropt2)]
  for (i in seq_along(tdrop)) {
    formula <- stats::update(formula, paste("~ . -", tdrop[i]))
  }
  formula
}

