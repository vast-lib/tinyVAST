#' Add vertex covariates to a mesh
#'
#' Interpolates covariate values from a data frame to mesh vertices using
#' inverse distance weighting (IDW). Uses \pkg{gstat} for exact IDW
#' interpolation by default, with an optional high-performance \pkg{RANN} method
#' for very large datasets.
#'
#' @param mesh A mesh object from \pkg{fmesher} or \pkg{sdmTMB} (for
#'   [sdmTMB::sdmTMB()] models only).
#' @param data A data frame with coordinate columns and covariate columns, or an
#'   \pkg{sf} object.
#' @param covariates Character vector of covariate column names to interpolate.
#' @param coords Character vector of coordinate column names. Ignored if data is
#'   an \pkg{sf} object.
#' @param power Numeric power parameter for inverse distance weighting (default
#'   `2`; Euclidean squared decay).
#' @param method Interpolation method. Options: `"gstat"` (default, exact
#'   inverse distance weighting using gstat package) or `"rann"` (fast k-nearest
#'   neighbours inverse distance weighting using \pkg{RANN} package for very
#'   large datasets).
#' @param k Number of nearest neighbours to use for `"rann"` method (default
#'   `10`). Ignored for `"gstat"` method.
#' @param barrier Optional \pkg{sf} polygon object defining barrier regions. If
#'   provided, adds a logical `barrier` column to `vertex_covariates` where
#'   `TRUE` indicates vertices inside the barrier polygon and `FALSE` indicates
#'   vertices outside. E.g., in the case of modelling fish in the ocean, `TRUE`
#'   represents vertices over land and `FALSE` represents vertices over water.
#'
#' @return Modified mesh object with a `vertex_covariates` element added and
#'   class `vertex_cov` added.
#' @export
#'
#' @examplesIf requireNamespace("RANN", quietly = TRUE)
#' library(sdmTMB)
#' library(sf)
#'
#' # Regular data frame
#' mesh <- fmesher::fm_mesh_2d(pcod[, c("X", "Y")], cutoff = 10)
#' mesh_with_covs <- add_vertex_covariates(
#'   mesh,
#'   data = qcs_grid,
#'   covariates = c("depth_scaled", "depth_scaled2"),
#'   coords = c("X", "Y")
#' )
#' head(mesh_with_covs$vertex_covariates)
#'
#' # Piped version
#' mesh_with_covs <- fmesher::fm_mesh_2d(pcod[, c("X", "Y")], cutoff = 10) |>
#'   add_vertex_covariates(
#'     qcs_grid,
#'     covariates = c("depth_scaled", "depth_scaled2"),
#'     coords = c("X", "Y")
#'   )
#'
#' # With sf objects (coords automatically extracted)
#' pcod_sf <- st_as_sf(pcod, coords = c("X", "Y"))
#' grid_sf <- st_as_sf(qcs_grid, coords = c("X", "Y"))
#' mesh_sf <- fmesher::fm_mesh_2d(pcod_sf, cutoff = 10) |>
#'   add_vertex_covariates(grid_sf, c("depth"))
#'
#' # With sdmTMB mesh (coordinate names and mesh automatically detected)
#' mesh <- make_mesh(pcod, c("X", "Y"), cutoff = 10) |>
#'   add_vertex_covariates(qcs_grid, c("depth"))
#'
#' # Use RANN method for very large datasets (much faster)
#' mesh_fast <- fmesher::fm_mesh_2d(pcod[, c("X", "Y")], cutoff = 10) |>
#'   add_vertex_covariates(
#'     qcs_grid,
#'     covariates = c("depth_scaled", "depth_scaled2"),
#'     coords = c("X", "Y"),
#'     method = "rann",
#'     k = 15
#'   )
add_vertex_covariates <-
  function(mesh,
           data,
           covariates,
           coords,
           power = 2,
           method = c("gstat", "rann"),
           k = 10,
           barrier = NULL) {
    method <- match.arg(method)

    if (inherits(mesh, "sdmTMBmesh")) {
      fm_mesh <- mesh$mesh
      coords <- mesh$xy_cols
    } else if (inherits(mesh, "fm_mesh_2d")) {
      fm_mesh <- mesh
    } else {
      stop("mesh must be a fm_mesh_2d or sdmTMBmesh object", call. = FALSE)
    }

    if (!is.data.frame(data)) {
      stop("data must be a data frame", call. = FALSE)
    }

    # Handle sf objects
    if (inherits(data, "sf")) {
      coords_matrix <- sf::st_coordinates(data)
      data_coords_df <- data.frame(
        X = coords_matrix[, 1],
        Y = coords_matrix[, 2]
      )
      # Combine with attribute data (drop geometry)
      data <- cbind(sf::st_drop_geometry(data), data_coords_df)
      coords <- c("X", "Y")
    }

    if (!all(coords %in% names(data))) {
      stop("Coordinate columns ", paste(coords, collapse = ", "), " not found in data", call. = FALSE)
    }

    if (!all(covariates %in% names(data))) {
      stop("Covariate columns ", paste(covariates, collapse = ", "), " not found in data", call. = FALSE)
    }

    # Extract mesh vertex coordinates (first 2 columns are X, Y)
    mesh_vertices <- fm_mesh$loc[, 1:2, drop = FALSE]

    # Initialize result matrix
    vertex_covs <- matrix(NA, nrow = nrow(mesh_vertices), ncol = length(covariates))
    colnames(vertex_covs) <- covariates

    if (method == "gstat") {
      vertex_covs <- .interpolate_gstat(mesh_vertices, data, covariates, coords, power)
    } else if (method == "rann") {
      vertex_covs <- .interpolate_rann(mesh_vertices, data, covariates, coords, power, k)
    }

    mesh$vertex_covariates <- as.data.frame(vertex_covs)

    # Add barrier detection if barrier polygon is provided
    if (!is.null(barrier)) {
      if (!inherits(barrier, "sf")) {
        stop("barrier must be an sf object", call. = FALSE)
      }

      # Convert mesh vertices to sf points
      vertex_coords <- data.frame(
        X = mesh_vertices[, 1],
        Y = mesh_vertices[, 2]
      )
      vertex_sf <- sf::st_as_sf(vertex_coords, coords = c("X", "Y"))

      # Set CRS to match barrier if barrier has one
      if (!is.na(sf::st_crs(barrier))) {
        sf::st_crs(vertex_sf) <- sf::st_crs(barrier)
      }

      intersected <- sf::st_intersects(vertex_sf, barrier)
      barrier_col <- lengths(intersected) > 0

      mesh$vertex_covariates$barrier <- barrier_col
    }

    class(mesh) <- c("vertex_coords", class(mesh))
    mesh
  }

.prepare_interpolation_data <- function(mesh_vertices, data, covariates, coords) {
  # Initialize result matrix
  vertex_covs <- matrix(NA, nrow = nrow(mesh_vertices), ncol = length(covariates))
  colnames(vertex_covs) <- covariates

  # Validate and prepare data for each covariate
  prepared_data <- list()
  for (j in seq_along(covariates)) {
    cov <- covariates[j]
    valid_rows <- !is.na(data[[cov]])

    if (sum(valid_rows) == 0) {
      warning("No non-missing values for covariate ", cov, call. = FALSE)
      prepared_data[[j]] <- NULL
      next
    }

    data_subset <- data[valid_rows, ]
    prepared_data[[j]] <- list(
      coords = as.matrix(data_subset[, coords, drop = FALSE]),
      values = data_subset[[cov]],
      subset = data_subset
    )
  }

  return(list(vertex_covs = vertex_covs, prepared_data = prepared_data))
}

.interpolate_gstat <- function(mesh_vertices, data, covariates, coords, power) {
  prep <- .prepare_interpolation_data(mesh_vertices, data, covariates, coords)
  vertex_covs <- prep$vertex_covs

  pred_coords <- data.frame(
    x = mesh_vertices[, 1],
    y = mesh_vertices[, 2]
  )
  pred_coords_sf <- sf::st_as_sf(pred_coords, coords = c("x", "y"))

  for (j in seq_along(covariates)) {
    if (is.null(prep$prepared_data[[j]])) next

    cov <- covariates[j]
    data_subset <- prep$prepared_data[[j]]$subset

    data_sf <- sf::st_as_sf(data_subset, coords = coords)

    # Create IDW model
    idw_model <- gstat::gstat(
      formula = as.formula(paste(cov, "~ 1")),
      data = data_sf,
      set = list(idp = power)
    )
    predicted <- predict(idw_model, pred_coords_sf)

    # Extract predicted values
    vertex_covs[, j] <- predicted$var1.pred
  }

  return(vertex_covs)
}

.interpolate_rann <- function(mesh_vertices, data, covariates, coords, power, k) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("Package 'RANN' is required for method='rann'. Please install it.", call. = FALSE)
  }

  prep <- .prepare_interpolation_data(mesh_vertices, data, covariates, coords)
  vertex_covs <- prep$vertex_covs

  for (j in seq_along(covariates)) {
    if (is.null(prep$prepared_data[[j]])) next

    data_coords <- prep$prepared_data[[j]]$coords
    data_values <- prep$prepared_data[[j]]$values

    # Use RANN for fast k-nearest neighbour search
    nn_result <- RANN::nn2(data_coords, mesh_vertices, k = min(k, nrow(data_coords)))

    # Apply IDW using only k nearest neighbors for each vertex
    for (i in seq_len(nrow(mesh_vertices))) {
      neighbor_indices <- nn_result$nn.idx[i, ]
      neighbor_distances <- nn_result$nn.dists[i, ]

      # Handle exact matches
      if (any(neighbor_distances == 0)) {
        vertex_covs[i, j] <- data_values[neighbor_indices[which(neighbor_distances == 0)[1]]]
      } else {
        # IDW with k nearest neighbors
        weights <- 1 / (neighbor_distances^power)
        neighbor_values <- data_values[neighbor_indices]
        vertex_covs[i, j] <- sum(weights * neighbor_values) / sum(weights)
      }
    }
  }

  return(vertex_covs)
}
