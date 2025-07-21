#' Add vertex covariates to a mesh
#'
#' Interpolates covariate values from a data frame to mesh vertices using
#' inverse distance weighting.
#'
#' @param mesh A mesh object from \pkg{fmesher} or \pkg{sdmTMB} (for
#'   [sdmTMB::sdmTMB()] models only).
#' @param data A data frame with coordinate columns and covariate columns, or an
#'   \pkg{sf} object.
#' @param covariates Character vector of covariate column names to interpolate.
#' @param coords Character vector of coordinate column names.
#'   Ignored if data is an \pkg{sf} object.
#' @param power Numeric power parameter for inverse distance weighting
#'  (default `2`; Euclidean squared decay).
#'
#' @return Modified mesh object with a `vertex_covariates` element added and
#' class `vertex_cov` added.
#' @export
#'
#' @examples
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
#'     coords  = c("X", "Y")
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

add_vertex_covariates <-
function( mesh,
          data,
          covariates,
          coords,
          power = 2) {

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

  for (j in seq_along(covariates)) {
    cov <- covariates[j]

    # Remove rows with missing covariate values
    valid_rows <- !is.na(data[[cov]])
    data_coords <- as.matrix(data[valid_rows, coords, drop = FALSE])
    data_values <- data[valid_rows, cov]

    if (length(data_values) == 0) {
      warning("No non-missing values for covariate ", cov, call. = FALSE)
      next
    }

    # Calculate distances and interpolate for each vertex
    for (i in seq_len(nrow(mesh_vertices))) {
      dx <- data_coords[, 1] - mesh_vertices[i, 1]
      dy <- data_coords[, 2] - mesh_vertices[i, 2]
      distances_sq <- dx^2 + dy^2

      # Handle case where vertex coincides exactly with a data point
      if (any(distances_sq == 0)) {
        vertex_covs[i, j] <- data_values[which(distances_sq == 0)[1]]
      } else {
        # Inverse distance weighting
        distances <- sqrt(distances_sq)
        weights <- 1 / (distances^power)
        vertex_covs[i, j] <- sum(weights * data_values) / sum(weights)
      }
    }
  }

  mesh$vertex_covariates <- as.data.frame(vertex_covs)
  class(mesh) <- c("vertex_coords", class(mesh))
  mesh
}

