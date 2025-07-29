test_that("add_mesh_covariates works with barrier parameter", {
  library(sf)

  # Create a simple test mesh
  set.seed(123)
  coords <- data.frame(
    X = c(0, 1, 2, 0, 1, 2),
    Y = c(0, 0, 0, 1, 1, 1)
  )
  mesh <- fmesher::fm_mesh_2d(coords, cutoff = 0.5)

  # Create test covariate data
  test_data <- data.frame(
    X = c(0.5, 1.5),
    Y = c(0.5, 0.5),
    depth = c(10, 20)
  )

  # Create a barrier polygon that covers roughly the left half
  barrier_coords <- matrix(c(
    -0.5, -0.5,
    1.0, -0.5,
    1.0, 1.5,
    -0.5, 1.5,
    -0.5, -0.5
  ), ncol = 2, byrow = TRUE)

  barrier_polygon <- st_polygon(list(barrier_coords))
  barrier_sf <- st_sfc(barrier_polygon)
  barrier_sf <- st_sf(data.frame(id = 1), geometry = barrier_sf)

  # Test without barrier
  mesh_no_barrier <- add_mesh_covariates(
    mesh = mesh,
    data = test_data,
    covariates = "depth",
    coords = c("X", "Y")
  )

  # Test with barrier
  mesh_with_barrier <- add_mesh_covariates(
    mesh = mesh,
    data = test_data,
    covariates = "depth",
    coords = c("X", "Y"),
    barrier = barrier_sf
  )

  # Check that barrier column is added
  expect_true("barrier" %in% names(mesh_with_barrier$vertex_covariates))
  expect_false("barrier" %in% names(mesh_no_barrier$vertex_covariates))

  # Check that barrier column is logical
  expect_type(mesh_with_barrier$vertex_covariates$barrier, "logical")

  # Check that we have both TRUE and FALSE values (some vertices inside, some outside)
  barrier_values <- mesh_with_barrier$vertex_covariates$barrier
  expect_true(any(barrier_values)) # Some vertices should be TRUE (inside barrier)
  expect_true(any(!barrier_values)) # Some vertices should be FALSE (outside barrier)

  # Verify that vertices on the left should be TRUE (inside barrier)
  vertex_coords <- mesh$loc[, 1:2]
  left_vertices <- which(vertex_coords[, 1] < 1.0) # X < 1.0
  if (length(left_vertices) > 0) {
    expect_true(all(barrier_values[left_vertices]))
  }
})

test_that("add_mesh_covariates barrier error handling", {
  # Create a simple test mesh
  coords <- data.frame(X = c(0, 1, 0, 1), Y = c(0, 0, 1, 1))
  mesh <- fmesher::fm_mesh_2d(coords, cutoff = 0.5)

  test_data <- data.frame(X = 0.5, Y = 0.5, depth = 10)

  # Test with non-sf barrier object
  expect_error(
    add_mesh_covariates(
      mesh = mesh,
      data = test_data,
      covariates = "depth",
      coords = c("X", "Y"),
      barrier = data.frame(x = 1, y = 1) # Not an sf object
    ),
    "barrier must be an sf object"
  )
})

test_that("add_mesh_covariates barrier with different CRS", {
  library(sf)

  # Create a simple test mesh
  coords <- data.frame(X = c(0, 1, 0, 1), Y = c(0, 0, 1, 1))
  mesh <- fmesher::fm_mesh_2d(coords, cutoff = 0.5)

  test_data <- data.frame(X = 0.5, Y = 0.5, depth = 10)

  # Create barrier with CRS
  barrier_coords <- matrix(c(-0.5, -0.5, 1.5, -0.5, 1.5, 1.5, -0.5, 1.5, -0.5, -0.5),
    ncol = 2, byrow = TRUE
  )
  barrier_polygon <- st_polygon(list(barrier_coords))
  barrier_sf <- st_sfc(barrier_polygon, crs = 4326) # WGS84
  barrier_sf <- st_sf(data.frame(id = 1), geometry = barrier_sf)

  # Should work without error (CRS gets transferred to vertices)
  expect_no_error({
    result <- add_mesh_covariates(
      mesh = mesh,
      data = test_data,
      covariates = "depth",
      coords = c("X", "Y"),
      barrier = barrier_sf
    )
  })

  # Should have barrier column
  result <- add_mesh_covariates(
    mesh = mesh,
    data = test_data,
    covariates = "depth",
    coords = c("X", "Y"),
    barrier = barrier_sf
  )
  expect_true("barrier" %in% names(result$vertex_covariates))
})

test_that("add_mesh_covariates basic covariate interpolation", {
  # Create a simple test mesh
  set.seed(123)
  coords <- data.frame(
    X = c(0, 2, 0, 2),
    Y = c(0, 0, 2, 2)
  )
  mesh <- fmesher::fm_mesh_2d(coords, cutoff = 0.5)

  # Create test data with known values at corners
  test_data <- data.frame(
    X = c(0, 2, 0, 2),
    Y = c(0, 0, 2, 2),
    depth = c(10, 20, 30, 40)
  )

  # Test basic interpolation
  result <- add_mesh_covariates(
    mesh = mesh,
    data = test_data,
    covariates = "depth",
    coords = c("X", "Y")
  )

  # Check that vertex_covariates is created
  expect_true("vertex_covariates" %in% names(result))
  expect_true("depth" %in% names(result$vertex_covariates))
  expect_true(is.data.frame(result$vertex_covariates))

  # Check that all vertices have interpolated values
  expect_equal(nrow(result$vertex_covariates), nrow(mesh$loc))
  expect_true(all(!is.na(result$vertex_covariates$depth)))

  # Check that vertices at exact data points have exact values
  vertex_coords <- mesh$loc[, 1:2]
  for (i in 1:nrow(test_data)) {
    # Find vertices very close to data points
    distances <- sqrt((vertex_coords[, 1] - test_data$X[i])^2 +
      (vertex_coords[, 2] - test_data$Y[i])^2)
    exact_match <- which(distances < 1e-10)
    if (length(exact_match) > 0) {
      expect_equal(result$vertex_covariates$depth[exact_match[1]],
        test_data$depth[i],
        tolerance = 1e-10
      )
    }
  }
})

test_that("add_mesh_covariates multiple covariates", {
  # Create a simple test mesh
  coords <- data.frame(X = c(0, 1, 0, 1), Y = c(0, 0, 1, 1))
  mesh <- fmesher::fm_mesh_2d(coords, cutoff = 0.5)

  # Create test data with multiple covariates
  test_data <- data.frame(
    X = c(0.5, 1.5),
    Y = c(0.5, 0.5),
    depth = c(10, 20),
    temperature = c(5, 15),
    salinity = c(35, 30)
  )

  # Test with multiple covariates
  result <- add_mesh_covariates(
    mesh = mesh,
    data = test_data,
    covariates = c("depth", "temperature", "salinity"),
    coords = c("X", "Y")
  )

  # Check that all covariates are present
  expect_true(all(c("depth", "temperature", "salinity") %in% names(result$vertex_covariates)))

  # Check dimensions
  expect_equal(ncol(result$vertex_covariates), 3)
  expect_equal(nrow(result$vertex_covariates), nrow(mesh$loc))

  # Check that all values are finite
  expect_true(all(is.finite(result$vertex_covariates$depth)))
  expect_true(all(is.finite(result$vertex_covariates$temperature)))
  expect_true(all(is.finite(result$vertex_covariates$salinity)))
})

test_that("add_mesh_covariates with sf objects", {
  library(sf)

  # Create a simple test mesh
  coords <- data.frame(X = c(0, 1, 0, 1), Y = c(0, 0, 1, 1))
  mesh <- fmesher::fm_mesh_2d(coords, cutoff = 0.5)

  # Create sf test data
  test_points <- st_sf(
    depth = c(10, 20),
    temperature = c(5, 15),
    geometry = st_sfc(st_point(c(0.5, 0.5)), st_point(c(1.5, 0.5)))
  )

  # Test with sf input (coords should be ignored)
  result <- add_mesh_covariates(
    mesh = mesh,
    data = test_points,
    covariates = c("depth", "temperature"),
    coords = c("ignored", "also_ignored") # Should be ignored for sf objects
  )

  # Check results
  expect_true("vertex_covariates" %in% names(result))
  expect_true(all(c("depth", "temperature") %in% names(result$vertex_covariates)))
  expect_equal(nrow(result$vertex_covariates), nrow(mesh$loc))
  expect_true(all(is.finite(result$vertex_covariates$depth)))
  expect_true(all(is.finite(result$vertex_covariates$temperature)))
})

test_that("add_mesh_covariates handles missing values", {
  # Create a simple test mesh
  coords <- data.frame(X = c(0, 1, 0, 1), Y = c(0, 0, 1, 1))
  mesh <- fmesher::fm_mesh_2d(coords, cutoff = 0.5)

  # Create test data with missing values
  test_data <- data.frame(
    X = c(0.5, 1.5, 0.5),
    Y = c(0.5, 0.5, 1.5),
    depth = c(10, NA, 30), # One missing value
    temperature = c(5, 15, NA) # Different missing pattern
  )

  # Test that function handles missing values gracefully
  expect_warning(
    result <- add_mesh_covariates(
      mesh = mesh,
      data = test_data,
      covariates = c("depth", "temperature"),
      coords = c("X", "Y")
    ),
    NA # Expect no warnings for now, but could change implementation
  )

  # Check that interpolation still works with remaining values
  expect_true("vertex_covariates" %in% names(result))
  expect_true(all(c("depth", "temperature") %in% names(result$vertex_covariates)))
  expect_true(all(is.finite(result$vertex_covariates$depth)))
  expect_true(all(is.finite(result$vertex_covariates$temperature)))
})

test_that("add_mesh_covariates inverse distance weighting", {
  # Create a test mesh with vertex at (1, 1)
  coords <- data.frame(X = c(0, 2, 1), Y = c(0, 0, 1))
  mesh <- fmesher::fm_mesh_2d(coords, cutoff = 0.5)

  # Create test data - two points equidistant from (1,1)
  test_data <- data.frame(
    X = c(1, 1), # Same X coordinate
    Y = c(0.5, 1.5), # Equidistant from Y=1
    depth = c(10, 30) # Should average to 20 at (1,1)
  )

  # Test with default power (should be 2)
  result <- add_mesh_covariates(
    mesh = mesh,
    data = test_data,
    covariates = "depth",
    coords = c("X", "Y"),
    power = 2
  )

  # Find vertex closest to (1, 1)
  vertex_coords <- mesh$loc[, 1:2]
  distances_to_center <- sqrt((vertex_coords[, 1] - 1)^2 + (vertex_coords[, 2] - 1)^2)
  center_vertex <- which.min(distances_to_center)

  # The interpolated value should be close to the average (20) due to equal distances
  interpolated_value <- result$vertex_covariates$depth[center_vertex]
  expect_true(interpolated_value > 15 && interpolated_value < 25)
})

test_that("add_mesh_covariates error conditions", {
  # Create a simple test mesh
  coords <- data.frame(X = c(0, 1, 0, 1), Y = c(0, 0, 1, 1))
  mesh <- fmesher::fm_mesh_2d(coords, cutoff = 0.5)

  test_data <- data.frame(X = 0.5, Y = 0.5, depth = 10)

  # Test missing coordinate columns
  expect_error(
    add_mesh_covariates(
      mesh = mesh,
      data = test_data,
      covariates = "depth",
      coords = c("missing_X", "missing_Y")
    ),
    "Coordinate columns missing_X, missing_Y not found in data"
  )

  # Test missing covariate columns
  expect_error(
    add_mesh_covariates(
      mesh = mesh,
      data = test_data,
      covariates = "missing_depth",
      coords = c("X", "Y")
    ),
    "Covariate columns missing_depth not found in data"
  )

  # Test non-data.frame input
  expect_error(
    add_mesh_covariates(
      mesh = mesh,
      data = "not_a_dataframe",
      covariates = "depth",
      coords = c("X", "Y")
    ),
    "data must be a data frame"
  )

  # Test invalid mesh type
  expect_error(
    add_mesh_covariates(
      mesh = "not_a_mesh",
      data = test_data,
      covariates = "depth",
      coords = c("X", "Y")
    ),
    "mesh must be a fm_mesh_2d or sdmTMBmesh object"
  )
})

# Tests for interpolation methods gstat and RANN
test_that("add_mesh_covariates method parameter validation", {
  # Create a simple test mesh
  coords <- data.frame(X = c(0, 1, 0, 1), Y = c(0, 0, 1, 1))
  mesh <- fmesher::fm_mesh_2d(coords, cutoff = 0.5)

  test_data <- data.frame(X = 0.5, Y = 0.5, depth = 10)

  # Test valid method arguments
  expect_no_error(
    add_mesh_covariates(
      mesh = mesh,
      data = test_data,
      covariates = "depth",
      coords = c("X", "Y"),
      method = "gstat"
    )
  )

  expect_no_error(
    add_mesh_covariates(
      mesh = mesh,
      data = test_data,
      covariates = "depth",
      coords = c("X", "Y"),
      method = "rann"
    )
  )

  # Test invalid method argument
  expect_error(
    add_mesh_covariates(
      mesh = mesh,
      data = test_data,
      covariates = "depth",
      coords = c("X", "Y"),
      method = "invalid_method"
    ),
    "'arg' should be one of"
  )
})

test_that("add_mesh_covariates gstat method produces correct results", {
  skip_if_not_installed("sp")

  # Create a simple test mesh
  set.seed(123)
  coords <- data.frame(X = c(0, 2, 0, 2), Y = c(0, 0, 2, 2))
  mesh <- fmesher::fm_mesh_2d(coords, cutoff = 0.5)

  # Create test data with known values
  test_data <- data.frame(
    X = c(0, 2, 0, 2),
    Y = c(0, 0, 2, 2),
    depth = c(10, 20, 30, 40)
  )

  # Test gstat method
  result_gstat <- add_mesh_covariates(
    mesh = mesh,
    data = test_data,
    covariates = "depth",
    coords = c("X", "Y"),
    method = "gstat"
  )

  # Check basic structure
  expect_true("vertex_covariates" %in% names(result_gstat))
  expect_true("depth" %in% names(result_gstat$vertex_covariates))
  expect_equal(nrow(result_gstat$vertex_covariates), nrow(mesh$loc))
  expect_true(all(is.finite(result_gstat$vertex_covariates$depth)))

  # Check that interpolated values are within reasonable range
  depth_values <- result_gstat$vertex_covariates$depth
  expect_true(all(depth_values >= min(test_data$depth)))
  expect_true(all(depth_values <= max(test_data$depth)))
})

test_that("add_mesh_covariates RANN method produces correct results", {
  skip_if_not_installed("RANN")

  # Create a simple test mesh
  set.seed(123)
  coords <- data.frame(X = c(0, 2, 0, 2), Y = c(0, 0, 2, 2))
  mesh <- fmesher::fm_mesh_2d(coords, cutoff = 0.5)

  # Create test data with known values
  test_data <- data.frame(
    X = c(0, 2, 0, 2),
    Y = c(0, 0, 2, 2),
    depth = c(10, 20, 30, 40)
  )

  # Test RANN method with different k values
  result_rann_k5 <- add_mesh_covariates(
    mesh = mesh,
    data = test_data,
    covariates = "depth",
    coords = c("X", "Y"),
    method = "rann",
    k = 5
  )

  result_rann_k2 <- add_mesh_covariates(
    mesh = mesh,
    data = test_data,
    covariates = "depth",
    coords = c("X", "Y"),
    method = "rann",
    k = 2
  )

  # Check basic structure for both
  for (result in list(result_rann_k5, result_rann_k2)) {
    expect_true("vertex_covariates" %in% names(result))
    expect_true("depth" %in% names(result$vertex_covariates))
    expect_equal(nrow(result$vertex_covariates), nrow(mesh$loc))
    expect_true(all(is.finite(result$vertex_covariates$depth)))

    # Check that interpolated values are within reasonable range
    depth_values <- result$vertex_covariates$depth
    expect_true(all(depth_values >= min(test_data$depth)))
    expect_true(all(depth_values <= max(test_data$depth)))
  }

  # Results with different k should be similar but not identical
  correlation <- cor(
    result_rann_k5$vertex_covariates$depth,
    result_rann_k2$vertex_covariates$depth
  )
  expect_true(correlation > 0.8) # Should be highly correlated
})

test_that("add_mesh_covariates methods comparison", {
  skip_if_not_installed("RANN")

  # Create a larger test dataset for meaningful comparison
  set.seed(123)
  n_points <- 20
  coords <- data.frame(
    X = runif(n_points, 0, 10),
    Y = runif(n_points, 0, 10)
  )
  mesh <- fmesher::fm_mesh_2d(coords, cutoff = 1)

  # Create test data
  test_data <- data.frame(
    X = runif(50, 0, 10),
    Y = runif(50, 0, 10),
    depth = rnorm(50, 25, 5)
  )

  # Compare gstat and RANN methods
  result_gstat <- add_mesh_covariates(
    mesh = mesh,
    data = test_data,
    covariates = "depth",
    coords = c("X", "Y"),
    method = "gstat"
  )

  result_rann <- add_mesh_covariates(
    mesh = mesh,
    data = test_data,
    covariates = "depth",
    coords = c("X", "Y"),
    method = "rann",
    k = 10
  )

  # Both should produce finite results
  expect_true(all(is.finite(result_gstat$vertex_covariates$depth)))
  expect_true(all(is.finite(result_rann$vertex_covariates$depth)))

  # Results should be reasonably correlated
  correlation <- cor(
    result_gstat$vertex_covariates$depth,
    result_rann$vertex_covariates$depth
  )
  expect_true(correlation > 0.7) # Should be reasonably correlated

  # Both should have similar ranges
  gstat_range <- range(result_gstat$vertex_covariates$depth)
  rann_range <- range(result_rann$vertex_covariates$depth)

  # Ranges should overlap significantly
  expect_true(max(gstat_range[1], rann_range[1]) < min(gstat_range[2], rann_range[2]))
})

test_that("add_mesh_covariates RANN k parameter edge cases", {
  skip_if_not_installed("RANN")

  # Create a simple test mesh
  coords <- data.frame(X = c(0, 1, 0, 1), Y = c(0, 0, 1, 1))
  mesh <- fmesher::fm_mesh_2d(coords, cutoff = 0.5)

  # Create test data with only 2 points
  test_data <- data.frame(
    X = c(0.5, 1.5),
    Y = c(0.5, 0.5),
    depth = c(10, 20)
  )

  # Test with k larger than number of data points
  result <- add_mesh_covariates(
    mesh = mesh,
    data = test_data,
    covariates = "depth",
    coords = c("X", "Y"),
    method = "rann",
    k = 10 # More than the 2 available points
  )

  # Should still work (function limits k to available points)
  expect_true("vertex_covariates" %in% names(result))
  expect_true(all(is.finite(result$vertex_covariates$depth)))

  # Test with k = 1
  result_k1 <- add_mesh_covariates(
    mesh = mesh,
    data = test_data,
    covariates = "depth",
    coords = c("X", "Y"),
    method = "rann",
    k = 1
  )

  expect_true("vertex_covariates" %in% names(result_k1))
  expect_true(all(is.finite(result_k1$vertex_covariates$depth)))
})

test_that("add_mesh_covariates creates triangle_covariates", {
  # Create a simple test mesh
  set.seed(123)
  coords <- data.frame(
    X = c(0, 2, 0, 2),
    Y = c(0, 0, 2, 2)
  )
  mesh <- fmesher::fm_mesh_2d(coords, cutoff = 0.5)

  # Create test data with known values
  test_data <- data.frame(
    X = c(0, 2, 0, 2),
    Y = c(0, 0, 2, 2),
    depth = c(10, 20, 30, 40)
  )

  # Test that triangle_covariates is created
  result <- add_mesh_covariates(
    mesh = mesh,
    data = test_data,
    covariates = "depth",
    coords = c("X", "Y")
  )

  # Check that triangle_covariates exists
  expect_true("triangle_covariates" %in% names(result))
  expect_true("depth" %in% names(result$triangle_covariates))
  expect_true(is.data.frame(result$triangle_covariates))

  # Check that number of rows matches number of triangles
  n_triangles <- nrow(mesh$graph$tv)
  expect_equal(nrow(result$triangle_covariates), n_triangles)

  # Check that all triangle centers have interpolated values
  expect_true(all(!is.na(result$triangle_covariates$depth)))
  expect_true(all(is.finite(result$triangle_covariates$depth)))

  # Check that interpolated values are within reasonable range
  depth_values <- result$triangle_covariates$depth
  expect_true(all(depth_values >= min(test_data$depth)))
  expect_true(all(depth_values <= max(test_data$depth)))
})

test_that("add_mesh_covariates triangle_covariates with multiple covariates", {
  # Create a simple test mesh
  coords <- data.frame(X = c(0, 1, 0, 1), Y = c(0, 0, 1, 1))
  mesh <- fmesher::fm_mesh_2d(coords, cutoff = 0.5)

  # Create test data with multiple covariates
  test_data <- data.frame(
    X = c(0.5, 1.5),
    Y = c(0.5, 0.5),
    depth = c(10, 20),
    temperature = c(5, 15),
    salinity = c(35, 30)
  )

  # Test with multiple covariates
  result <- add_mesh_covariates(
    mesh = mesh,
    data = test_data,
    covariates = c("depth", "temperature", "salinity"),
    coords = c("X", "Y")
  )

  # Check that all covariates are present in triangle_covariates
  expect_true(all(c("depth", "temperature", "salinity") %in% names(result$triangle_covariates)))

  # Check dimensions
  expect_equal(ncol(result$triangle_covariates), 3)
  expect_equal(nrow(result$triangle_covariates), nrow(mesh$graph$tv))

  # Check that all values are finite
  expect_true(all(is.finite(result$triangle_covariates$depth)))
  expect_true(all(is.finite(result$triangle_covariates$temperature)))
  expect_true(all(is.finite(result$triangle_covariates$salinity)))
})

test_that("add_mesh_covariates triangle_covariates with barrier", {
  library(sf)

  # Create a simple test mesh
  set.seed(123)
  coords <- data.frame(
    X = c(0, 1, 2, 0, 1, 2),
    Y = c(0, 0, 0, 1, 1, 1)
  )
  mesh <- fmesher::fm_mesh_2d(coords, cutoff = 0.5)

  # Create test covariate data
  test_data <- data.frame(
    X = c(0.5, 1.5),
    Y = c(0.5, 0.5),
    depth = c(10, 20)
  )

  # Create a barrier polygon that covers roughly the left half
  barrier_coords <- matrix(c(
    -0.5, -0.5,
    1.0, -0.5,
    1.0, 1.5,
    -0.5, 1.5,
    -0.5, -0.5
  ), ncol = 2, byrow = TRUE)

  barrier_polygon <- st_polygon(list(barrier_coords))
  barrier_sf <- st_sfc(barrier_polygon)
  barrier_sf <- st_sf(data.frame(id = 1), geometry = barrier_sf)

  # Test with barrier
  result <- add_mesh_covariates(
    mesh = mesh,
    data = test_data,
    covariates = "depth",
    coords = c("X", "Y"),
    barrier = barrier_sf
  )

  # Check that barrier column is added to triangle_covariates
  expect_true("barrier" %in% names(result$triangle_covariates))

  # Check that barrier column is logical
  expect_type(result$triangle_covariates$barrier, "logical")

  # Check that we have barrier values for all triangles
  expect_equal(length(result$triangle_covariates$barrier), nrow(mesh$graph$tv))

  # Check that barrier_proportion column is added
  expect_true("barrier_proportion" %in% names(result$triangle_covariates))

  # Check that barrier_proportion is numeric and between 0 and 1
  expect_type(result$triangle_covariates$barrier_proportion, "double")
  expect_true(all(result$triangle_covariates$barrier_proportion >= 0))
  expect_true(all(result$triangle_covariates$barrier_proportion <= 1))

  # Check that we have proportion values for all triangles
  expect_equal(length(result$triangle_covariates$barrier_proportion), nrow(mesh$graph$tv))
})

test_that("add_mesh_covariates barrier_proportion calculates correctly", {
  library(sf)

  # Create a simple rectangular mesh
  coords <- data.frame(
    X = c(0, 2, 0, 2),
    Y = c(0, 0, 2, 2)
  )
  mesh <- fmesher::fm_mesh_2d(coords, cutoff = 0.5)

  # Create test covariate data
  test_data <- data.frame(
    X = c(1, 1),
    Y = c(0.5, 1.5),
    depth = c(10, 20)
  )

  # Create a barrier polygon that covers exactly half of the mesh area
  # This should create different intersection proportions for different triangles
  barrier_coords <- matrix(c(
    0, 0,
    1, 0,
    1, 2,
    0, 2,
    0, 0
  ), ncol = 2, byrow = TRUE)

  barrier_polygon <- st_polygon(list(barrier_coords))
  barrier_sf <- st_sfc(barrier_polygon)
  barrier_sf <- st_sf(data.frame(id = 1), geometry = barrier_sf)

  # Test with barrier
  result <- add_mesh_covariates(
    mesh = mesh,
    data = test_data,
    covariates = "depth",
    coords = c("X", "Y"),
    barrier = barrier_sf
  )

  # Check that barrier_proportion values make sense
  proportions <- result$triangle_covariates$barrier_proportion

  # All proportions should be valid (between 0 and 1)
  expect_true(all(proportions >= 0))
  expect_true(all(proportions <= 1))

  # With this barrier design, we should have at least some non-zero proportions
  expect_true(any(proportions > 0))
})
