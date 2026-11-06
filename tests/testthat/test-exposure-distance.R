test_that("exposure mappings return the expected neighborhood treatments", {
  Di <- c(0, 1, 1, 0)
  neighbors <- list(c(2, 3), c(1, 4), c(1, 2, 4), integer())

  expect_equal(rddinterf:::one_treated(Di, neighbors), c(1, 0, 1, 0))
  expect_equal(rddinterf:::sum_treated(Di, neighbors), c(2, 0, 1, 0))
  expect_equal(rddinterf:::frac_treated(Di, neighbors), c(1, 0, 1 / 3, NA))
  expect_equal(rddinterf:::majority(Di, neighbors), c(1, 0, 0, NA))
})

test_that("distance calculations retain their reference values", {
  Xi <- c(-0.8, 0.2, -0.1, 0.7)
  Xj <- list(c(-0.5, -0.2), c(0.3, -0.4), c(0.1, 0.6), c(-0.7, 0.2))

  cases <- list(
    one.treated = list(contrast = c(1, 1, 0, 0), expected = c(
      0.824621125123532, 0.360555127546399, 0.616441400296898,
      0.728010988928052
    )),
    sum.treated = list(contrast = c(1, 2, 0, 1), expected = c(
      0.964365076099296, 0.447213595499958, 0.141421356237310,
      0.989949493661166
    )),
    frac.treated = list(contrast = c(1, 1, 0, 0), expected = c(
      0.964365076099296, 0.538516480713450, 0.616441400296898,
      1.009950493836210
    )),
    majority = list(contrast = c(1, 1, 0, 0), expected = c(
      0.824621125123532, 0.360555127546399, 0.616441400296898,
      0.728010988928052
    ))
  )

  for (exposure in names(cases)) {
    contrast <- cases[[exposure]]$contrast
    observed <- rddinterf:::compute_min_distance(
      exposure, Xi, Xj,
      contrast[1], contrast[2], contrast[3], contrast[4]
    )
    expect_equal(observed, cases[[exposure]]$expected, tolerance = 1e-12)
  }
})

test_that("boundaries cover each observed neighborhood size", {
  neighbors <- list(2:3, c(1, 3, 4), c(1, 2), 2:3)
  treatment <- list(c(1, 1), c(0, 0))

  boundary <- get_boundary(neighbors, treatment, "one.treated")

  expect_s3_class(boundary, "rdd_boundary")
  expect_named(boundary, c("2 -dimensional", "3 -dimensional"))
  expect_match(attr(boundary[[1]], "description"), "Codimension: 2")
})
