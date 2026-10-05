x <- c(-2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2)

test_that("bounded profiles are identifiable with an interpolated interval", {
  d <- c(40, 12, 3, 0.8, 0, 0.8, 3, 12, 40)
  r <- classify_profile_points(x, d, 3.84, 1)
  expect_identical(r$status, "identifiable")
  expect_identical(r$shape, "bounded")
  expect_equal(r$upper, 0.5 + (sqrt(3.84) - sqrt(3)) / (sqrt(12) - sqrt(3)) * 0.5)
  expect_equal(r$lower, -r$upper)
  # exact for a quadratic profile: dOFV = (x / 0.2)^2 crosses 3.84 at 0.392
  q <- classify_profile_points(x, (x / 0.2)^2, 3.84, 1)
  expect_equal(q$upper, 0.2 * sqrt(3.84))
})

test_that("one-sided profiles report the open direction", {
  d <- c(40, 12, 3, 0.8, 0, 0.2, 0.3, 0.35, 0.36)
  r <- classify_profile_points(x, d, 3.84, 1)
  expect_identical(r$status, "non_identifiable")
  expect_identical(r$shape, "one-sided")
  expect_identical(r$open, "upper")
  expect_true(is.finite(r$lower) && is.na(r$upper))
  r2 <- classify_profile_points(rev(-x), rev(d), 3.84, 1)
  expect_identical(r2$open, "lower")
})

test_that("flat and shallow profiles are non-identifiable", {
  r <- classify_profile_points(x, c(0.4, 0, 0.2, -0.3, 0, 0.1, 0.5, 0.3, 0.2), 3.84, 1)
  expect_identical(c(r$status, r$shape, r$open), c("non_identifiable", "flat", "both"))
  r <- classify_profile_points(x, c(2, 1, 0.5, 0.1, 0, 0.1, 0.5, 1, 2), 3.84, 1)
  expect_identical(c(r$status, r$shape), c("non_identifiable", "shallow"))
})

test_that("failed fits give unresolved", {
  r <- classify_profile_points(x, c(40, NA, 3, 0.8, 0, 0.8, 3, 12, 40), 3.84, 1)
  expect_identical(c(r$status, r$shape), c("unresolved", "failed fit"))
})

test_that("a lower objective found while profiling becomes the reference", {
  # the original fit (x = 0) is 2 units above the minimum at x = 0.5
  d <- c(40, 20, 6, 1, 0, -1, -2, 3, 30)
  r <- classify_profile_points(x, d, 3.84, 1)
  expect_equal(r$drop, 2)
  expect_equal(r$centre, 0.5)
  expect_identical(r$status, "identifiable")
  expect_true(r$lower < 0 && r$upper > 0.5 && r$upper < 1)
})

test_that("profiles are interpreted against the structural result", {
  expect_identical(interpret_profile("identifiable", "identifiable"), "determined by these data")
  expect_identical(interpret_profile("non_identifiable", "identifiable"),
                   "structurally identifiable but not determined by these data")
  expect_identical(interpret_profile("non_identifiable", "combination_only"), "structurally non-identifiable")
  expect_match(interpret_profile("identifiable", "non_identifiable"), "^conflict")
  expect_identical(interpret_profile("identifiable", "unresolved"), "unresolved")
})

test_that("invalid arguments give coded errors", {
  skip_if_not_installed("nlmixr2extra")
  expect_error(classify_profiles(list()), class = "pkident_PKI009")
})

test_that("an nlmixr2 fit is profiled and classified", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2est")
  skip_if_not_installed("nlmixr2extra")
  fit <- fit_profile_case("oral", nid = 12)
  expect_error(classify_profiles(fit, which = "nope"), class = "pkident_PKI009")
  r <- classify_profiles(fit, which = "tcl", steps = c(0.25, 1),
                         structural = structural_identifiability(reference_model("A3"), method = "lie"),
                         map = c(tcl = "CL"))
  expect_identical(r$parameters$status, "identifiable")
  expect_identical(r$parameters$interpretation, "determined by these data")
  expect_true(r$parameters$lower < log(2) && log(2) < r$parameters$upper)
  expect_s3_class(as_table(r, "profiles"), "tbl_df")
  expect_output(print(r), "explored range")
})
