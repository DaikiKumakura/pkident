check_against_expected <- function(id, seed = 1L, method = "both") {
  m <- reference_model(id)
  ex <- attr(m, "expected")
  r <- structural_identifiability(m, method = method, seed = seed)
  p <- r$parameters
  list(
    result = r,
    identifiable = sort(p$parameter[p$status == "identifiable"]),
    not_alone = sort(p$parameter[p$status %in% c("non_identifiable", "combination_only")]),
    unresolved = p$parameter[p$status == "unresolved"],
    expected = ex
  )
}

test_that("decisions agree with the expected results for every reference model", {
  for (id in reference_models()$id) {
    x <- check_against_expected(id)
    expect_identical(x$result$decision, "decided", label = id)
    expect_true(all(x$result$agreement$agree), label = paste(id, "methods agree"))
    expect_length(x$unresolved, 0)
    expect_identical(x$identifiable, sort(x$expected$identifiable), label = paste(id, "identifiable"))
    expect_identical(x$not_alone, sort(x$expected$non_identifiable), label = paste(id, "not identifiable alone"))
    expect_identical(nrow(x$result$combinations), length(x$expected$combinations), label = paste(id, "number of combinations"))
    for (cmb in x$expected$combinations) {
      expect_true(combination_identifiable(x$result, cmb), label = paste(id, cmb))
    }
  }
})

test_that("decisions do not depend on the random points", {
  for (id in reference_models()$id) {
    for (seed in 2:4) {
      x <- check_against_expected(id, seed, method = "lie")
      expect_identical(x$identifiable, sort(x$expected$identifiable), label = paste(id, "seed", seed))
      expect_identical(x$not_alone, sort(x$expected$non_identifiable), label = paste(id, "seed", seed))
    }
  }
})

test_that("results are reproducible for a fixed seed and leave the global seed unchanged", {
  set.seed(42); before <- stats::runif(1)
  set.seed(42)
  r1 <- structural_identifiability(reference_model("A9"), method = "lie", seed = 7)
  after <- stats::runif(1)
  expect_identical(before, after)
  r2 <- structural_identifiability(reference_model("A9"), method = "lie", seed = 7)
  expect_identical(r1$parameters, r2$parameters)
  expect_identical(r1$rank, r2$rank)
})

test_that("every result states that the decision is local", {
  r <- structural_identifiability(reference_model("A3"), method = "lie")
  expect_true(r$local)
  expect_output(print(r), "LOCAL")
})

test_that("non-identifiable combinations are not reported as identifiable", {
  r <- structural_identifiability(reference_model("A2"))
  expect_false(combination_identifiable(r, "V"))
  expect_false(combination_identifiable(r, "CL*F"))
  expect_true(combination_identifiable(r, "ka"))
  expect_true(combination_identifiable(r, "V/F"))
})

test_that("an insufficient search is reported as unresolved, not as non-identifiable", {
  r <- structural_identifiability(reference_model("A4"), method = "lie", max_order = 1)
  expect_identical(r$decision, "unresolved")
  expect_false(any(r$parameters$status == "non_identifiable"))
  expect_match(r$reason, "rank still increasing")

  r2 <- structural_identifiability(reference_model("A10"), method = "lie", max_seconds = 0)
  expect_identical(r2$decision, "unresolved")
  expect_false(any(r2$parameters$status == "non_identifiable"))
})

test_that("tables are returned as tibbles", {
  r <- structural_identifiability(reference_model("A7"))
  expect_true(all(c("lie_rank", "sensitivity_rank", "agree") %in% names(as_table(r, "agreement"))))
  expect_s3_class(as_table(r, "parameters"), "tbl_df")
  expect_identical(as_table(r, "combinations")$combination, "Rtot*ke")
  expect_true(all(c("order", "point", "rank") %in% names(as_table(r, "rank"))))
})

test_that("the sensitivity method alone agrees with the expected results", {
  for (id in reference_models()$id) {
    x <- check_against_expected(id, method = "sensitivity")
    expect_identical(x$result$decision, "decided", label = id)
    expect_identical(x$identifiable, sort(x$expected$identifiable), label = paste(id, "identifiable"))
    expect_identical(x$not_alone, sort(x$expected$non_identifiable), label = paste(id, "not identifiable alone"))
    for (cmb in x$expected$combinations) {
      expect_true(combination_identifiable(x$result, cmb), label = paste(id, cmb))
    }
  }
})

test_that("disagreement between the methods gives unresolved, not non-identifiable", {
  # Two time points cannot determine four parameters: the sensitivity rank is
  # too low while the Lie method finds full rank.
  r <- structural_identifiability(reference_model("A4"), times = c(0, 0.001))
  expect_identical(r$decision, "unresolved")
  expect_false(any(r$agreement$agree))
  expect_match(r$reason, "disagree")
  expect_false(any(r$parameters$status == "non_identifiable"))
  expect_output(print(r), "agree at 0 of 5")
})

test_that("invalid input gives a coded error", {
  expect_error(structural_identifiability(list()), class = "pkident_PKI009")
  r <- structural_identifiability(reference_model("A1"), method = "lie")
  expect_error(as_table(r, "agreement"), class = "pkident_PKI009")
  expect_error(structural_identifiability(reference_model("A1"), method = "sensitivity", times = c(0, Inf)),
               class = "pkident_PKI009")
})
