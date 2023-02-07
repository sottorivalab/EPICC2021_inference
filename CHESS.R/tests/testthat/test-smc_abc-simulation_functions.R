test_that("helper functions work", {

  expect_equal(calc_ess(rep(0.1, 10)), 10)
  expect_equal(calc_ess(rep(0.1, 5)), 5)
  expect_equal(calc_ess(c(1,1,0,0)), 2)
  
  expect_error(calc_ess(c(1, 1, "A")))
  expect_error(calc_ess(c(NA, 1, 1)))
  expect_error(calc_ess(c(Inf, 1, 1)))
  expect_error(calc_ess(c(-1, 1, 1)))

})

