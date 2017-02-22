context("genepairFunctions")


test_that("uniquePairPerGeneBySim runs correclty", {

  gp <- data.frame(
    g1=c(1, 1, 2, 2),
    g2=c(2, 3, 3, 4)
  )

  sim <- c(2, 1, 3, 4)

  uniqPairs <- uniquePairPerGeneBySim(gp, sim)

  # test that all columns are present
  expect_equal(ncol(gp), ncol(uniqPairs))

  # test that only uniq genes are left
  expect_equal(
    length(c(uniqPairs[,1], uniqPairs[,2])),
    length(unique(c(uniqPairs[,1], uniqPairs[,2])))
    )

})


