
context("genepairFunctions")


# some general test instances
testGR <- GenomicRanges::GRanges(
  c("chr1", "chr1", "chr2"),
  IRanges::IRanges(
    c(10, 50,  10),
    c(20, 100, 20)
  ),
  strand=c("+", "-", "+")
)

testGP <- data.frame(
  g1=c(1, 2, 3, 1),
  g2=c(1, 1, 1, 2)
)

test_that("getPairIDsorted() retunrs correct IDs on testGP dataset.", {

  ids <- getPairIDsorted(testGP)
  expect_equal(ids, c("1_1", "1_2", "1_3", "1_2"))

})

test_that("getPairIDsorted() retunrs same id for permuted order on random
          example with indexes and letters.", {

  randGP <- data.frame(
    g1=sample.int(100, 10^4, replace=TRUE),
    g2=sample.int(100, 10^4, replace=TRUE)
  )

  lettersGP <- data.frame(
    g1=sample(letters, 10^3, replace=TRUE),
    g2=sample(letters, 10^3, replace=TRUE),
    stringsAsFactors = FALSE
  )

  expect_equal(getPairIDsorted(randGP), getPairIDsorted(randGP[,2:1]))
  expect_equal(getPairIDsorted(lettersGP), getPairIDsorted(lettersGP[,2:1]))

})


test_that("containsGenePairs works on custom example", {

  A <- data.frame(
    g1=c(2, 1, 3),
    g2=c(1, 2, 3)
  )

  B <- data.frame(
    g1=c(1, 3),
    g2=c(2, 4)
  )

  expect_equal(containsGenePairs(A, B), c(TRUE, TRUE, FALSE))
  expect_equal(containsGenePairs(B, A), c(TRUE, FALSE))

})

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


test_that("filterForCisPairs works on example", {

  cisGP <- filterForCisPairs(testGP, testGR)

  expect_equal(nrow(cisGP), 3)
  expect_equal(ncol(cisGP), ncol(testGP))
})


