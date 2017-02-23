require(GenomicRanges)

context("genepairFunctions")


# some general test instances
testGR <- GRanges(
  c("chr1", "chr1", "chr2"),
  IRanges(
    c(10, 50,  10),
    c(20, 100, 20)
  ),
  strand=c("+", "-", "+")
)

testGP <- data.frame(
  g1=c(1, 2, 3, 1),
  g2=c(1, 1, 1, 2)
)

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
  expect_eqaul(ncol(cisGP), ncol(testGP))
})


