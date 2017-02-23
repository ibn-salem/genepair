context("applyToClosePairs")


test_that("applyToClosePairs runs on small example", {

  rangesGR <- GenomicRanges::GRanges(
    rep(c("chr1", "chr2"), c(3,2)),
    IRanges::IRanges(
      c(100, 200, 300, 100, 200),
      c(150, 250, 350, 150, 250)
    ))

  GenomeInfoDb::seqlengths(rangesGR) <- c(1000, 1000)

  gp <- data.frame(
    g1=c(1,4,2,1,4),
    g2=c(2,4,3,3,5)
  )

  # order gp
  gp <- gp[order(gp[,1], gp[,2]),]

  datamat <- rbind(
    c(10, 20, 30),
    c(15, 26, 40),
    c(100, 2, 0),
    c(10, 20, 30),
    c(15, 26, 40)
  )

  corVal <- applyToClosePairs(gp, rangesGR, datamat, maxDist=1000)

  expect_equal(length(corVal), 5)
  expect_equal(corVal[4], 1)
  expect_equal(corVal[1], cor(datamat[1,], datamat[2,]))
  realCor <- cor(t(datamat))[as.matrix(gp[,1:2])]
  expect_equal(corVal, realCor)
  expect_equal(
    corVal,
    cor(t(datamat))[as.matrix(gp[,1:2])]
  )

})

test_that("applyToClosePairs runs on large exmaple dataset", {

  corVal <- applyToClosePairs(loopDF, ancGR, datamat, fun=cor, maxDist=10^6)

  expect_equal(length(corVal), nrow(loopDF))
  expect_equal(
      corVal[20],
      cor(datamat[loopDF[20,1],], datamat[loopDF[20,2],])
      )
  expect_equal(
      corVal[100000:100010],
      sapply(100000:100010, function(i) cor(datamat[loopDF[i,1],], datamat[loopDF[i,2],]))
      )
})


