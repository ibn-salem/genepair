context("applyToCisPairs")


test_that("applyToCisPairs runs on small example", {

  rangesGR <- GenomicRanges::GRanges(
    rep(c("chr1", "chr2"), c(3,2)),
    IRanges::IRanges(
      c(100, 200, 300, 100, 200),
      c(150, 250, 350, 150, 250)
    ))

  gp <- data.frame(
   g1=c(1,4,2,1,4),
   g2=c(2,4,3,3,5)
  )

  datamat <- rbind(
    c(10, 20, 30),
    c(15, 26, 40),
    c(100, 2, 0),
    c(10, 20, 30),
    c(15, 26, 40)
  )

  corVal <- applyToCisPairs(gp, rangesGR, datamat)

  expect_equal(length(corVal), 5)
  expect_equal(corVal[2], 1)
  expect_equal(corVal[1], cor(datamat[1,], datamat[2,]))
  expect_equal(corVal, sapply(1:nrow(gp),
                              function(i){
                                cor(datamat[gp[i, 1], ], datamat[gp[i, 2], ])
                              }))
})


