#
#'
#'Load data from loop_predicton v08 project
#'

load("data-raw/loopPred.v08.ancGR.Rdata")
load("data-raw/loopPred.v08.loopDF.tmp.Rdata")
load("data-raw/loopPred.v08.Rad21.w1001.b10.datamat.Rdata")

devtools::use_data(ancGR)
devtools::use_data(loopDF)
devtools::use_data(datamat)
