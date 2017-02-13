#
#'
#'Load data from loop_predicton v08 project
#'

load("data-raw/loopPred.v08.ancGR.Rdata")
load("data-raw/loopPred.v08.loopDF.tmp.Rdata")
load("data-raw/loopPred.v08.Rad21.w1001.b10.datamat.Rdata")

# remove duplicates
loopDF <- loopDF[!duplicated(loopDF[,1:2]),]

# sort rows according to first two columns
loopDF <- loopDF[order(loopDF[,1], loopDF[,2]),]

devtools::use_data(ancGR)
devtools::use_data(loopDF)
devtools::use_data(datamat)
