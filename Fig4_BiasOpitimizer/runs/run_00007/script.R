rm(list=ls())

library('lhs')
library('SWIFT')
library('sm')
library('ks')
library('GenSA')

setwd('../../')
source('./R/dataprepStandardPara.R')
source('./R/dataprepRestrictionRanges.R')
source('./R/functionVarMatrix.R')
source('./R/functionIsoSpace.R')
source('./R/SWIFT.R')
source('./R/functionPSIvariance.R')
source('./R/functionRandomDataToIsospace.R')
source('./R/functionLogLik.R')
source('./R/functionOptim.R')

itterations = 100
Btrue = 0.95024468133226
RunForWhichIsotope =  'Both'
scenario =  3
run_id =  7

folder <-paste0('./runs/run_',sprintf('%05i',run_id))
setwd(folder)
FD_file <- list.files(getwd(),pattern = 'FD_*')
FD <- read.csv(file.path(getwd(),FD_file))

setwd('../../')
Scenario <- switch(scenario,'Sc1','Sc2','Sc3','Sc4')
B_Sc1 <- optim(par=0.95, fn=LogLikOptim, lower = 0.905, upper = 0.995, method = 'Brent')

DF <- data.frame(Btrue=Btrue,Bsc=B_Sc1$par)
write.csv(DF, file = file.path(getwd(),paste0('./runs/run_',sprintf('%05i',run_id)),'results.csv'))
