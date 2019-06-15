create_Rscript_SWIFT <- function(file,itterations=10,RunForWhichIsotope='Both',scenario=1,Btrue,run_id,N_FD,iBs) {
  
  writeLines("rm(list=ls())",con = file)
  
  write("",file=file,append=TRUE)
  write("library('lhs')",file=file,append=TRUE)
  write("library('SWIFT')",file=file,append=TRUE)
  write("library('sm')",file=file,append=TRUE)
  write("library('ks')",file=file,append=TRUE)
  write("library('GenSA')",file=file,append=TRUE)
  
    # load standard model parameters
  write("",file=file,append=TRUE)
  write(paste("setwd('../../')"),file=file,append=TRUE)
  write("source('./R/dataprepStandardPara.R')",file=file,append=TRUE)
  write("source('./R/dataprepRestrictionRanges.R')",file=file,append=TRUE)
  write("source('./R/functionVarMatrix.R')",file=file,append=TRUE)
  write("source('./R/functionIsoSpace.R')",file=file,append=TRUE)
  write("source('./R/SWIFT.R')",file=file,append=TRUE)
  write("source('./R/functionPSIvariance.R')",file=file,append=TRUE)
  write("source('./R/functionRandomDataToIsospace.R')",file=file,append=TRUE)
  write("source('./R/functionLogLik.R')",file=file,append=TRUE)
  write("source('./R/functionOptim.R')",file=file,append=TRUE)
  
  write("",file=file,append=TRUE)
  write(paste("itterations =",itterations),file=file,append=TRUE)
  write(paste("Btrue =",Btrue),file=file,append=TRUE)
  write(paste("RunForWhichIsotope = ",paste0("'",RunForWhichIsotope,"'")),file=file,append=TRUE)
  write(paste("scenario = ",scenario),file=file,append=TRUE)
  write(paste("run_id = ",run_id),file=file,append=TRUE)
  write(paste("N_FD = ",N_FD),file=file,append=TRUE)
  write(paste("iBs = ",iBs),file=file,append=TRUE)  

  write("",file=file,append=TRUE)
  write(paste("folder <-paste0('./runs/run_',sprintf('%05i',run_id))"),file=file,append=TRUE)
  write(paste("setwd(folder)"),file=file,append=TRUE)
  #write("FD_file <- list.files(getwd(),pattern = 'FD_*')",file=file,append=TRUE)  
  #write("FD <- read.csv(file.path(getwd(),FD_file))",file=file,append=TRUE)
  write("FD_file <- file.path(paste0('FD_Btrue_',iBs,'_FD_',N_FD,'_scenario_',scenario,'.csv'))",file=file,append=TRUE)  
  write("FD <- read.csv(FD_file)",file=file,append=TRUE)
  
  write("",file=file,append=TRUE)
  write(paste("setwd('../../')"),file=file,append=TRUE)
  write("Scenario <- switch(scenario,'Sc1','Sc2','Sc3','Sc4')",file=file,append=TRUE)
  write("B_Sc1 <- optim(par=0.95, fn=LogLikOptim, lower = 0.905, upper = 0.995, method = 'Brent')",file=file,append=TRUE)
  
  write("",file=file,append=TRUE)
  write("DF <- data.frame(Btrue=Btrue,Bsc=B_Sc1$par)",file=file,append=TRUE)
  write("write.csv(DF, file = file.path(getwd(),paste0('./runs/run_',sprintf('%05i',run_id)),'results.csv'))",file=file,append=TRUE)
  #print("SC1_Done")

}
