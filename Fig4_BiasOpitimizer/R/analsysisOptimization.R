#
#     OPTIMIZATION USING THE GENSA FUNCTION
#
#############################################################
rm(list=ls())

library('lhs')
library("SWIFT")
library("GenSA")
library('sm')
library('ks')

# load standard model parameters
source('./R/dataprepStandardPara.R')
source('./R/dataprepRestrictionRanges.R')

source("./R/functionVarMatrix.R")
source("./R/functionIsoSpace.R")
source("./R/SWIFT.R")
source('./R/functionPSIvariance.R')
source('./R/functionRandomDataToIsospace.R')
source("./R/functionLogLik.R")



# select a true beta value
Bs=10  # number of itterations over beta treus
Btrues= runif(Bs, min = 0.905, max = 0.995)
FDtotal=c(5,25,50)

# run for different sizes of collected data
for(FF in 1:length(FDtotal)){
  
  
  B_est_Sc1=B_est_solo_Sc1=
    B_est_Sc2=B_est_solo_Sc2=
    B_est_Sc3=B_est_solo_Sc3=
    B_est_Sc4=B_est_solo_Sc4= rep(NA, Bs)
  
  for (hh in 1:Bs){
    
    
    # CREATION OF TRUE FIELD DATA SAMPLES
    #-------------------------------------------------------------------------------
    FDitter=FDtotal[FF]   # number of samples to generate
    scenario_FD='Sc4'      # nature equals strongest scenario
    Btrue=Btrues[hh]
    
    # start.time <- Sys.time()
    FD<-RandomDataToIsospace(FDitter, B, scenario_FD,Btrue,"Both",Z, relSF,dZ, TCOR, t, 
                             tF, Meissner, n) 
    # FDsolo<-RandomDataToIsospace(1, scenario_FD,Btrue,"D2H",Z, relSF,dZ, TCOR, t, 
    #                              tF, Meissner, n) 
    # end.time <- Sys.time()
    # time.taken <- end.time - start.time
    # time.taken
    
    
    #----------------------------------------------  
    #----------------------------------------------  
    
    # OPTIMIZATION
    #---------------
    
    # define/load the optimization function
    itterations=250
    RunForWhichIsotope='Both'
    
    # OPTIMIZATION WITH GenSA
    # start.time <- Sys.time()
    # test<-GenSA(par=0.95, fn=LogLikOptim, lower = 0.905, upper = 0.995,
    #             control=list(maxit=3, simple.function=TRUE))
    # end.time <- Sys.time()
    # time.taken <- end.time - start.time
    
    # OPTIMIZATION WITH OPTIM "BRENT"
    
    Scenario="Sc1"
      source("./R/functionOptim.R")
      B_Sc1<-optim(par=0.95, fn=LogLikOptim, lower = 0.905, upper = 0.995, 
                 method = "Brent")
      print("SC1_Done")
    rm(Scenario)
    
    
    Scenario="Sc2"
      source("./R/functionOptim.R")
      B_Sc2<-optim(par=0.95, fn=LogLikOptim, lower = 0.905, upper = 0.995, 
                 method = "Brent")
      print("SC2_Done")
    rm(Scenario)
    
    
    Scenario="Sc3"
      source("./R/functionOptim.R")
      B_Sc3<-optim(par=0.95, fn=LogLikOptim, lower = 0.905, upper = 0.995, 
                 method = "Brent")
      print("SC3_Done")
    rm(Scenario)
    
    
    Scenario="Sc4"
      source("./R/functionOptim.R")
      B_Sc4<-optim(par=0.95, fn=LogLikOptim, lower = 0.905, upper = 0.995, 
                 method = "Brent")
      print("SC4_Done")
    rm(Scenario)
    
    
    B_est_Sc1[hh]     <-B_Sc1$par
    B_est_Sc2[hh]     <-B_Sc2$par
    B_est_Sc3[hh]     <-B_Sc3$par
    B_est_Sc4[hh]     <-B_Sc4$par

  
  DF<-data.frame(
    Btrue=Btrues,
    Bsc1=B_est_Sc1,
    Bsc2=B_est_Sc2,
    Bsc3=B_est_Sc3,
    Bsc4=B_est_Sc4 #,
  )
  
  write.csv(DF, file = paste0( "./csvFiles/bias_",FDitter,"_",hh,".csv"))

  }
  
print(paste0('Yippie, sample size of ',FDitter,' completed'))
  
  
  
}


source('./R/plotBoxPlot.R')

