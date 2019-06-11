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

setwd(getwd())
maindir <- (getwd())

source('./R/dataprepStandardPara.R')
source('./R/dataprepRestrictionRanges.R')
source("./R/functionVarMatrix.R")
source("./R/functionIsoSpace.R")
source("./R/SWIFT.R")
source('./R/functionPSIvariance.R')
source('./R/functionRandomDataToIsospace.R')
source("./R/functionLogLik.R")
source("./R/functionOptim.R")
source("./R/create_Rscript_SWIFT.r")

###########################################################
## Parameters
# select a true beta value
Bs=1  # number of itterations over beta trues
Btrues= runif(Bs, min = 0.905, max = 0.995)
FDtotal=c(5,10,25)
scenarios <- 4
    
# define/load the optimization function
itterations=100
RunForWhichIsotope='Both'
    
scenario_FD='Sc4'      # nature equals strongest scenario

# Submission parameters
args <- c('-l walltime=12:00:00','-l nodes=1:ppn=16')
run_per_nodes <- 100

###########################################################

# CREATION OF TRUE FIELD DATA SAMPLES
#-------------------------------------------------------------------------------
compt <- 1
for (iBs in seq(Bs)){
  
  Btrue=Btrues[iBs]
  for (iFD in seq(FDtotal)){
    FDitter=FDtotal[iFD]   # number of samples to generate
    
    # start.time <- Sys.time()
    FD<-RandomDataToIsospace(FDitter, B, scenario_FD,Btrue,"Both",Z, relSF,dZ, TCOR, t, tF, Meissner, n) 
    
    for (iscenar in seq(scenarios)){
        
        folder <- file.path(getwd(),'runs',paste0('run_',sprintf('%05i',compt)))
        dir.create(folder)
        write.csv(FD, file = file.path(folder,paste0( "FD_Btrue_",iBs,"_FD_",FDtotal[iFD],"_scenario_",iscenar,".csv")))
       
        script_file <- file.path(folder,paste0("script",'.R'))
        create_Rscript_SWIFT(script_file,itterations=itterations,RunForWhichIsotope=RunForWhichIsotope,scenario=iscenar,Btrue = Btrue,run_id=compt)
        compt <- compt + 1
    }
  }
}

# Creating Submission files
Nruns <- Bs*scenarios*length(FDtotal)
folder_all <- seq(1,Nruns,run_per_nodes)

for (ifolder in seq(folder_all)){
  folder <- file.path(getwd(),'runs',paste0('run_',sprintf('%05i',folder_all[ifolder])))
  simus <- seq(folder_all[ifolder],min(Nruns,folder_all[ifolder]+run_per_nodes-1))
  
  job_list_file <- file.path(folder,'joblist.txt')
  writeLines(paste0("echo ",'"',"source('script.R')",'"',"| R --vanilla"),con = job_list_file)
  
  cmd_l <- lapply(simus,function(sim){
  line <- file.path(dirname(folder),paste0('run_',sprintf('%05i',sim)))
  write(line,file=job_list_file,append=TRUE)})
  
  launcher_file <- file.path(folder,'launcher.sh')
  writeLines("#!/bin/bash",con = launcher_file)
  write("ml R/3.4.4-intel-2018a-X11-20180131",file=launcher_file,append=TRUE)
  write(paste("mpirun",file.path(maindir,"modellauncher","modellauncher"),job_list_file),file=launcher_file,append=TRUE)
}

# Submitting files
cmd <- "qsub"
for (ifolder in seq(folder_all)){
  folder <- file.path(getwd(),'runs',paste0('run_',sprintf('%05i',folder_all[ifolder])))
  launcher_file <- file.path(folder,'launcher.sh')
  out <- system2(cmd, c(args,launcher_file), stdout = TRUE, stderr = TRUE)
}






