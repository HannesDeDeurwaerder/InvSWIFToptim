#
#     DATAPREP --> MAKE VARIABLE MATRIX AND SOIL HETEROGENEITY MATRIX
#
################################################################################

#-------------------------------------------------------------------------------
#                     Function all (bio) physical plant variables
#===============================================================================

VarMatrix <- function(itterations=NULL, scenario=NULL, beta=NULL){
  
  # load all the ranges and restriction schemes
  source('./R/dataprepRestrictionRanges.R')
  
  # check which scenario is selected
  if(scenario=="Sc1"){rest=Scenarios[,1]}
  if(scenario=="Sc2"){rest=Scenarios[,2]}
  if(scenario=="Sc3"){rest=Scenarios[,3]}
  if(scenario=="Sc4"){rest=Scenarios[,4]}
  if(scenario=="Sc5"){rest=Scenarios[,5]}

  # create empty data matrix
  Datamatrix=matrix(NA,itterations,ncol=11)
  Names <- c("DBH","Area_alg", "SFtot_alg","ARtot",
             "LA", "Asapwood","Kr","SFdaily",'beta','tstud','hom') 
  colnames(Datamatrix)<-Names 
  
  # Fill in the datamatrix
  
  # 1. DBH
  Datamatrix[,1] <- runif(itterations, min = DBHrange[1], max = DBHrange[2])
  
  # 2. Area
  AreaT=pi*(Datamatrix[,1]/2)^2   # area of the tree in cm
  Datamatrix[,2] <- AreaT
  cm2_to_m2=1/10000   # the model needs area in m²!!
  
  
  # 3. SFtot_alg: algoritmic derived SFtot via coupling with DBH
  Datamatrix[,3] <- 10.7442*exp(0.0406*Datamatrix[,1]) * 10^-3
  # cf. algoritm of Christiano et al, 2015
  
  # 4. ARtot 
  ARtotal = exp(0.88*log(Datamatrix[,2])-2)
  Datamatrix[,4] <- ARtotal * runif(itterations, 
                                    min = ARtotrange[rest[1],1], max = ARtotrange[rest[1],2])
  
  
  #------
  
  # 5. LA: 
  Datamatrix[,5] <- runif(itterations, 
                          min = LArange[rest[2],1], max = LArange[rest[2],2])
  
  # 8. Asap:
  Asap <- 1.582*(Datamatrix[,1])^1.764    # Algoritm of MEIZNER et al 2001
  Datamatrix[,6] <- cm2_to_m2*Asap*runif(itterations, 
                                         min = Asaprange[rest[3],1], max = Asaprange[rest[3],2])
  # make sure that sapwood is never bigger then tree area
  for (i in 1: length(Datamatrix[,6])){
    if(Datamatrix[i,6]>=(Datamatrix[i,2]*cm2_to_m2)){
      Datamatrix[i,6]<-(Datamatrix[i,2]*cm2_to_m2)}
  }
  
  # Kr:
  Datamatrix[,7] <- runif(itterations,
                          min = Krrange[rest[4],1], max = Krrange[rest[4],2])
  # SFDAILY:  
  Datamatrix[,8]<- Datamatrix[,3]*runif(itterations,
                                        min = SFDAILYrange[rest[5],1], max = SFDAILYrange[rest[5],2])
  
  # BETA
  Datamatrix[,9] <- beta          
  
  
  # Tstud
  # NOT RESTRICTED, else do not fill anything
  if(rest[6]==2){
    Datamatrix[,10] <- round(runif(itterations,
                                   min=tstudrange[1]+TCOR, max=tstudrange[2]+TCOR),0)}
  
  #Hom
  Datamatrix[,11]<- runif(itterations,
                          min=Homrange[rest[7],1], max=Homrange[rest[7],2])  
  
  
  return(Datamatrix)
}


#-------------------------------------------------------------------------------
#                     Datasets Soil Heterogeneity
#===============================================================================


SoilHeterogeneity <-function(itterations=NULL, scenario=NULL, DataType=NULL,
                             Z=NULL, SaveData=NULL,Meissner=NULL){
  
  # check which scenario is selected
  source('./R/dataprepRestrictionRanges.R')
  SDfactor=ifelse(scenario=="Sc1",sdMrange[2],sdMrange[1])
  
  
  # allocate correct variables
  avgD<-Meissner$avgDepth
  avgDeut<-Meissner$avgDeuterium;
  sdDeut<-Meissner$sdDeuterium*SDfactor
  avgOxy<-Meissner$avgOxygen;          
  sdOxy<-Meissner$sdOxygen*SDfactor
  avgPSI<-Meissner$avgPsi;                
  sdPSI<-Meissner$sdPsi*SDfactor
  
  
  # to be repeated for all itterations
  Datamatrix=matrix(NA,itterations,length(avgD))
  
  for (ii in 1:length(avgD)){
    if(DataType=="PSI"){
      Datamatrix[,ii]<-rnorm(itterations, avgPSI[ii], (sdPSI[ii]/3))}
    if(DataType=="D2H"){
      Datamatrix[,ii]<-rnorm(itterations, avgDeut[ii], (sdDeut[ii]/3))}
    if(DataType=="D18O"){
      Datamatrix[,ii]<-rnorm(itterations, avgOxy[ii], (sdOxy[ii]/3))}
  }
  
  
  if(DataType=="PSI"){Profiles <- PSImatrix(Datamatrix,Z,scenario,SaveData,avgD)}  
  if(DataType=="D2H"){Profiles <- D2Hsoilmatrix(Datamatrix,Z,scenario,SaveData, avgD)}
  if(DataType=="D18O"){Profiles <- D18Osoilmatrix(Datamatrix,Z, scenario,SaveData, avgD)}
  
  return(Profiles)
}
