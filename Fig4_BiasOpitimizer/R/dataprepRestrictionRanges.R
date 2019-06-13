#
#   Dataprep Restriction scheme - Dataranges
#
#######################################################
# variable notation is Xr and Xnr
# where r means restricted range
# where nr means non-restricted range

rownames=c('restricted','non-restricted')

# DBH (in cm)
DBHrange = c(10,60)

# ARi
ARtotr=c(0.9,1.1)
ARtotnr=c(0.5,1.5)
ARtotrange=rbind(ARtotr, ARtotnr)
row.names(ARtotrange)<-rownames

# LA
LAr=c(0.10,0.16)
LAnr=c(0.0411,0.451)
LArange=rbind(LAr, LAnr)
row.names(LArange)<-rownames

# Asap
Asapr=c(0.9,1.1)   
Asapnr=c(0.5,1.5)
Asaprange=rbind(Asapr, Asapnr)
row.names(Asaprange)<-rownames

# Kr
Krr=c(9*10^-10,11*10^-10)
Krnr=c(1*10^-10,14*10^-10)
Krrange=rbind(Krr, Krnr)
row.names(Krrange)<-rownames

# SFDAILY
SFDAILYr=c(0.9,1.1)   
SFDAILYnr=c(0.5,1.5)
SFDAILYrange=rbind(SFDAILYr, SFDAILYnr)
row.names(SFDAILYrange)<-rownames
# tstud
tstudrange=c(8*tF, tF*14)  # only 1, otherwhise SWIFT approach!!
# Hom
Homr=c(1.3,1.3)
Homnr=c(0,25)
Homrange=rbind(Homr, Homnr)
row.names(Homrange)<-rownames

# Soil heterogeneity
sdMr=1
sdMnr=3
sdMrange=rbind(sdMr, sdMnr)
row.names(sdMrange)<-rownames




# Define the scenarios  with restriction
#----------------------------------------------
rownames=c('ARtot','LA','Asap','Kr','SFDAILY','tstud','Hom','Meissner SD')
# no restriction at all: --> 2 is non restricted, 1 is restricted
sc1=c(2,2,2,2,2,2,2,2)
sc2=c(2,2,2,2,2,2,2,1)
sc3=c(2,2,2,2,2,1,1,1)
sc4=c(1,1,1,1,1,1,1,1)
sc5=c(1,1,1,1,1,2,2,1)
Scenarios=cbind(sc1,sc2,sc3,sc4, sc5)
colnames(Scenarios)<-c('Sc1','Sc2','Sc3','Sc4','Sc5')
row.names(Scenarios)<-rownames
