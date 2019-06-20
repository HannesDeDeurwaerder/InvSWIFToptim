
#-------------------------------------------------------------------------------
#     1. Defining non-plant specific Parameters
#-------------------------------------------------------------------------------
n= 20             # Number of days studied
tF= 60            # Time frequence of measurements per hour [in measurments 
# per h] 
t= seq(0,24*n,length.out = 24*tF*n)     # Discrete time vector [in h]
dZ=0.001          # Thickness of sampled layer [in m]
L=1               # maximum soil depth [in m]
Z=seq(dZ,L,dZ)    # Discrete depth vector centered [in m]
CTpsi=101.97      # Conversion factor between MPa and m H2O
TCOR=24*60*(n-2)     # time correction

# relative sapflow day --> this is standard in the SWIFT MODEL CODE
data(SFday) # standard from SWIFT model
uch <- (60)*60*10^3 
SFdayA=SFday/uch
SF  <- rep(SFdayA,n)		# repetition of SF day over n prefered days
relSF <- (SF/ sum(SFdayA))/tF   # multiply by n to have the entire period but 
#devide by tF to get the fequency correct

# open document of data Meissner
Meissner<-read.csv("./data/DeuteriumPsi.csv")

# B is kept standard over all runs!!
Bbeta=0.95
BR0=438688
B<-Bprep(Bbeta, BR0, Z)
