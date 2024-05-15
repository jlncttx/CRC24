#---------------------------------------------------------------------------------
# Source libraries / functions
#---------------------------------------------------------------------------------
source("scan_source.R")
  
#---------------------------------------------------------------------------------
# Namelist of the event
#---------------------------------------------------------------------------------
  
# Variable of interest
VAR="GSAT"
  
# Years used for data samples
YEARS=1940:2024

# Scenario used for initial estimate of the forced response
#SCEN="45"
##################################################################################
for (SCEN in c("26","45","85")){ ### LOOP ON SCENARIO
##################################################################################

# Years to scan (possibly a subset of YEARS)
SCAN.YRS=YEARS

# Time windows for searching the event (durations in days)
WINDOWS=c(1,3,5,7,10,15,20,25,30,40,50,60,75,90,120,150,180,210,240,270,300,365)

# Location
LOC="Global"

# Year of reference for counter-factual world
YREF=1940

# Include first/last years for annual-max/min samples (are Txnday or Tnnday representative?)
FIRST.YEAR.MAX=FIRST.YEAR.MIN=T
LAST.YEAR.MAX=(max(YEARS)==2023)
LAST.YEAR.MIN=(max(YEARS)==2023)

# First calendar day for annual blocks (e.g. 1 for Jan-Dec, 182 for Jul-Jun)
BLOCK.DAY1.MAX=1
BLOCK.DAY1.MIN=182

# Smoothing parameters (nb of dof) for non-stationary normals
nDF.f=18  # annual cycle (365d)
nDF.g=10  # multi-model mean of yearly tas (250y) -- not used here
nDF.h=9   # change in annual cycle (365d)

# Key to force a zero trend
FORCE.ZERO.TREND=F

# Output directory
ODIR=paste0(SCANDIR,"scans_1d/tglo_",SCEN,"/")
system(paste("mkdir -p",ODIR))

#---------------------------------------------------------------------------------
# Get daily obs data. We need:
#  - XLOC   = data.frame $date $var with daily values at location LOC
#  - PERIOD = vector of years corresponding to XLOC
#---------------------------------------------------------------------------------
# Data for 1D-analysis
XLOC=myno("../data/era5/tglo_era5_day_1940-2024.nc")
XLOC=data.frame(date=mdy2yyyymmdd(time2mdy(XLOC$time,XLOC$timeu)),var=XLOC$var-273.15)

# Restrict to YEARS
XLOC=XLOC[which(yyyymmdd2mdy(XLOC$date)$y %in% YEARS),]

# Add a year of NAs before/after (edge effects) + delete feb29s
PERIOD=(min(YEARS)-1):(max(YEARS)+1)
dum=data.frame(date=mdy2yyyymmdd(makedates(y=PERIOD,avg="day",caltype="noleap")),var=NA)
dum$var[which(dum$date %in% XLOC$date)]=XLOC$var[which(XLOC$date %in% dum$date)]
XLOC=dum

#---------------------------------------------------------------------------------
# Get yearly forced response. We need:
#  - XMMM        = vector of initial estimate (e.g. raw multi-model mean)
#  - CMIP.PERIOD = vector of years corresponding to XMMM
#  - TREND       = vector of adjusted estimate (e.g. smoothed XMMM) restricted to PERIOD 
#---------------------------------------------------------------------------------
# Here we use forced response from Aurelien's constraints
CONST=read.table(paste0("../scripts_Aurel/ts_p",SCEN,"_cons.txt"),header=T)
CMIP.PERIOD=CONST[,1]
XMMM=CONST[,3]
TREND=CONST[which(CONST[,1] %in% PERIOD),3]

#---------------------------------------------------------------------------------
# Scan
#---------------------------------------------------------------------------------
# Load functions
source("scan_1d_sub.R")

if (2==1){ ###
# Loop on methods
#for (m in c("calend","annmax","annmin","locmax","locmin",
#            "annmax1","annmax2","annmax3",
#            "locmax15","locmax30")){
for (m in c("calend")){
  # Method to use
  METHOD=m
  # Source scan subroutine
  source("scan_1d_exe.R")
  # Save result
  save(list="SCAN",file=paste0(ODIR,"scan_",LOC,"_",VAR,"_",METHOD,".Rdata"))
}

} ###

#---------------------------------------------------------------------------------
# Redo the scan but on anomalies wrt. non-stationary normals
#---------------------------------------------------------------------------------
#if (2==1){ ###

# Non-stationary normals
NORM=fit_onto_forced_response(XLOC$var,TREND,nDF.f,nDF.h)

# Save things that are modified below
save=list(VAR=VAR,XLOC=XLOC,nDF.f=nDF.f,nDF.h=nDF.h)

# Modify things for the new scan
VAR=paste0(VAR,"-ano")
XLOC$var=XLOC$var-NORM
nDF.f=1
nDF.h=1

# Specific options for scanning anomalies
#LAST.YEAR.MAX=(max(YEARS)==2024)
BLOCK.DAY1.MAX=182
FORCE.ZERO.TREND=T

# New scan
#for (m in c("calend","annmax","annmin","locmax","locmin",
#            "annmax1","annmax2","annmax3",
#            "locmax15","locmax30")){
#for (m in c("calend","annmax","annmax3","locmax","locmax15","locmax30")){
#for (m in c("calend","annmax3","locmax15")){
for (m in c("annmax2")){
  METHOD=m
  source("scan_1d_exe.R")
  save(list="SCAN",file=paste0(ODIR,"scan_",LOC,"_",VAR,"_",METHOD,".Rdata"))
}

# Reinitialize things
VAR=save$VAR
XLOC=save$XLOC
nDF.f=save$nDF.f
nDF.h=save$nDF.h

#} ###

#---------------------------------------------------------------------------------
# Redo the scan but on anomalies wrt. a reference (fixed) climatology
#---------------------------------------------------------------------------------
if (2==1){ ###

# Stationary normals
CLIM=apply(matrix(XLOC$var,nrow=365)[,which(PERIOD %in% 1991:2020)],1,mean,na.rm=T)
CLIM=mysmoothy(CLIM,df=nDF.f,per=T)

# Save things that are modified below
save=list(VAR=VAR,XLOC=XLOC,nDF.f=nDF.f,nDF.h=nDF.h)

# Modify things for the new scan
VAR=paste0(VAR,"-anoref")
XLOC$var=XLOC$var-rep(CLIM,length(PERIOD))

# Specific options
FORCE.ZERO.TREND=F

# New scan
#for (m in c("calend","annmax","annmin","locmax","locmin",
#            "annmax1","annmax2","annmax3",
#            "locmax15","locmax30")){
for (m in c("calend","annmax3","locmax15")){
#for (m in "calend"){
  METHOD=m
  source("scan_1d_exe.R")
  save(list="SCAN",file=paste0(ODIR,"scan_",LOC,"_",VAR,"_",METHOD,".Rdata"))
}

# Reinitialize things
VAR=save$VAR
XLOC=save$XLOC

} ###


##################################################################################
} ### END LOOP ON SCENARIO
##################################################################################

