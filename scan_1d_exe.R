#-------------------------------------------------------------------------------------
# Computes scan for LOC / SCAN.YRS / METHOD
# Returns list 'SCAN' with results of compute.stat.1d() for all years + metadata
#-------------------------------------------------------------------------------------
VERBOSE=paste(VAR,LOC,paste(range(SCAN.YRS),collapse=" "),paste(METHOD,collapse=" "))

#-------------------------------------------------------------------------------------
# Arguments of compute.stat.1d() corresponding to METHOD
#-------------------------------------------------------------------------------------
# Retrieve default parameters from generic table
args=read.table("scan_methods.txt",header=T)
args=args[which(args$method==METHOD),]

# Additional parameters specified in the main script
fiyr=layr=T ; bld1=NA ; fztr=FORCE.ZERO.TREND
if (substr(METHOD,1,6)=="annmax"){
  fiyr=FIRST.YEAR.MAX
  layr=LAST.YEAR.MAX
  bld1=BLOCK.DAY1.MAX
}
if (substr(METHOD,1,6)=="annmin"){
  fiyr=FIRST.YEAR.MIN
  layr=LAST.YEAR.MIN
  bld1=BLOCK.DAY1.MIN
}

# Save parameters for output
params=list(long.name=args$details,annual.maxima=args$annm,calendar.flex=args$calf,
            distrib=args$dist,fixed.ksi=args$fksi,smooth.fit.sdf=args$smfs,exceedance=args$excd,
            first.year=fiyr,last.year=layr,block.day1=bld1,
            force.zero.trend=fztr,ndf.f=nDF.f,ndf.g=nDF.g,ndf.h=nDF.h)

#-------------------------------------------------------------------------------------
# Data time series ts (length 365 x nYRS)
#-------------------------------------------------------------------------------------
print(paste(VERBOSE,": preparing data"))
print(".. time series")
ts=XLOC$var

#-------------------------------------------------------------------------------------
# Forced response fr (length nYRS)
#-------------------------------------------------------------------------------------
print(".. forced response")
fr=TREND

#-------------------------------------------------------------------------------------
# Trend tr (size depends on METHOD)
#-------------------------------------------------------------------------------------
print(".. trend")

# Trend time series tr with seasonal cycle for calendar method (length 365 x nYRS)
if (substr(METHOD,1,6) %in% c("calend","locmax","locmin"))
  tr=fit_onto_forced_response(ts,fr,nDF.f,nDF.h)

# Trend time series tr for annmax/annmin methods, with a dependance to the time
# window, i.e. TX1d vs. TX10d have different trends (matrix nWIN x nYRS)
if (substr(METHOD,1,6) %in% c("annmax","annmin")){
  # - compute samples
  tsx=compute.array.of.data.samples(ts,ann=args$annm,exc=args$excd,fir=fiyr,las=layr,blo=bld1)[1,,]
  #  - fit values onto forced response
  cfr=fr-mean(fr)
  f_tsx=t(apply(tsx,1,function(x) lm(x~cfr)$coef))
  #  - smooth scaling factors over time windows
  sf_tsx=mysmoothy(f_tsx[,2],sdf=2)
  #  - computes tr matrices
  tr=matrix(rep(f_tsx[,1],nYRS)+rep(sf_tsx,nYRS)*rep(cfr,each=nWIN),nWIN,nYRS)
}

#-------------------------------------------------------------------------------------
# Initialize scan
#-------------------------------------------------------------------------------------
print(paste(VERBOSE,": initializing scan"))
YEVE=PERIOD[1] ;  iYEVE=which(PERIOD==YEVE) ; iDEVE=which(DATES$y==YEVE)

#-------------------------------------------------------------------------------------
# Data values (dum) samples (datasp) and trends (datatr) 
#-------------------------------------------------------------------------------------
print(".. computes compute.stat.1d() matrix inputs")

# Values for all years (faster than using compute.matrix.of.event.values in the loop below)
dum=compute.array.of.data.samples(ts)

# Sample (depends on METHOD through parameters)
{
if (!args$annm & args$calf==0) datasp=dum
else datasp=compute.array.of.data.samples(ts,args$annm,args$calf,args$excd,fiyr,layr,bld1)
}

# Trend (depends on METHOD through tr computation above)
datatr=compute.array.of.trend.values(ts,tr,force.zero.trend=fztr)

#-------------------------------------------------------------------------------------
# Scan
#-------------------------------------------------------------------------------------
print(paste(VERBOSE,": starting scan"))

# Initialize SCAN
SCAN=list(var=VAR,location=LOC,years=YEARS,period=PERIOD,scan.years=SCAN.YRS,yref=YREF,windows=WINDOWS,
          method=METHOD,params=params,data=datasp,trend=datatr,cmip.period=CMIP.PERIOD,cmip.mmm=XMMM,
          forced.response=TREND,scan=list())

# Loop on requested years
for (YEVE in SCAN.YRS){
  print(paste(LOC,": scan",METHOD,YEVE))
  # Indices of year/days of the considered event (hard coded..)
  iYEVE=which(PERIOD==YEVE) ; iDEVE=which(DATES$y==YEVE)
  # Matrix of event values
  #dateve=compute.matrix.of.event.values(ts)
  dateve=dum[,,iYEVE]
  # Compute stat 1d
  out=compute.stat.1d(input=list(dateve=dateve,datasp=datasp,datatr=datatr),
		      ann=args$annm,cal=args$calf,dis=args$dist,
		      fix=args$fksi,smo=args$smfs,exc=args$excd)
  # Save in SCAN
  SCAN$scan[[paste0("Y",YEVE)]]=list(p1=out$p1,p0=out$p0,value=out$value,fit1=out$fit1,fit0=out$fit0)
}

