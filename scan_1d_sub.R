#---------------------------------------------------------------------------------
# A few helpful assigns used in scan_1d_sub.R and scan_1d_exe.R
#---------------------------------------------------------------------------------
nWIN=length(WINDOWS)                        # number of time windows
nYRS=length(PERIOD)                         # number of years
nDAY=365                                    # number of calendar days
DATES=yyyymmdd2mdy(XLOC$date)               # dates as data.frame $m $d $y
iYREF=which(PERIOD==YREF)                   # index of the year of reference

#---------------------------------------------------------------------------------
# Subroutine of compute.stat.1d() that computes the nDAY x nWIN matrix
# of event values from ts (vector)
#---------------------------------------------------------------------------------
compute.matrix.of.event.values=function(ts){
  # Initialize output
  out=matrix(NA,nDAY,nWIN)
  # Loop on time windows
  for (iw in 1:nWIN){ w=WINDOWS[iw]
    tsw=myrunavg(ts,w)
    # Loop on calendar days
    for (id in 1:nDAY){ d=iDEVE[id]
      out[id,iw]=tsw[d]
    }
  }
  out
}

#---------------------------------------------------------------------------------
# Subroutine of compute.stat.1d() that computes the nDAY x nWIN x nYRS array
# of data samples from ts (vector)
#---------------------------------------------------------------------------------
compute.array.of.data.samples=function(ts,annual.maxima=FALSE,calendar.flex=0,
                                       exceedance="above",first.year=TRUE,last.year=TRUE,
                                       block.day1=1){
  # Initialize output
  out=array(NA,dim=c(nDAY,nWIN,nYRS))
  # Loop on time windows
  for (iw in 1:nWIN){ w=WINDOWS[iw]
    # Running mean of ts over w
    tsw=myrunavg(ts,w)
    # If annual.maxima, same sample for all calendar days (only id=1 is filled)
    if (annual.maxima){
      if (block.day1>1) tsw=c(tsw[1:nDAY],tsw[(nDAY+block.day1):length(tsw)],rep(NA,block.day1))
      for (y in PERIOD) out[1,iw,which(PERIOD==y)]=mymax(tsw[which(DATES$y==y)],exceedance)
    }
    # Else, loop on calendar days
    if (!annual.maxima){ for (id in 1:nDAY){ d=iDEVE[id]
      idays=which(DATES$m==DATES$m[d] & DATES$d==DATES$d[d])
      if (calendar.flex==0) out[id,iw,]=tsw[idays]
      else {
        flex=-calendar.flex:calendar.flex
        idays=rep(idays,each=length(flex))+rep(flex,nYRS)
        idays[idays<1]=1 ; idays[idays>nrow(DATES)]=nrow(DATES)
        out[id,iw,]=apply(matrix(tsw[idays],nrow=length(flex)),2,mymax,exceedance)
      }
    }}
    # Remove first and/or last years if requested (e.g. for annual.maxima when the year is incomplete)
    if (!first.year) out[,,2]=NA
    if (!last.year)  out[,,nYRS-1]=NA
  }
  out
}

#---------------------------------------------------------------------------------
# Subroutine of compute.stat.1d() that computes the nDAY x nWIN x nYRS array
# of trend values from tr (vector of length ts or nYRS or matrix nWIN x nYRS)
#---------------------------------------------------------------------------------
compute.array.of.trend.values=function(ts,tr,annual.maxima=FALSE,force.zero.trend=F){
  # If zero is forced, easy
  if (force.zero.trend) out=array(0,dim=c(nDAY,nWIN,nYRS))
  else {
    # Initialize output
    out=array(NA,dim=c(nDAY,nWIN,nYRS))
    # If annual.maxima and daily trend, compute annual means of tr
    if (annual.maxima & length(ts)==length(tr)) tr=apply(matrix(tr,nrow=nDAY),2,mean,na.rm=T)
    # If annual trend (no seasonal cycle), repeat annual value
    if (length(tr)==nYRS){ for (i in 1:nYRS) out[,,i]=tr[i] }
    # If annual trend with window dependence, repeat annual values for each window
    if (length(tr)==nWIN*nYRS) {for (i in 1:nWIN){ for (j in 1:nYRS) out[,i,j]=tr[i,j]}}
    # If daily trend (seasonal cycle), compute running averages of tr
    if (length(tr)==length(ts)){
      for (iw in 1:nWIN){ w=WINDOWS[iw]
        trw=myrunavg(tr,w)
        for (id in 1:nDAY){ d=iDEVE[id]
          idays=which(DATES$m==DATES$m[d] & DATES$d==DATES$d[d])
          out[id,iw,]=trw[idays]
        }
      }
    }
  }
  out
}

#---------------------------------------------------------------------------------
# Function that returns the list of event attribution stats (p1, p0, pr).
# Inputs:
#   > ts and tr = vectors of data + trend
#     or input  = a list with pre-computed arrays dateve, datasp and datatr
#   > methodological parameters for computing stats with myperclevel()
#   > a few parameters for computing arrays with above functions if input=NULL
#
# Ouput:
#   < statistics, input data and fits
#
# Warning: this version uses a very simple detrending method (tr subtracted from ts)
#          which is relevant for temperatures but maybe not for other variables
#---------------------------------------------------------------------------------
compute.stat.1d=function(ts=NULL,tr=NULL,input=NULL,annual.maxima=FALSE,calendar.flex=0,
			 distrib="gauss",fixed.ksi=NA,smooth.fit.sdf=NA,exceedance="above",
			 first.year=TRUE,last.year=TRUE,block.day1=1,force.zero.trend=F){
  #---- tr must be of size ts or nYRS
  if (! length(tr) %in% c(length(ts),nYRS)) return("Error: tr size must be ts, nYRS or nWINxnYRS.")
  #---- fixed.ksi must be of size 1 or nWIN
  if (! length(fixed.ksi) %in% c(1,nWIN)) return("Error: fixed.ksi must be size 1 or nWIN.")
  #---- Building arrays datasp (samples), dateve (event values) and datatr (trends)
  if (!is.null(input)) {datasp=input$datasp; dateve=input$dateve; datatr=input$datatr}
  else {
    datasp=compute.array.of.data.samples(ts,annual.maxima,calendar.flex,exceedance,first.year,last.year,block.day1)
    dateve=compute.matrix.of.event.values(ts)
    datatr=compute.array.of.trend.values(ts,tr,annual.maxima,force.zero.trend)
  }
  #---- Computing arrays for detrended data wrt. t=YEVE (data1, for p1) and t=YREF (data0, for p0)
  data1=datasp-datatr+rep.abind(datatr[,,iYEVE],nYRS)
  data0=datasp-datatr+rep.abind(datatr[,,iYREF],nYRS)
  #---- Fitting parameters of the requested distribution
  if (distrib %in% c("gev","gevmin") & !is.na(fixed.ksi)){
    fit1=data1[,,1:3]; fit0=data0[,,1:3]
    ksi=rep(fixed.ksi,length(WINDOWS))[1:length(WINDOWS)]
    for (iw in 1:nWIN){ for (id in 1:nDAY){ if (any(!is.na(datasp[id,iw,]))){
      fit1[id,iw,]=myfitparams(data1[id,iw,],paste0(distrib,"fix"),ksi=ksi[iw])
      fit0[id,iw,]=myfitparams(data0[id,iw,],paste0(distrib,"fix"),ksi=ksi[iw])
    }}}
  }
  else {
    fit1=aperm(apply(data1,1:2,myfitparams,distrib),c(2,3,1))
    fit0=aperm(apply(data0,1:2,myfitparams,distrib),c(2,3,1))
  }
  #---- Smoothing fits if requested
  if (!is.na(smooth.fit.sdf)){
    # Calendar: smoothing over days
    if (!annual.maxima){
      fit1=apply(fit1,2:3,mysmoothy,sdf=smooth.fit.sdf,periodic=(nDAY==365))
      fit0=apply(fit0,2:3,mysmoothy,sdf=smooth.fit.sdf,periodic=(nDAY==365))
    }
    # Annual max: smoothing over windows
    # EDIT for GEV: take ksi as the mean over windows and re-fit mu and sigma with ksi constant
    if (annual.maxima){
      if (! distrib %in% c("gev","gevmin")){
        fit1[1,,]=apply(fit1[1,,],2,mysmoothy,sdf=smooth.fit.sdf)
        fit0[1,,]=apply(fit0[1,,],2,mysmoothy,sdf=smooth.fit.sdf)
      }
      if (distrib %in% c("gev","gevmin")){
        ksi=mean(fit1[1,,3],na.rm=T)
        for (iw in 1:nWIN){ for (id in 1:nDAY){ if (any(!is.na(datasp[id,iw,]))){
          fit1[id,iw,]=myfitparams(data1[id,iw,],paste0(distrib,"fix"),ksi)
          fit0[id,iw,]=myfitparams(data0[id,iw,],paste0(distrib,"fix"),ksi)
        }}}
      }
    }
  }
  #---- If annual.maxima, fits are repeated over all days
  if (annual.maxima){ for (id in 2:nDAY){ fit1[id,,]=fit1[1,,] ; fit0[id,,]=fit0[1,,] }}  
  #---- Computing p1 (percentile level of dateve within fit1)
  out.p1=abind(dateve,fit1,along=3)
  out.p1=apply(out.p1,1:2,function(x) myperclevel(x[1],distrib,x[-1],exceedance))
  out.p1=as.data.frame(out.p1); names(out.p1)=WINDOWS
  #---- Computing p0 (percentile level of dateve within fit0)
  out.p0=abind(dateve,fit0,along=3)
  out.p0=apply(out.p0,1:2,function(x) myperclevel(x[1],distrib,x[-1],exceedance))
  out.p0=as.data.frame(out.p0); names(out.p0)=WINDOWS
  #---- Computing probability ratio (warning: p0 and p1 are percentile levels here!)
  out.pr=(1-out.p1)/(1-out.p0)
  #---- Output
  list(p1=out.p1,p0=out.p0,pr=out.pr,
       data=datasp,value=dateve,trend=datatr,data1=data1,data0=data0,
       fit1=fit1,fit0=fit0,distrib=distrib)
}
