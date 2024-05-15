
#------------------------------------------
# Figure non-stationary normals
#------------------------------------------
fig.ns.norms=function(s,yrs1=1940:1941,yrs2=2022:2024,title=NULL){
  ts=c(s$data[,1,which(s$period %in% yrs1)],rep(NA,200),s$data[,1,which(s$period %in% yrs2)])
  tr=c(s$trend[,1,which(s$period %in% yrs1)],rep(NA,200),s$trend[,1,which(s$period %in% yrs2)])
  ylim=range(ts,na.rm=T); ylim=ylim+c(-0.1,0.1)*(ylim[2]-ylim[1])
  id=c(rep(1:365,length(yrs1)),rep(0,200),rep(1:365,length(yrs2)))
  matplot(cbind(ts,tr),type="l",lwd=c(1.5,2),lty=1,col=c(1,4),axes=F,frame.plot=T,ylim=ylim,
          main="",xlab="Days from 1 to 365 each year",ylab=var.name(s$var))
  polygon(range(which(id==0))[c(1,2,2,1)],c(-999,999)[c(1,1,2,2)],bg="white",border=NA,density=10)
  abline(v=unique(c(which(id==1)-0.5,which(id==365)+0.5)),lty=3)
  axis(1,at=which(id %in% c(1,90,180,270)),lab=rep(c(1,90,180,270),length(c(yrs1,yrs2))))
  axis(3,at=which(id %in% c(1,90,180,270)),lab=F); axis(2);axis(4,lab=F)
  text(which(id==182),ylim[1],c(yrs1,yrs2),font=2)
  legend("topleft",bg="white",box.lty=0,ncol=2,inset=0.02,
         leg=c(var.shortname(s$var),"Normals"),lwd=c(1.5,2),lty=1,col=c(1,4))
  if (is.null(title)) title=s$var
  title(title,adj=0)
}

#------------------------------------------
# Figure 2023 vs. recent years
#------------------------------------------
fig.comp.yrs=function(s,yrs=1991:2024,yr0=c(2016,2023,2024),ref.period=1991:2020,
                           cols0=c(2,4,"violet"),lty0=c(1,1,1),
                           yrnrm=c(1940,yr0),colsnrm=c(3,cols0),
                           title=NULL,force.anomaly=F){
  # Plot data
  ylim=range(s$data[,1,which(s$period %in% c(yrs,yr0))],na.rm=T)
  ylim=ylim+c(0,0.2)*(ylim[2]-ylim[1])
  matplot(s$data[,1,which(s$period %in% yrs)],type="l",lty=1,lwd=0.8,col="gray",
          main="",xlab="Days",ylab=var.name(s$var),axes=F,frame.plot=T,ylim=ylim,xlim=c(0,385))
  # Highlight yr0
  matplot(s$data[,1,which(s$period %in% yr0)],type="l",lty=lty0,col=cols0,lwd=2,add=T)
  # Get normals for ynrm + ref
  nrm=s$trend[,1,which(s$period %in% yrnrm)]
  ref=apply(s$trend[,1,which(s$period %in% ref.period)],1,mean)
  if (force.anomaly) ref=ref*0
  # Plot
  matplot(nrm,type="l",lty=5,col=colsnrm,lwd=2,add=T)
  lines(ref,col=1,lwd=2,lty=5)
  vals=round(apply(nrm,2,mean,na.rm=T),2)
  for (i in 1:length(vals)){ if (vals[i]>0) vals[i] = paste0("+",vals[i])}
  text(365,nrm[365,],paste0(vals,"°C"),pos=4,col=c(3,cols0),font=2,cex=0.9)
  # Axes
  my.calendar.axis(xlim=c(1,366)); my.calendar.axis(3,lab=F,xlim=c(1,366))
  axis(2,seq(-999,999,0.3)); axis(2,seq(-999,999,0.1),lab=F,tcl=0.25)
  # Legends
  legend("topleft",bg="white",box.lty=0,inset=0.02,ncol=2,
         leg=c(paste0(min(yrs),"-",max(yrs)),yr0),
         lty=c(1,lty0),lwd=c(1,2,2,2),col=c("gray",cols0))
  legend("topright",bg="white",box.lty=0,inset=0.02,ncol=2,
         leg=c(paste0("Ref. ",min(ref.period),"-",max(ref.period)),paste("Normals",yrnrm)),
         lty=5,lwd=1,col=c(1,colsnrm))
  # Title
  if (is.null(title)) title=s$var
  title(title,adj=0)
}

#------------------------------------------
# Figure yearly averages
#------------------------------------------
fig.yr.avg=function(s,yr0=2023,current.year=F,month=NULL,title=NULL){
  # Compute yearly averages
  # - complete year
  if (is.null(month)){
    ts=apply(s$data[,1,],2,mean,na.rm=T)
    fr=apply(s$trend[,1,],2,mean,na.rm=T)
  }
  # - individual month if specified
  if (!is.null(month)){
    nd=c(31,28,31,30,31,30,31,31,30,31,30,31)
    im1=cumsum(c(1,nd[-12]))
    im2=cumsum(nd)
    ts=apply(s$data[(im1[month]):(im2[month]),1,],2,mean,na.rm=T)
    fr=apply(s$trend[(im1[month]):(im2[month]),1,],2,mean,na.rm=T)
  }
  # Remove current year when incomplete
  if (!current.year) ts[length(ts)-1]=NA
  fr[c(1,length(fr))]=NA
  # Graphical params
  xlim=range(s$period); xlim=xlim+c(0,0.03)*(xlim[2]-xlim[1])
  ylim=range(c(ts,fr),na.rm=T); ylim=ylim+c(0,0.15)*(ylim[2]-ylim[1])
  # Plot
  plot(s$period,ts,main="",xlab="Years",ylab=var.name(s$var),
       axes=F,frame.plot=T,xlim=xlim,ylim=ylim)
  lines(s$period,fr,col=4); points(s$period,fr,pch=16,col=4)
  # Highlight yr0 value
  points(yr0,ts[which(s$period==yr0)],pch=21,bg=1,cex=1.2)
  text(yr0,ts[which(s$period==yr0)],yr0,cex=0.9,pos=4,font=2)
  # Highlight normal in yr0
  points(yr0,fr[which(s$period==yr0)],pch=21,bg=4,cex=1.5)
  text(yr0,ylim[1],paste(yr0,"normal"),srt=90,col=4,adj=c(0,1.5),font=2)
  # Highlight record previous to yr0
  rec=max(ts[which(s$period<yr0)],na.rm=T); irec=which(ts==rec)
  points(s$period[irec],rec,pch=21,bg=2,cex=1.5)
  text(s$period[irec],ylim[1],paste(s$period[irec],"record"),srt=90,col=2,adj=c(0,-0.5),font=2)
  # Add lines
  abline(h=c(rec,fr[which(s$period==yr0)]),lwd=1.5,col=c(2,4))
  abline(v=c(s$period[irec],yr0),lwd=1.5,col=c(2,4))
  # Axes
  ref=mean(ts[which(s$period %in% 1991:2020)])
  for (i in c(1,3)) axis(i,seq(1940,2030,10),tcl=0.5,lab=(i==1))
  for (i in c(1,3)) axis(i,seq(1940,2024,2),tcl=0.2,lab=F)
  axis(2,at=seq(10,20,0.5),tcl=0.5); axis(2,at=seq(10,20,0.1),tcl=0.25,lab=F)
  axis(4,at=ref+seq(-6,6,0.3),tcl=0.5,lab=seq(-6,6,0.3))
  axis(4,at=ref+seq(-6,6,0.1),tcl=0.25,lab=F)
  abline(h=ref,lty=3)
  # Legend
  legend("topleft",bg="white",box.lty=0,ncol=2,inset=0.02,
         leg=c(var.shortname(s$var),"Normals"),lty=c(0,1),pch=c(1,16),col=c(1,4))
  # Title
  if (is.null(title)) title=s$var
  title(title,adj=0)
}
