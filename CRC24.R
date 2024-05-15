setwd("/home/cattiaux/Articles/SCAN2021/R")

source("scan_source.R")
source("scan_figures.R")
source("scan_functions.R")
source("CRC24_figures.R")

figdir="../article_tglo/NCC/figures/"
tabdir="../article_tglo/NCC/tables/"

# Get data
a=get.scan.1d("GSAT",loc="Global",set="tglo")
b=get.scan.1d("GSAT-anoref",loc="Global",set="tglo")
c=get.scan.1d("GSAT-ano",loc="Global",set="tglo")
d=get.scan.1d("GSAT-ano",meth="locmax15",loc="Global",set="tglo")
e=get.scan.1d("GSAT-ano",meth="annmax3",loc="Global",set="tglo")
f=get.scan.1d("GSAT-ano",meth="annmax3",loc="Global",set="tglo_85")

# Figure 1
pdf(paste0(figdir,"fig1.pdf"),6.5,9.5)
layout(matrix(1:3,3,1))
par(mar=c(2.5,2.5,1.5,1.5),mgp=c(1.2,0.1,0),las=0,tcl=0.5)
fig.ns.norms(a,title="a) ERA5 daily global surface air temperature")
fig.yr.avg(a,title="b) ERA5 yearly global surface air temperature")
fig.comp.yrs(b,yrnrm=c(1940,2016,2023),title="c) ERA5 daily global surface air temperature | Anomalies wrt. 1991-2020",force.anomaly=T)
dev.off()

# For the text
mean(b$trend[,1,which(b$period%in%2023)])  # normal in 2023 vs 1991-2020 - yearly average
mean(b$trend[,1,which(b$period%in%2022)])

range(b$trend[,1,which(b$period%in%2023)]) # normal in 2023 vs 1991-2020 - min/max values
which.max(b$trend[,1,which(b$period%in%2023)])
which.min(b$trend[,1,which(b$period%in%2023)])

cbind(a$period,apply(a$data[,1,],2,mean,na.rm=T))[which(a$period>1991),]  # yearly data
cbind(b$period,apply(b$data[,1,],2,mean,na.rm=T))[which(b$period>1991),]  # yearly anoref
cbind(c$period,apply(c$data[,1,],2,mean,na.rm=T))[which(c$period>1991),]  # yearly ano (top10 below)
cbind(c$period,apply(c$data[,1,],2,mean,na.rm=T))[order(apply(c$data[,1,],2,mean,na.rm=T),decreasing=T)[1:10],]

# Figure 2
# Max events 2023/24 to add on the plot
dx=get.max.events(d,years=2023:2024,only.max=T)$events
ex=get.max.events(e,years=2023:2024,only.max=T)$events
# A function to add an event on myfig.eve
myarrow=function(x,y,col,txt){
  segments(x[1],y,x[2],y,col=col)
  text(x[1],y,"<",adj=c(0,0.3),col=col)
  text(x[2],y,">",adj=c(1,0.3),col=col)
  text(mean(x),y,txt,adj=c(0.5,-0.7),col=col,cex=0.9)
}
# Panel plot
#pdf(paste0(figdir,"fig2.pdf"),7.5,7.5)
#layout(rbind(c(1,4),c(2,3)),widths=c(1,0.2))
pdf(paste0(figdir,"fig2.pdf"),6.5,9.5)
layout(rbind(c(1,6),c(2,3),c(4,5)),widths=c(1,0.2))
myfig.eve(c,y=2023:2024,xlim=c(1,456),id=dx$id,iw=dx$iw,title=NULL,force.anom=T,margins=F,rtitle=F)
myarrow(c(dx$day1,dx$day2)+365,-0.4,"darkviolet",dates2leg(dx$date1,dx$date2,short=T))
abline(v=c(ex$day1,ex$day2)+365)
myarrow(c(ex$day1,ex$day2)+365,-0.55,"black",dates2leg(ex$date1,ex$date2,short=T))
title("a) GSAT anomaly 2023/24",adj=0)
myfig.1d(d,y=2023:2024,xlim=c(1,456),title=NULL,rtitle=F,colorbar=F,lev=BRKSM,show.events=3)
title("b) GSAT anomaly 2023/24 | Calendar p1",adj=0)
myfig.1d(d,image=F,colorbar=T,double.cb=F)
myfig.1d(e,y=2023:2024,xlim=c(1,456),title=NULL,rtitle=F,colorbar=F,lev=BRKSM)
title("c) GSAT anomaly 2023/24 | All-year p1",adj=0)
myfig.1d(e,image=F,colorbar=T,double.cb=F)
dev.off()

1-e$scan$Y2023$p1[dx$id,dx$iw]
1/(1-e$scan$Y2023$p1[dx$id,dx$iw])

dxx=get.max.events(d,years=2023:2024)$events[1:3,]
1/(1-e$scan$Y2023$p1[dxx$id[2],dxx$iw[2]])
1/(1-e$scan$Y2023$p1[dxx$id[3],dxx$iw[3]])

1-ex$p1
1/(1-ex$p1)

1/(1-f$scan$Y2023$p1[dx$id,dx$iw])
1/(1-f$scan$Y2023$p1[ex$id,ex$iw])

# Table 1
hot1=get.max.events(d,min.value=0.95)
tab1=events2tab.anom(hot1,10,p1=T)
write.table(tab1,paste0(tabdir,"tab1.tex"),quote=F,row.names=F,col.names=F)

# Table 2
hot2=get.max.events(e,min.value=0.75)
tab2=events2tab.anom(hot2,10,p1=T)
write.table(tab2,paste0(tabdir,"tab2.tex"),quote=F,row.names=F,col.names=F)

# Supplementary

for (y in c(1983,1995,1998,2016)){
  yleg=paste0(y-1,"/",mymod(y,100))
  neve=1; if (y==1998) neve=4
  pdf(paste0(figdir,"suppfig2_",y,".pdf"),9.5,7.5)
  layout(rbind(c(1,2),c(3,4)),widths=c(1,0.2))
  myfig.1d(d,y=(y-1):y,xlim=c(1,730),title=NULL,rtitle=F,colorbar=F,lev=BRKSM,show.events=neve)
  title(paste("a) GSAT anomaly",yleg,"| Calendar p1"),adj=0)
  myfig.1d(d,image=F,colorbar=T,double.cb=F)
  myfig.1d(e,y=(y-1):y,xlim=c(1,730),title=NULL,rtitle=F,colorbar=F,lev=BRKSM)
  title(paste("b) GSAT anomaly",yleg,"| All-year p1"),adj=0)
  myfig.1d(e,image=F,colorbar=T,double.cb=F)
  dev.off()
}

for (scen in c("SSP1-2.6","SSP2-4.5","SSP5-8.5")){
  exp=paste0(substr(scen,6,6),substr(scen,8,8))
  aa=get.scan.1d("GSAT",loc="Global",set=paste0("tglo_",exp))
  bb=get.scan.1d("GSAT-anoref",loc="Global",set=paste0("tglo_",exp))
  pdf(paste0(figdir,"suppfig1_",exp,".pdf"),9,9)
  layout(matrix(1:2,2,1))
  par(mar=c(2.5,2.5,1.5,1.5),mgp=c(1.2,0.1,0),las=0,tcl=0.5)
  fig.yr.avg(aa,title=paste("a) ERA5 yearly GSAT | Normals with",scen))
  fig.comp.yrs(bb,yrnrm=c(1940,2016,2023),title=paste("b) ERA5 daily GSAT | Anomalies wrt. 1991-2020 | Normals with",scen),force.anomaly=T)
  dev.off()
}
