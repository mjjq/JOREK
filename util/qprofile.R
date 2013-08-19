# This is an R script for plotting the q-profile data
# in nice vector graphics.
# Can be run in the command line with:
# R CMD BATCH ./qprofile.R
# Written by Jane Pratt, August 2013. email: jane.pratt@ipp.mpg.de

#Uncomment the following line, and insert your working directory!
#setwd("/insert-working-directory-here/")

###### load q profile data
onam=paste("qprofile.dat",sep="")
all <- scan(file = onam, what = double(0), nmax = -1, dec = ".", skip= 0, nlines = 0, na.strings = "NA", flush = FALSE,
fill = FALSE, strip.white = FALSE, quiet = FALSE, blank.lines.skip = TRUE, multi.line = TRUE, comment.char = "")

p<-c()
q<-c()
r<-c()

elen=length(all)/3 -1
for(ii in c(0:elen)){
p<-append(p,all[3*ii+1],length(p))
q<-append(q,all[3*ii+2],length(q))
r<-append(r,all[3*ii+3],length(r))
}

## Chop off data outside the plasma radius for nice plotting.
## This is only necessary because qprofile.dat contains data outside
## the plasma, for which q=0.
## (This block can be commented out, as you prefer.)
psiq<-c()
qpro<-c()
rad<-c()
plen=length(p)
for(ii in c(1:plen)) {
pii=p[ii]
qii=q[ii]
rii=r[ii]
if(pii<1){
psiq<-append(psiq,pii,length(psiq))
qpro<-append(qpro,qii,length(qpro))
rad <-append(rad ,rii,length(rad))
}
}

########### plot q vs psi
fnam=paste("q_vs_psi.pdf",sep="")
pdf(fnam,width=10,height=8,family = "Times", pointsize=20)
par(cex = 1.0,mar=c(4.5,5,1,2))
tvec<-c(psiq)
pvec<-c(qpro)
plot(tvec,pvec,type="n", pch="o", cex=.3, cex.main=1.0, cex.lab=1.4, cex.axis=1.4,main="", xlab=expression(paste(psi[n])),ylab=expression(paste("q")), font=9, col = "purple", font=5, axes=FALSE, frame.plot=TRUE,lwd=1,xlim=c(0,1))

lines(psiq,qpro,type="l", pch=".", cex = 1.0, col = "deeppink",lty="solid",lwd=2)

### uncomment the following lines to draw lines to help pinpoint
### where the resonant surface is in psi.
#lines(c(-100,100),c(2.0,2.0),type="l", pch="o", cex = 0.3, col = "black",lty="solid",lwd=1)
#lines(c(.46,.46),c(-100,100),type="l", pch="o", cex = 0.3, col = "black",lty="solid",lwd=1)
#lines(c(.5,.5),c(-100,100),type="l", pch="o", cex = 0.3, col = "red",lty="dashed",lwd=1)
#lines(c(.6,.6),c(-100,100),type="l", pch="o", cex = 0.3, col = "black",lty="solid",lwd=1)

axis(1, tick=TRUE, lwd=1)
axis(2, tick=TRUE, lwd=1)
dev.off()


fnam=paste("q_vs_rad.pdf",sep="")
pdf(fnam,width=10,height=8,family = "Times", pointsize=20)
par(cex = 1.0,mar=c(4.5,5,1,2))
tvec<-c(rad)
pvec<-c(qpro)
plot(tvec,pvec,type="n", pch="o", cex=.3, cex.main=1.0, cex.lab=1.4, cex.axis=1.4,main="", xlab=expression(paste(r,"(m)")),ylab=expression(paste("q")), font=9, col = "purple", font=5, axes=FALSE, frame.plot=TRUE,lwd=1,xlim=c(0,1))

lines(rad,qpro,type="l", pch=".", cex = 1.0, col = "deeppink",lty="solid",lwd=2)

### uncomment the following lines to draw lines to help pinpoint 
### where the resonant surface is in radius.
#lines(c(.46,.46),c(-100,100),type="l", pch="o", cex = 0.3, col = "black",lty="solid",lwd=1)
#lines(c(.53,.53),c(-100,100),type="l", pch="o", cex = 0.3, col = "red",lty="dashed",lwd=1)
#lines(c(.6,.6),c(-100,100),type="l", pch="o", cex = 0.3, col = "black",lty="solid",lwd=1)
#lines(c(-100,100),c(2.0,2.0),type="l", pch="o", cex = 0.3, col = "black",lty="solid",lwd=1)

axis(1, tick=TRUE, lwd=1)
axis(2, tick=TRUE, lwd=1)
dev.off()


