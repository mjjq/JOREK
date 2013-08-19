# This is an R script for plotting (in nice vector graphics)
# the value of psi at a given point vs time data.
# The raw data from JOREK is processed with the post-processing program
# diagnostics/jordel.f90, and this R script plots the data output.

# Written by Jane Pratt, Feb 2013. email: jane.pratt@ipp.mpg.de

# R scripts can be run in the command line with:
# R CMD BATCH ./psifreqplot.R
# Running R in batch this way automatically produces an error file
# named psifreqplot.Rout.

# Uncomment, and put the appropriate directory in the following line:
#setwd("/insert-your-home-directory-here/")

###############
## import time and global energy data
###############
## 
## Time data is taken from macroscopic_vars.dat, produced with the command:
## ./extract_live_data.sh energies energies.dat
## Strip the first line (the characters) from energies.dat before plotting
## note that this data input is for n_tor=7, for different n_tor, 
## the number of energy vectory needs to be changed

onam=paste("energies.dat",sep="")
all <- scan(file = onam, what = double(0), nmax = -1, dec = ".", skip= 0, nlines = 0, na.strings = "NA", flush = FALSE,
fill = FALSE, strip.white = FALSE, quiet = FALSE, blank.lines.skip = TRUE, multi.line = TRUE, comment.char = "")

#array for time data
dt<-c()

#arrays for magnetic energies
em00<-c()
em01<-c()
em02<-c()
em03<-c()

#arrays for kinetic energies
ek00<-c()
ek01<-c()
ek02<-c()
ek03<-c()

elen=length(all)/9 -1
for(ii in c(0:elen)){
dt  <-append(dt  ,all[9*ii+1],length(dt))
em00<-append(em00,all[9*ii+2],length(em00))
em01<-append(em01,all[9*ii+3],length(em01))
em02<-append(em02,all[9*ii+4],length(em02))
em03<-append(em03,all[9*ii+5],length(em03))
ek00<-append(ek00,all[9*ii+6],length(ek00))
ek01<-append(ek01,all[9*ii+7],length(ek01))
ek02<-append(ek02,all[9*ii+8],length(ek02))
ek03<-append(ek03,all[9*ii+9],length(ek03))
}

########################
## load psi data produced with the jordel program
#########################

psi<-c()
s<-c()
t<-c()
dexs<-c()

# insert in sel the number of files from deltas that you want to plot
sel <-c(1:1912)
for(ii in sel){
lnum=paste("0000",ii,sep="")
if(ii>9) lnum=paste("000",ii,sep="")
if(ii>99)lnum=paste("00",ii,sep="")
if(ii>999)lnum=paste("0",ii,sep="")
if(ii>9999)lnum=paste("",ii,sep="")

onam=paste("deltas/dr",lnum,".jdat",sep="")
dr <- scan(file = onam, what = double(0), nmax = -1, dec = ".", skip= 0, nlines = 0, na.strings = "NA", flush = FALSE, fill = FALSE, strip.white = FALSE, quiet = FALSE, blank.lines.skip = TRUE, multi.line = TRUE, comment.char = "")

onam=paste("deltas/dz",lnum,".jdat",sep="")
dz <- scan(file = onam, what = double(0), nmax = -1, dec = ".", skip= 0, nlines = 0, na.strings = "NA", flush = FALSE, fill = FALSE, strip.white = FALSE, quiet = FALSE, blank.lines.skip = TRUE, multi.line = TRUE, comment.char = "")

onam=paste("deltas/psi",lnum,".jdat",sep="")
dp <- scan(file = onam, what = double(0), nmax = -1, dec = ".", skip= 0, nlines = 0, na.strings = "NA", flush = FALSE, fill = FALSE, strip.white = FALSE, quiet = FALSE, blank.lines.skip = TRUE, multi.line = TRUE, comment.char = "")

oa<-order(dr)
dr<-dr[oa]
dz<-dz[oa]
dp<-dp[oa]

ds<-sqrt(dr^2+dz^2)
sdex=1
for(jj in c(1:length(ds)) ) {
# choose a point along the line of data that jordel produces
# psi at this point will be plotted in each time step
if(ds[jj]>.501 && ds[jj]<.505) sdex=jj
}

ds1=dr[sdex]

dexs<-append(dexs,sdex,length(dexs))
s<-append(s,ds1,length(s))
psi<-append(psi,dp[sdex],length(psi))

} # ii loop

psi3<-psi/max(psi)
t3<-dt[sel]
# This last line may need to be altered, depending on what time steps
# you actually want to plot.

#######################
## produce plot
#######################

fnam=paste("psifreqplot.pdf",sep="")
pdf(fnam,width=10,height=8,family = "Times", pointsize=20)
par(cex = 1.0,mar=c(4.5,5,1,2))
tvec<-c(t3)
pvec<-c(psi3)
plot(tvec,pvec,type="n", pch="o", cex=.3, cex.main=1.0, cex.lab=1.4, cex.axis=1.4,main=expression(paste("")), xlab=expression(paste(t,"(jorek units)")),ylab=expression(paste(psi[n=1], "(jorek units)")), font=9, col = "purple", font=5, axes=FALSE, frame.plot=TRUE,lwd=1)

lines(t3,psi3,type="l", pch=".", cex = 1.0, col = "green4",lty="solid",lwd=2)
lines(t3,psi3,type="p", pch="o", cex = 0.3, col = "green4",lty="solid",lwd=1)

axis(1, tick=TRUE, lwd=1)
axis(2, tick=TRUE, lwd=1)
dev.off()


