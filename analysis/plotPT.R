
dat <- read.table("RatesMC.PT")

i <- 0
j <- 0

myX11 <- function(...) 
{ 
    grDevices::X11(...) 
    par(cex.axis=1.3, cex.lab=1.5,   # Font sizes
        las=1,                       # Always horisontal text
        lwd=2,                       # Line width
        mar=c(5,5,3,2)+0.1,          # Margins
        pch=19,                      # Point type (solid circles)
        tcl=0.5,
        mgp=c(3,0.5,0)) 
}

while(dev.cur()>1)dev.off()
myX11(height=8,width=5)
par(mfrow=c(2,1))

cut1 <- dat[,2]==i
cut2 <- dat[,2]==j
h1 <- hist(log(dat[cut1,4]),plot=F,breaks=50)
h2 <- hist(log(dat[cut2,4]),plot=F,breaks=50)

plot(h1$mids,h1$density,xlab=paste("Resonance",i),ylab="Freq",
     type='s',xaxt='n')
aX <- axTicks(1)
axis(1,at=aX,labels=aX)
#,
#     xlim=range(c(h1$breaks,h2$breaks)),
#     ylim=range(c(h1$density,h2$density)))
abline(v=log(1e-3),col="red",lty=2)

plot(h2$mids,h2$density,xlab=paste("Resonance",j),ylab="Freq",
     type='s',xaxt='n')
aX <- axTicks(1)
axis(1,at=aX,labels=aX)
#,
#     xlim=range(c(h2$breaks,h2$breaks)),
#     ylim=range(c(h2$density,h2$density)))
##lines(h2$mids,h2$density,col="darkgreen",type='s')
abline(v=log(4.5e-3),col="red",lty=2)

cat("Res: ",i," PT = ",mean(dat[dat[,1]==i,3]),
    " +/- ",sd(dat[dat[,1]==i,3])," Fac. Unc. = ",
    sd(dat[dat[,1]==i,3])/mean(dat[dat[,1]==i,3])," \n")
cat("         median = ",median(dat[dat[,1]==i,3]),"\n")

cat("Res: ",j," PT = ",mean(dat[dat[,1]==j,3]),
    " +/- ",sd(dat[dat[,1]==j,3])," Fac. Unc. = ",
    sd(dat[dat[,1]==j,3])/mean(dat[dat[,1]==j,3])," \n")
cat("         median = ",median(dat[dat[,1]==j,3]),"\n")

##----------------------------------------------------------------------

myX11(height=8,width=5)
par(mfrow=c(2,1))

d1 <- dat[cut1,4]
d2 <- dat[cut2,4]
b <- 10^seq(from=log10(min(c(d1,d2))),to=log10(max(c(d1,d2))),length.out=1000)
h1 <- hist(d1,plot=F,breaks=b)
h2 <- hist(d2,plot=F,breaks=b)

plot(h1$mids,h1$density,xlab=paste("Resonance",i),ylab="Freq",
     type='s',xaxt='n',log='xy',xlim=c(1e-10,10),ylim=c(1e-5,1e10))
aX <- axTicks(1)
axis(1,at=aX,labels=aX)
#,
#     xlim=range(c(h1$breaks,h2$breaks)),
#     ylim=range(c(h1$density,h2$density)))
abline(v=4.5e-3,col="red",lty=2)
med <- median(d1)
abline(v=med,col="darkgreen",lty=2)

plot(h2$mids,h2$density,xlab=paste("Resonance",j),ylab="Freq",
     type='s',xaxt='n',log='xy',xlim=c(1e-10,10),ylim=c(1e-5,1e10))
aX <- axTicks(1)
axis(1,at=aX,labels=aX)
#,
#     xlim=range(c(h2$breaks,h2$breaks)),
#     ylim=range(c(h2$density,h2$density)))
##lines(h2$mids,h2$density,col="darkgreen",type='s')
abline(v=4.5e-3,col="red",lty=2)
med <- median(d2)
abline(v=med,col="darkgreen",lty=2)

cat("\n--------------------------------------------------\n")
cat("Res: ",i," PT = ",mean(dat[dat[,1]==i,4]),
    " +/- ",sd(dat[dat[,1]==i,4])," \n")
cat("         median = ",median(dat[dat[,1]==i,4]),"\n")

cat("Res: ",j," PT = ",mean(dat[dat[,1]==j,4]),
    " +/- ",sd(dat[dat[,1]==j,4])," \n")
cat("         median = ",median(dat[dat[,1]==j,4]),"\n")
