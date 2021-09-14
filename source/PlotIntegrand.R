while(dev.cur()>1)dev.off()
xlim <- NULL ##c(-1e-5,1e-5)
ylim <- NULL ##c(1e-22,1e2)

data <- readLines("test.dat")
breaks <- c(0,which(data=="!"))

for(i in 1:(length(breaks)-1)){
    myX11()
    par(mar=0.1+c(5,6,0.5,0.5))
    di <- data[(breaks[i]+1):(breaks[i+1]-1)]
    dat <- matrix(as.numeric(unlist(strsplit(di, " "))),ncol=2, byrow=TRUE)
    data2 <- dat[order(dat[,1]),]
    plot(data2, ylim=ylim, xlim=xlim, log="y", type='l',
         xlab="Energy (MeV)", ylab="Reaction rate integrand",
         col=add.alpha("black",1))

}
##myX11()
##    par(mar=0.1+c(5,6,0.5,0.5))
   
##    xlim <- NULL
##    ylim <- NULL
##data <- read.table("test.dat")
##data2 <- data[order(data[,1]),]
##plot(data2,type='l',log='y',xlim=xlim,ylim=ylim,pch=19,col=add.alpha("black",1),
##     xlab="Energy (MeV)", ylab="Reaction rate integrand")
