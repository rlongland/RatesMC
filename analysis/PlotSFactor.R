######################################################################
## PlotSFactor.R
## Plot the broad resonance s-factor produced by the integration
## routines
##
## Need to add vertical lines for narrow resonances
######################################################################

RatesMCFile   <- "RatesMC.in"
IntegrandFile <- "RatesMC.integ"
SFactorFile   <- "RatesMC.sfact"

## Set limits to NA for auto plotting
xlim <- c(0,1.5)
ylim <- c(1e-10,1e15)

######################################################################
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
## Read the proper s-factor file
data <- read.table("RatesMC.sfact",header=TRUE)

nparts <- dim(data)[2]-1

## The list of colours
lc <- nparts
hue <- c(seq(from=0,to=1-2/lc,length.out=floor(lc/2)),
         seq(from=1/lc,to=1-1/lc,length.out=ceiling(lc/2)))
sat <- 0.8
val <- 0.7
col <- hsv(hue,sat,val,alpha=1)
icol <- 1

while(dev.cur()>1)dev.off()
myX11()


## First make the plotting pane
ylim.default <- range(data[,2:dim(data)[2]])
if(is.na(ylim[1]))ylim <- ylim.default
xlim.default <- c(0,10)
if(is.na(xlim[1]))xlim <- xlim.default
 
plot(xlim,ylim,ylim=ylim,type='n',log='y',
     xlab="E (MeV)", ylab="SFactor (MeV b)")
 
for(i in 1:nparts){

    lines(data[,1],data[,i+1],col=col[icol])
    icol <- icol+1
}

leg <- c("Non-res-1",
	 "Non-res-2",
	 paste("Res",1:(nparts-2)))

legend(x="topright",legend=leg,col=col,lty=1)


## Pick a temperature. This doesn't affect the astrophysical s-factor,
## but does affect the integration routines that produce the data
## file. Thus, pick a good representative energy
##Temp <- 0.3

##   ## Read all lines to look for temperature
##   lines <- readLines(IntegrandFile)
##   TLines <- grep("Temp",lines)
##   ## Which temperature entry is the one we want 
##   Tentry <- which(read.table(text=lines[TLines])[,3] == Temp)
##   ## and corresponds to line:
##   goodTLine <- TLines[Tentry]
##   ## The next temperature entry is
##   nextTLine <- TLines[Tentry+1]
##   
##   if(nextTLine - goodTLine == 2){
##       stop("There is no broad resonances!!")
##   }
##   
##   ## The list of lines with integrand data
##   ilines <- (goodTLine+1):(nextTLine-2)
##   
##   ## Now each resonance is separated by spaces
##   isep <- c(1,grep("^$",lines[ilines]))
##   ## Now we can read the integrand for each resonance we encounter
##   integ <- lapply(2:length(isep), function(x){
##       data <- read.table(text=lines[ilines[isep[x - 1]:(isep[x]-1)]])
##       data
##   })
##   
##   ## Close any open windows
##   while(dev.cur()>1)dev.off()
##   myX11()
##   
##   ## First make the plotting pane
##   ylim.default <- range(unlist(lapply(integ, function(x)range(x[,3]))))
##   if(is.na(ylim[1]))ylim <- ylim.default
##   xlim.default <- c(0,10)
##   if(is.na(xlim[1]))xlim <- xlim.default
##   
##   plot(xlim,ylim,ylim=ylim,type='n',log='y',
##        xlab="E (MeV)", ylab="SFactor (arb. units)")
##   
##   ## The list of colours
##   lc <- length(integ)
##   hue <- c(seq(from=0,to=1-2/lc,length.out=floor(lc/2)),
##            seq(from=1/lc,to=1-1/lc,length.out=ceiling(lc/2)))
##   sat <- 0.8
##   val <- 0.7
##   col <- hsv(hue,sat,val,alpha=1)
##   icol <- 1
##   
##   ## Plot each line
##   t <- lapply(integ, function(data){
##       ## Sort data by energy
##       ord <- order(data[,1])
##       T9 <- data[ord,1]
##       y <- data[ord,3]
##       ## Finally plot
##       lines(T9,y,col=col[icol])
##       icol <<- icol+1
##       y
##   })
##   
##   legend(x="topright",legend=paste("Res",1:length(integ)),col=col,lty=1)


##----------------------------------------------------------------------
