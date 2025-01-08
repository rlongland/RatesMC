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
xlim <- NA  ## Use NA for 0-1 MeV or enter range c(xmin,xmax) in MeV
ylim <- NA  ## Use NA for auto scaling or enter range c(1e-10,1e20)

## Draw the narrow resonances?
drawNarrow <- TRUE

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
mypdf <- function(file="output.pdf",...)
  {
    grDevices::pdf(file=file,...)
    par(cex.axis=1.3, cex.lab=1.5,   # Font sizes
        las=1,                       # Always horisontal text
        lwd=2,                       # Line width
        mar=c(5,5,3,2)+0.1,          # Margins
        pch=19,                      # Point type (solid circles)
        tcl=0.5,
        mgp=c(3,0.5,0))     
  }
axTexpr <- function(side, at = axTicks(side, axp=axp, usr=usr, log=log),
                    axp = NULL, usr = NULL, log = NULL, pad=FALSE,
                    pre=FALSE,dot=TRUE)
{
  ## Purpose: Do "a 10^k" labeling instead of "a e<k>"
  ##	      this auxiliary should return 'at' and 'label' (expression)
  ## -------------------------------------------------------------
  ## Arguments: as for axTicks()
  ##            pre determines if you should put the 3 in 3 \times 10^3
  ## --------------------------------------------------------------
  ## Author: Martin Maechler, Date:  7 May 2004, 18:01
  eT <- round(log10(abs(at)))# at == 0 case is dealt with below
  mT <- at / 10^eT 
  if(pad){
    eT[sign(eT)>-1] <- paste("+",eT[sign(eT)>-1],sep="")
    temp <- eT[nchar(tmp)<3]
    eT[nchar(eT)<3] <- paste(substr(temp,1,1),"0",substr(temp,2,2),sep="")
  }
  ss <- lapply(seq(along = at),
               function(i) if(at[i] == 0) quote(0) else
               if(pre){
                 if(dot)substitute(A%.%10^E, list(A=mT[i], E=eT[i])) else
                 substitute(A%*%10^E, list(A=mT[i], E=eT[i]))
               }else{
                 substitute(10^E, list(A=mT[i], E=eT[i]))
               })
   do.call("expression", ss) 
}


## Read and format the reaction name
ReacName = as.character(read.table(RatesMCFile,skip=0,header=FALSE,nrows=1)$V1)
# Make superscripts, subscripts etc
ReacName <- sub(",g\\)",",gamma\\)",ReacName)
ReacName <- sub("\\(g,","\\(gamma,",ReacName)
ReacName <- sub(",a\\)",",alpha\\)",ReacName)
ReacName <- sub("\\(a,","\\(alpha,",ReacName)
ReacName <- sub(",","*','*",ReacName)
##ReacName <- sub("\\(a*","\\(alpha",ReacName)
ReacName <- sub("\\)","\\)*",ReacName)
ReacName <- gsub("([[:digit:]]+)", "phantom()^{\\1}*", ReacName)

## Find the beginning of the resonances
skip <- 29
while(substr(x <- try(unlist(read.table(RatesMCFile,skip=skip,nrows=1)[1]
			     )),1,3) != "Ecm"){
				 skip <- skip+1
			     }

## Read in the resonance energies.
## This procedure should read in all resonance energies without erroring!
skip <- skip+1
Energies <- numeric()
while(
  ## Try to read a line and see if an actual number was read.
    substr(x <- try(unlist(read.table(RatesMCFile,skip=skip,nrows=1)[1]
			   )),1,2) != "**"){
			       ## append the resonance energy list and skip to the next line
			       skip <- skip+1
			       if(is.numeric(x))
				   Energies <- c(Energies,x)
			       if(x == "+")
				   Energies <- c(Energies,Energies[length(Energies)])
			   }
is.UL <- rep(0,length(Energies))

## Now read the upper limit resonances in the same way
skip <- skip+5
while(
  substr(x <- try(unlist(read.table(RatesMCFile,skip=skip,nrows=1)[1]
                             )),1,2) != "**"){
  skip <- skip+1
  if(is.numeric(x))
    Energies <- c(Energies,x)
  if(x == "+")
      Energies <- c(Energies,Energies[length(Energies)])
}
Energies <- format(Energies,digits=1,trim=TRUE)
is.UL <- c(is.UL,rep(1,length(Energies)-length(is.UL)))


## Read the proper s-factor file
data <- read.table("RatesMC.sfact",header=TRUE)

## Cut out the data to only include s-factors
is.good <- !apply(data,2,function(x)all(x == 0))

## Extract the narrow resonance energies
Energies.narrow <- as.double(Energies[!is.good[4:(3+length(Energies))]])/1000

## Cut the broad resonance data, Energies, and count
data <- data[,is.good]
is.UL.good <- is.UL[is.good[4:(3+length(Energies))]]
Energies <- Energies[is.good[4:(3+length(Energies))]]
nparts <- dim(data)[2]-1


## The list of colours
lc <- nparts
hue <- seq(from=0,to=1-1/lc,length.out=lc)
sat <- 0.8
val <- 0.7
col <- hsv(hue,sat,val,alpha=1)
icol <- 1

##while(dev.cur()>1)dev.off()
##myX11(width=8,height=6)
mypdf("GraphSFactor.pdf",width=8,height=6)
oldpar <- par(mar=0.1+c(5,5,1,1))

## First make the plotting pane
ylim.default <- range(data[,2:dim(data)[2]],na.rm=TRUE,finite=TRUE)
if(is.na(ylim[1]))ylim <- ylim.default
ylim[1] <- max(ylim[1],1e-10)
xlim.default <- c(0,1)
if(is.na(xlim[1]))xlim <- xlim.default



plot(xlim,ylim,ylim=ylim,type='n',log='y',
     xlab="Energy (MeV)", ylab="Astrophysical S-factor (MeV b)",
     yaxt='n',xaxs='i')
aY <- axTicks(2)
axis(2,at=aY,labels=axTexpr(2,at=aY))


## Draw all narrow resonances
if(drawNarrow){
    abline(v=Energies.narrow,lty=2,col="grey")
    yy <- grconvertY(0.9,from="npc",to="user")
    text(x=Energies.narrow,y=yy,
	 labels=paste(Energies.narrow*1000,"keV"),srt=90,pos=3)
}


## Draw all broad resonance
for(i in 1:nparts){
    lines(data[,1],data[,i+1],col=col[icol])
    icol <- icol+1
}

aRate.lab <- c("Non-res-1","Non-res-2")[is.good[2:3]]
ULstring <- ifelse(is.UL.good,"(UL)","")

## Make the interference labels
iInter <- grep("Int",names(data))
InterNumber <- as.double(substr(names(data)[iInter],4,6))
Inter.lab <- paste("Intf",InterNumber)

counter <- rep("",length(Energies))
if(length(Energies)>1){
    c <- 1
    for(i in 2:length(Energies)){
	if(Energies[i] == Energies[i-1]){
	    c <- c+1
	    counter[i] <- paste("#",c,sep="")
	} else {
	    c <- 1
	}
    }
}

leg <- aRate.lab
if(length(Energies)>0)
    leg <- c(aRate.lab,
	     ##paste("Res",1:(nparts-2)))
	     paste(Energies,"keV",counter,ULstring))
if(length(iInter)>0)
    leg <- c(leg,Inter.lab)

if(drawNarrow){
    leg <- c(leg,"Narrow res.")
    legend(x="topright",legend=leg,
	   col=c(col,"grey"),lty=c(rep(1,length(col)),2),
	   bg="white")
} else {
    legend(x="topright",legend=leg,
	   col=col,lty=1,
	   bg="white")
}

box()

dev.off()
##----------------------------------------------------------------------
