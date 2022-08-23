######################################################################
## TMatch.R
## Calculate the temperature where Hauser Feshbach rates need to be
## used, then use them!
## 
## J. R. Newton et al., PHYSICAL REVIEW C 78, 025805 (2008)
######################################################################

RatesMCFile      <- "RatesMC.in"
ContributionFile <- "RatesMC.cont"
HFFile           <- "Tmatch.HF"

## The parameter, n, to determine maximum of Gamow peak
## (see Newton paper, Eqn. 3)
n.Gamow <- 1

######################################################################
## You shouldn't need to touch anything below this line...
library("RColorBrewer")
cols <- brewer.pal(4,"Dark2")

## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
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


## Make sure the RatesMC.in title matches the Tmatch.HF one
title <- scan(RatesMCFile,what=character(),n=1)
titleHF <- scan(HFFile,what=character(),n=1)

if(title != titleHF){
    stop("\n",
        "Input file for reaction '",title,"'\ndoesn't match Tmatch file ",
        "for '",titleHF,"'\n",sep="")
    ##stop()
}
  
## Read in the contribution file
data <- read.table(ContributionFile,skip=1,header=FALSE)

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

## Read the masses and charges
Z0 <- scan(RatesMCFile,what=numeric(),skip=2,n=1)
Z1 <- scan(RatesMCFile,what=numeric(),skip=3,n=1)
M0 <- scan(RatesMCFile,what=numeric(),skip=5,n=1)
M1 <- scan(RatesMCFile,what=numeric(),skip=6,n=1)

## Read in the resonance energies.
## This procedure should read in all resonance energies without erroring!
skip <- 29
Energies <- numeric()
while(
  ## Try to read a line and see if an actual number was read.
  substr(x <- try(unlist(read.table(RatesMCFile,skip=skip,nrows=1)[1]
                             )),1,2) != "**"){
  ## append the resonance energy list and skip to the next line
  skip <- skip+1
  if(is.numeric(x))
    Energies <- c(Energies,x)
}
## Ignore Upper limit resonances for this calculation?
## Now read the upper limit resonances in the same way
skip <- skip+5
while(
  substr(x <- try(unlist(read.table(RatesMCFile,skip=skip,nrows=1)[1]
                             )),1,2) != "**"){
  skip <- skip+1
  if(is.numeric(x))
    Energies <- c(Energies,x)
}
##Energies <- format(Energies,digits=1,trim=TRUE)
## The maximum energy
maxE <- max(Energies)


## Look for resonances that contribute more than 5% in the specified
## temperature range
## Look only at the centroids
MedSeq <- seq(from=9,to=length(data[1,]),by=3)
T <- data[,1]
resCont <- data[,MedSeq]

Gamow <- function(Temp, plot=FALSE, ymax=1){
    kT <- 0.086171*Temp
    somer <- 0.989534*Z0*Z1*sqrt(M0*M1/(M0+M1))

    ## Find the Gamow peak E0 and DE
    E0 = 0.1220*((Z0*Z1)^2 * (M0*M1/(M0+M1)) *Temp^2)^(1/3)
    DE = (4/sqrt(3))*sqrt(E0*kT)

    ## The gamow peak
    gamow <- function(x) exp(-somer/sqrt(x))*exp(-(x)/kT)
    if(plot){
	E<-seq(0,max(Energies)/1000,length.out=1000)
	maxG <- max(gamow(E))
	lines(E*1000,0.8*gamow(E)*ymax/maxG,col=cols[2])
    }
    
    ## Return the Gamow peak and width
    c(E0,DE)
}

## For each temperature, plot the contributions of each resonance as a bar graph
mypdf(file="GraphTMatch.pdf",width=8.25,height=10.75)
layout(matrix(c(1:9), 3, 3, byrow=TRUE))

## For each temperature,
## 1. Make a plot of the contributions
## 2. Calculate the Gamow Window
## 3. Calculate the Effective Energy Window
Gamow.array <- numeric()
ETER.array <- numeric()

for(i in 1:length(T)){

    plot(range(Energies),c(0,max(resCont[i,])),type='n',xlab="E",ylab="Cont",main=T[i])
    segments(x0=Energies,y0=0,y1=as.double(resCont[i,]),col=cols[1])
 
    ## Gamow window method first
    ## Gamow window returns E0 and DE
    G <- 1000*Gamow(T[i], plot=TRUE, ymax=max(resCont[i,]))
    Gamow.array <- c(Gamow.array,G[1]+n.Gamow*G[2])
    segments(x0=G[1]-n.Gamow*G[2], x1=G[1]+n.Gamow*G[2], y0=max(resCont[i,])*0.85,
	     col=cols[2],lwd=10,lend=1)
    text(x=G[1],y=max(resCont[i,])*0.85, "Gamow")

    ## Now do Joe's algorithm
    ## First the energies have to be in order, so sort them
    o <- order(Energies)
    Energies.o <- Energies[o]
    cont <- as.numeric(resCont[i,])[o]   ## These are the individual contributions

    ## Now do a linear interpolation. cubic spline acts too crazy!
    interp.ETER <- approxfun(y=Energies.o,x=cumsum(cont/sum(cont)))
    ETER <- interp.ETER(c(0.08,0.5,0.92))
    ETER.array <- c(ETER.array,ETER[2] + (ETER[3]-ETER[1]))

    segments(x0=ETER[1], x1=ETER[3], y0=max(resCont[i,])*0.92,
	     col=cols[1],lwd=10,lend=1)
    arrows(x0=ETER[2], x1=ETER[2] + (ETER[3]-ETER[1]), y0=max(resCont[i,])*0.92,
	     col=cols[1],lwd=2,length=0.1)
    text(x=(ETER[1]+ETER[3])/2,y=max(resCont[i,])*0.92, "ETER",adj=0.5)
    text(x=ETER[3],y=max(resCont[i,])*0.92, "Median + ETER",pos=3)
    ##plot(y=Energies.o,x=cumsum(cont/sum(cont)))
    ##abline(v=c(0.08,0.5,0.92))
    ##abline(h=interp.ETER(c(0.08,0.5,0.92)))
    
}    

## Now do the interpolation
## First Gamow peak method
interp.Gammow <- splinefun(x=Gamow.array, y=T)
TMatch.Gamow <- interp.Gammow(maxE)
cat("TMatch from Gamow peak method: ",TMatch.Gamow,"\n")

## Then ETER method
## Cut out low temperatures
cut <- T>0.1
interp.ETER <- splinefun(x=ETER.array[cut], y=T[cut])
TMatch.ETER <- interp.ETER(maxE)
cat("TMatch from ETER method:       ",TMatch.ETER,"\n")


dev.off()

## Next:
## Find the normalization needed for HF rate
## Normalize
## Format output

## Then come back and calculate narrow resonances

## END
######################################################################
