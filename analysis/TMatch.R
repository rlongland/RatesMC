######################################################################
## TMatch.R
## Calculate the temperature where Hauser Feshbach rates need to be
## used, then use them!
## 
## J. R. Newton et al., PHYSICAL REVIEW C 78, 025805 (2008)
######################################################################

RatesMCFile      <- "RatesMC.in"
RatesMCOutput    <- "RatesMC.out"
ContributionFile <- "RatesMC.cont"
HFFile           <- "Tmatch.HF"

OutputFile       <- "TMatch.out"
LaTeXFile        <- "TMatch.latex"

## Read the contributions file? (otherwise use narrow resonances)
ReadCont         <- FALSE

## The parameter, n, to determine maximum of Gamow peak
## (see Newton paper, Eqn. 3)
n.Gamow          <- 1

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
title <- scan(RatesMCFile,what=character(),n=1,quiet=TRUE)
titleHF <- scan(HFFile,what=character(),n=1,quiet=TRUE)

if(title != titleHF){
    stop("\n",
        "Input file for reaction '",title,"'\ndoesn't match Tmatch file ",
        "for '",titleHF,"'\n",sep="")
    ##stop()
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

## Read the masses and charges
Z0 <- scan(RatesMCFile,what=numeric(),skip=2,n=1,quiet=TRUE)
Z1 <- scan(RatesMCFile,what=numeric(),skip=3,n=1,quiet=TRUE)
M0 <- scan(RatesMCFile,what=numeric(),skip=5,n=1,quiet=TRUE)
M1 <- scan(RatesMCFile,what=numeric(),skip=6,n=1,quiet=TRUE)
J0 <- scan(RatesMCFile,what=numeric(),skip=8,n=1,quiet=TRUE)
J1 <- scan(RatesMCFile,what=numeric(),skip=9,n=1,quiet=TRUE)
mu <- M0*M1/(M0+M1)

## The RatesMC rates
Rates <- read.table(RatesMCOutput,skip=4)
T <- Rates[,1]


if(ReadCont){
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

    ## Read in the contribution file
    data <- read.table(ContributionFile,skip=1,header=FALSE)

    ## Look for resonances that contribute more than 5% in the specified
    ## temperature range
    ## Look only at the centroids
    MedSeq <- seq(from=9,to=length(data[1,]),by=3)
    resCont <- data[,MedSeq]

} else {

    ## Read all resonances from the RatesMC input file
    data <- numeric()
    skip <- 30
    ## Keep looping until we're done
    while(TRUE){
	## Read a line
	x <- try(read.table(RatesMCFile,skip=skip,nrows=1))
	## Break the loop if we see *s
	if(substr(x[,1],1,2) == "**") break
	## Is the first entry a '+' or '!'? Ignore if so
	if(x[1] == "+" || x[1] == "!" || sign(x[1]) == -1) {
	    skip <- skip+1
	    next
	}
	## Otherwise pull the information
	data <- rbind(data,x[c(1,3,5,6,9,12)])

	## Move onto next line and start again!
	skip <- skip+1
    }

    ## The energies and resonance strengths
    Energies <- data[,1]
    ## If resonance strength is unknown, calculate it from the given partial widths
    wgUnknown <- data[,2]==0
    w <- (2*data[,3]+1)/((2*J0 + 1)*(2*J1 + 1))
    G1 <- data[,4]
    G2 <- data[,5]
    GT <- G1+G2+data[,6]
    wg <- w*G1*G2/GT
    ## Then back-fill the unknown values
    data[wgUnknown,2] <- wg[wgUnknown]
    ## And grab it all back into a convenient vector
    wg <- data[,2]

    ## Now that we know the resonance energies and strengths, the rate
    ## can be calculated for each resonance at each temperature
    allRate <- sapply(T,function(Ti)(1.5399e11/(mu*Ti)^(1.5))*wg*1e-6*exp(-0.011605*Energies/Ti))
    resCont <- t(apply(allRate,2,function(t)t/sum(t)))

    ##resCont[is.nan(resCont)] <- NA
}

## The maximum energy
maxE <- max(Energies)

## A function to calculate the Gamow peak and width, and plot it if desired!
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
def.par <- par(no.readonly = TRUE)

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
    csum <- cumsum(cont/sum(cont))
    if(length(unique(csum)) > 1){
	interp.ETER <- approxfun(y=Energies.o,x=cumsum(cont/sum(cont)))
	ETER <- interp.ETER(c(0.08,0.5,0.92))
	ETER.array <- c(ETER.array,ETER[2] + (ETER[3]-ETER[1]))

	segments(x0=ETER[1], x1=ETER[3], y0=max(resCont[i,])*0.92,
		 col=cols[1],lwd=10,lend=1)
	arrows(x0=ETER[2], x1=ETER[2] + (ETER[3]-ETER[1]), y0=max(resCont[i,])*0.92,
	       col=cols[1],lwd=2,length=0.1)
	text(x=(ETER[1]+ETER[3])/2,y=max(resCont[i,])*0.92, "ETER",adj=0.5)
	text(x=ETER[3],y=max(resCont[i,])*0.92, "Median + ETER",pos=3)
	##plot(y=Energies.o,x=cumsum(cont/sum(cont)),xlim=c(0,1))
	##abline(v=c(0.08,0.5,0.92))
	##abline(h=interp.ETER(c(0.08,0.5,0.92)))
    } else {
	ETER.array <- c(ETER.array,NA)
    }
    
}    

##----------------------------------------------------------------------
## Now do the interpolation to find the temperature where maximum energy resonance corresponds
## to the end of the Gamow peak or effective burning energy
cat("\n")

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


######################################################################
## Now match the HF rate
HF <- read.table(HFFile, skip=3)

## Make an interpolation spline for Rates and HF Rate
interp.Rate.med  <- splinefun(Rates[,1],Rates[,3])
interp.Rate.low  <- splinefun(Rates[,1],Rates[,2])
interp.Rate.high <- splinefun(Rates[,1],Rates[,4])
interp.HF <- splinefun(HF)

Rates.MatchT.med  <- interp.Rate.med(TMatch.ETER)
Rates.MatchT.low  <- interp.Rate.low(TMatch.ETER)
Rates.MatchT.high <- interp.Rate.high(TMatch.ETER)
HF.MatchT <-    interp.HF(TMatch.ETER)
## Normalization (multiplicative) to apply to HF rate
Norm.med  <- Rates.MatchT.med/HF.MatchT
Norm.low  <- Rates.MatchT.low/HF.MatchT
Norm.high <- Rates.MatchT.high/HF.MatchT

## Now the matched HF rate!
HF.matched.med  <- HF[,2]*Norm.med
HF.matched.low  <- HF[,2]*Norm.low
HF.matched.high <- HF[,2]*Norm.high
HF <- cbind(HF, HF.matched.low, HF.matched.med, HF.matched.high)

## We need to interpolate the matched HF rate at the RatesMC temperatures
interp.HF.low  <- splinefun(HF[,1], HF[,3])
interp.HF.med  <- splinefun(HF[,1], HF[,4])
interp.HF.high <- splinefun(HF[,1], HF[,5])
HF.final.low  <- interp.HF.low(Rates[,1])
HF.final.med  <- interp.HF.med(Rates[,1])
HF.final.high <- interp.HF.high(Rates[,1])

## Plot the HF and, scaled HF, and interpolated HF to check it all went well!
layout(matrix(c(1:2), 2, 1, byrow=TRUE))
ylim <- range(c(Rates[,2],HF[,2:5]))

plot(Rates[,1], Rates[,3], type='l', ylim=ylim,
     xlab="Temperature", ylab="Rate", col=cols[1], log="",
     xaxs='i',main=paste("TMatch = ",format(TMatch.ETER,digits=2,nsmall=2),"GK"))
abline(v=TMatch.ETER,lty=3,col="grey")
lines(HF[,1],HF[,2],col=cols[2],lty=2)
lines(HF[,1],HF[,4],col=cols[3],lty=1)
lines(Rates[,1],HF.final.med,col=cols[4],lty=2)

legend("topleft",legend=c("Rec.","HF","Matched HF", "Matched HF at RatesMC gridpoints"),
       lty=c(1,2,1,2),col=cols[c(1,2,3,4)])

##----------------------------------------------------------------------
## Create the final rates table
## Logical vectors for above and below TMatch
cut <- T<TMatch.ETER
cut.length <- sum(cut)
acut <- T >= TMatch.ETER
acut.length <- sum(acut)

Rates.final <- Rates[cut,]
Rates.final <- rbind(Rates.final,
		     cbind(Rates[acut,1],
			   HF.final.low[acut], HF.final.med[acut], HF.final.high[acut],
			   rep(Rates[cut.length,5],acut.length)))

## Plot the final rate in a nice simple figure
ylim <- range(Rates.final[,2:4])
plot(Rates[,1], Rates[,3], type='l', ylim=ylim,
     xlab="Temperature", ylab="Rate", col=cols[1], log="",
     xaxs='i')
abline(v=TMatch.ETER,lty=3,col="grey")
lines(Rates.final[,1], Rates.final[,2], lwd=1, col=cols[3])
lines(Rates.final[,1], Rates.final[,3], lwd=2, col=cols[3])
lines(Rates.final[,1], Rates.final[,4], lwd=1, col=cols[3])
legend("topleft",legend=c("Unmatched Rate","Matched Rate"),
       lty=1,col=cols[c(1,3)])

dev.off()

######################################################################
## Finally write the output!

## Read in the header from RatesMC.out
header <- readLines(RatesMCOutput,n=4)
writeLines(header, OutputFile)

## Output the Rates up to the matching temperature
cut <- T<TMatch.ETER
OP.T <- formatC(Rates.final[cut,1],"f",width=6, digits=3)
OP.Low  <- formatC(Rates.final[cut,2],"E",width=12, digits=3,flag="#")
OP.Med  <- formatC(Rates.final[cut,3],"E",width=16, digits=3,flag="#")
OP.High <- formatC(Rates.final[cut,4],"E",width=16, digits=3,flag="#")
OP.fu   <- formatC(Rates.final[cut,5],"E",width=17, digits=3,flag="#")
OP.below <- cbind(OP.T, OP.Low, OP.Med, OP.High, OP.fu)
LaTeX.below <- cbind(OP.T, " & ", OP.Low, "  & ", OP.Med, "  &\n", "     ",
		     OP.High, "  &", OP.fu, "    \\\\")
write.table(OP.below,OutputFile,col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE, sep="")
write.table(LaTeX.below,LaTeXFile,col.names=FALSE, row.names=FALSE, quote=FALSE, append=FALSE, sep="")

## Repeat the process above the matching temperature
cut <- T>=TMatch.ETER
OP.T <- formatC(Rates.final[cut,1],"f",width=6, digits=3)
OP.Low  <- formatC(Rates.final[cut,2],"E",width=8, digits=3,flag="#")
OP.Med  <- formatC(Rates.final[cut,3],"E",width=8, digits=3,flag="#")
OP.High <- formatC(Rates.final[cut,4],"E",width=8, digits=3,flag="#")
OP.fu   <- formatC(Rates.final[cut,5],"E",width=8, digits=3,flag="#")
OP.above <- cbind(OP.T, "  (", OP.Low, ")", "       (", OP.Med,")", "     (", OP.High, ")","      (", OP.fu,")")
LaTeX.above <- cbind(OP.T, " &   (", OP.Low, ") &       (", OP.Med, ") &\n", "           ",
		     "(",OP.High, ") &       (", OP.fu, ")   \\\\")
write.table(OP.above,OutputFile,col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE, sep="")
write.table(LaTeX.above,LaTeXFile,col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE, sep="")

## Next:
## calculate narrow resonances
## 
## 

## END
######################################################################
