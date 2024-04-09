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
HFFile           <- "TMatch.HF"

massfile         <- "mass_1.mas20"
nubasefile       <- "nubase_3.mas20"

OutputFile       <- "TMatch.out"
LaTeXFile        <- "TMatch.latex"

## Read the contributions file? (otherwise use narrow resonances)
ReadCont         <- FALSE

## The parameter, n, to determine maximum of Gamow peak
## (see Newton paper, Eqn. 3)
n.Gamow          <- 1

## Constant scale of HF, or assume HF is correct at 10 GK and
## extrapolate to there?
extrapolateHF <- TRUE    ## TRUE="assume HF is correct at 10 GK

##sink("TMatch.log")

######################################################################
## You shouldn't need to touch anything below this line...
library("RColorBrewer")
cols <- brewer.pal(4,"Dark2")

## For warning messages
redtext <- "\033[0;31m"
normtext <- "\033[0m"

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

## Attempt to read the mass table. If unsuccessful, write a message
masses <- try(read.fortran(massfile,
			   skip=37,
			   format=c("4X","3I5","A3","84X","I3","I7","1X","I5"),
			   col.names=c("N","Z","A","El","Mass1","Mass2","Mass3")),
	      silent=TRUE)
if(!inherits(masses,"try-error")){
    ## If last digits don't exist, just make it zero
    masses$Mass3[is.na(masses$Mass3)] <- 0
    ## Calculate the total mass by combining the masses in the table
    masses$Mass <- masses$Mass1+masses$Mass2*1e-6+masses$Mass3*1e-11
    masses$El <- gsub(" ", "", masses$El)
} else {
    cat("\n",
	"WARNING: Mass file:",massfile,"not found\n")
}

spins <- try(read.fortran(nubasefile,
			  skip=25,
			  format=c("10X","A6","72X","A5"),
			  col.names=c("iso","Jpi")),
	     silent=TRUE)
if(!inherits(spins,"try-error")){
    spins$iso <- gsub(" ","",spins$iso)
    spins$Jpi <- gsub(" ","",spins$Jpi)
    spins$Jpi <- gsub("[#+-]","",spins$Jpi)
    spins$Jpi <- gsub("[*]","",spins$Jpi)
    spins$Jpi <- gsub("[(]","",spins$Jpi)
    spins$Jpi <- gsub("[)]","",spins$Jpi)
} else {
    cat("\n",
	"WARNING: NuBase file:",nubasefile,"not found\n")
}


## Read the masses and charges
getZ <- function(Read){
    A <- as.numeric(gsub("[^0-9.-]", "", Read))
    El <- gsub("[^A-Z,a-z]", "", Read)
    
    Z <- masses[masses$El==El,'Z'][1]
    Z
}
getM <- function(Read){
    A <- as.numeric(gsub("[^0-9.-]", "", Read))
    El <- gsub("[^A-Z,a-z]", "", Read)
    
    M <- masses[((masses$El==El) & (masses$A==A)),'A'][1]
    M
}
getJ <- function(Read){
    Jpi <- spins[spins$iso == Read,'Jpi']
    return <- NA
    if(length(grep("/",Jpi))>0){
	return <- sapply(strsplit(Jpi, split="/"),
	       function(x)as.numeric(x[1])/as.numeric(x[2]))
    } else {
	return <- as.numeric(Jpi)
    }
    return[1]
}

Read <- scan(RatesMCFile,what=character(),skip=2,n=1,quiet=TRUE)
Z0 <- suppressWarnings(as.double(Read))
if(is.na(Z0))Z0 <- getZ(Read)
Read <- scan(RatesMCFile,what=character(),skip=3,n=1,quiet=TRUE)
Z1 <- suppressWarnings(as.double(Read))
if(is.na(Z1))Z1 <- getZ(Read)
Read <- scan(RatesMCFile,what=character(),skip=5,n=1,quiet=TRUE)
M0 <- suppressWarnings(as.double(Read))
if(is.na(M0))M0 <- getM(Read)
Read <- scan(RatesMCFile,what=character(),skip=6,n=1,quiet=TRUE)
M1 <- suppressWarnings(as.double(Read))
if(is.na(M1))M1 <- getM(Read)
Read <- scan(RatesMCFile,what=character(),skip=8,n=1,quiet=TRUE)
J0 <- suppressWarnings(as.double(Read))
if(is.na(J0))J0 <- getJ(Read)
Read <- scan(RatesMCFile,what=character(),skip=9,n=1,quiet=TRUE)
J1 <- suppressWarnings(as.double(Read))
if(is.na(J1))J1 <- getJ(Read)

mu <- M0*M1/(M0+M1)

## Quit if masses of Jpi are NA
if(is.na(Z0) | is.na(Z1) | is.na(M0) | is.na(M1) | is.na(J0) | is.na(J1)){
    cat(paste0(redtext,"Unfortunately the AME/NuBase inputs could not be read properly.\nPlease enter the masses, spins, and J by hand into RatesMC.in\n",normtext))
    stop()
}
cat("\nIf you're using AME Mass or Nubase inputs,",
    "\nplease double-check the following:\n")
cat("Z0, M0, J0 = ",Z0,",",M0,",",J0,"\n")
cat("Z1, M1, J1 = ",Z1,",",M1,",",J1,"\n")


## The RatesMC rates
Rates <- read.table(RatesMCOutput,skip=4)
T <- Rates[,1]

## Find the start of the resonances section
lines <- readLines(RatesMCFile)
skip <- grep("Resonant Contribution",lines)+3

if(ReadCont){
    ## Read in the resonance energies.
    ## This procedure should read in all resonance energies without erroring!
    skip <- skip
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
    skip <- skip
    ## Keep looping until we're done
    while(TRUE){
	## Read a line
	x <- try(read.table(RatesMCFile,skip=skip,nrows=1))
	## Break the loop if we see *s
	if(substr(x[,1],1,2) == "**") break
	## Is the first entry a '+' or '!'? Ignore if so
	if(x[1] == "+" || x[1] == "!" || length(grep("!",x[1])) || sign(x[1]) == -1) {
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

    plot(c(0,max(Energies)),c(0,max(resCont[i,])),type='n',xlab="E",ylab="Cont",
	 main=paste(T[i]," (i=",i,")",sep=""))
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
    ##if(sum(cont>0.5) < 2){
    ##if(max(cont) < 0.8){   ## If a single resonance is contributing
			   ## most, the ETER algorithm doesn't work!

	interp.ETER <- approxfun(y=Energies.o,x=cumsum(cont/sum(cont)),ties="ordered")
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
cut <- T>0.1 & !is.na(ETER.array)
if(sum(cut, na.rm=TRUE)>0){
    interp.ETER <- approxfun(x=ETER.array[cut],
			     y=T[cut])
    TMatch.ETER <- interp.ETER(maxE)
    cat("TMatch from ETER method:       ",TMatch.ETER,"\n")

    ##if(ETER.array[length(ETER.array)] < maxE){
    if(!any(ETER.array[cut] > maxE)){
	cat(paste0(redtext,
		   "\nExperimental rate looks good up to 10 GK!\nCheck TMatch.pdf\n",
		   normtext))
	TMatch.ETER <- 11
##    dev.off()
##    stop()
    }
} else {
    TMatch.ETER <- NA
}

if(is.na(TMatch.ETER)){
##    dev.off()
    cat(paste0(redtext,
	       "\nSomething went wrong finding the ETER matching T!\n",
	       "There likely aren't enough resonances to use this method.\n",
	       "Check TMatch.pdf\n",
	       normtext))
    cat("Using Gamow window matching temperature\n\n")
    ##    stop()
    TMatch.ETER <- TMatch.Gamow
}

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

cat("The scaling factor to apply to HF: N =",Norm.med,"at",TMatch.ETER,"GK\n")

## Construct a normalization vector over the length of HF[T>TMatch]
if(extrapolateHF & TMatch.ETER<10){
    Norm.med.v <- rep(Norm.med,length(HF[,2]))
    Norm.low.v <- rep(Norm.low,length(HF[,2]))
    Norm.high.v <- rep(Norm.high,length(HF[,2]))
    cut <- HF[,1]>TMatch.ETER
    ## Using experimental uncertainty at TMatch
    ##Norm.med.v[cut] <- approx(x=c(TMatch.ETER,10),y=c(Norm.med,1),
    ##			      xout=HF[HF[,1]>TMatch.ETER,1])$y
    ##Norm.low.v[cut] <- approx(x=c(TMatch.ETER,10),y=c(Norm.low,Norm.low/Norm.med),
    ##			      xout=HF[HF[,1]>TMatch.ETER,1])$y
    ##Norm.high.v[cut] <- approx(x=c(TMatch.ETER,10),y=c(Norm.high,Norm.high/Norm.med),
    ##			       xout=HF[HF[,1]>TMatch.ETER,1])$y

    ## Assuming HF is correct at 10 GK with a factor of 10 uncertainty
    Norm.med.v[cut] <- approx(x=c(TMatch.ETER,10),y=c(Norm.med,1),
			      xout=HF[HF[,1]>TMatch.ETER,1])$y
    Norm.low.v[cut] <- Norm.med.v[cut]*approx(x=c(TMatch.ETER,10),
					       y=c(Norm.low/Norm.med,0.1),
			xout=HF[HF[,1]>TMatch.ETER,1])$y
    Norm.high.v[cut] <- Norm.med.v[cut]*approx(x=c(TMatch.ETER,10),
					       y=c(Norm.high/Norm.med,10),
			xout=HF[HF[,1]>TMatch.ETER,1])$y


    ## Now the matched HF rate!
    HF.matched.med  <- HF[,2]*Norm.med.v
    HF.matched.low  <- HF[,2]*Norm.low.v
    HF.matched.high <- HF[,2]*Norm.high.v
} else {
    HF.matched.med  <- HF[,2]*Norm.med
    HF.matched.low  <- HF[,2]*Norm.low
    HF.matched.high <- HF[,2]*Norm.high
}
HF <- cbind(HF, HF.matched.low, HF.matched.med, HF.matched.high)

## We need to interpolate the matched HF rate at the RatesMC temperatures
interp.HF.low  <- approxfun(HF[,1], HF[,3])
interp.HF.med  <- approxfun(HF[,1], HF[,4])
interp.HF.high <- approxfun(HF[,1], HF[,5])
HF.final.low  <- interp.HF.low(Rates[,1])
HF.final.med  <- interp.HF.med(Rates[,1])
HF.final.high <- interp.HF.high(Rates[,1])

## Plot the HF and, scaled HF, and interpolated HF to check it all went well!
layout(matrix(c(1:2), 2, 1, byrow=TRUE))
cut <- T > TMatch.ETER
cut[length(cut)-(2+sum(cut))] <- TRUE
cut.HF <- HF[,1] > TMatch.ETER
cut.HF[length(cut.HF)-(2+sum(cut.HF))] <- TRUE
##ylim <- range(c(Rates[,2],HF[,2]))
ylim <- range(c(Rates[cut,2],HF.final.med[cut],HF[cut.HF,2]))

plot(Rates[,1], Rates[,3], type='l', ylim=ylim,xlim=range(T[cut]),
     xlab="Temperature", ylab="Rate", col=cols[1], log="xy",
     xaxs='i',main=paste("TMatch = ",format(TMatch.ETER,digits=2,nsmall=2),"GK\n",
			 "Scale at TMatch = ",format(Norm.med,digits=2,nsmall=2)))
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

fu <- rowMeans(cbind(HF.final.med/HF.final.low,HF.final.high/HF.final.med))


Rates.final <- Rates[cut,]
Rates.final <- rbind(Rates.final,
		     cbind(Rates[acut,1],
			   HF.final.low[acut], HF.final.med[acut], HF.final.high[acut],
			   fu[acut]))

## Plot the final rate in a nice simple figure
cut <- T > TMatch.ETER
cut[length(cut)-(2+sum(cut))] <- TRUE
ylim <- range(Rates.final[cut,2:4])
plot(Rates[,1], Rates[,3], type='l', ylim=ylim,xlim=range(T[cut]),
     xlab="Temperature", ylab="Rate", col=cols[1], log="xy",
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
##OP.above <- cbind(OP.T, "  (", OP.Low, ")", "       (", OP.Med,")", "     (", OP.High, ")","      (", OP.fu,")")
OP.above <- cbind(OP.T, "   ", OP.Low, "       ", OP.Med, "       ", OP.High, " ","       ", OP.fu," ")
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
