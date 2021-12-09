# Plot the uncertainty bands
# Author: Richard L. Longland
########################################

outputfile <- "RatesMC.out"
litfilename <- ""#"RatesMC-comparison.out"
histfile <- "RatesMC.hist"
sampfile <- "RatesMC.samp"
headerskip <- 3  ## 1 for Normal literature
                 ## or 3 for RatesMC type files

## Include another comparison rate?
## Include this as a separate [T, rate] file
ExtraRate <- FALSE
extrafilename <- "testing.dat" #"Tmatch.HF"

## Temeprature range to plot
TMin <- 0.01
TMax <- 10
## Y-axis range
#YRangeUser <- c(0.1,10)

## Control to change the axis and tick label scales
axislabelscale <- 1
ticklabelscale <- 1

## Number of colours
nColours <- 100

## Linear color scale?
LinScale <- FALSE

## Type of figure output ("pdf" or "png")
outputtype <- "png"

######################################################
# Function to print the displayed 2D curve to a file
lognorm <- function(x,mu,sigma){
  (1/(x*sigma*sqrt(2*3.142))) * exp(-((log(x)-mu)^2)/(2*sigma^2))
}
mypng <- function(file="output.png",...)
  {
    grDevices::png(file=file,units="in",res=150,...)
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
    grDevices::cairo_pdf(file=file,...)
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
  eT <- floor(log10(abs(at)))# at == 0 case is dealt with below
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

library("gsl")

## Read and format the reaction name
ReacName = as.character(read.table(outputfile,skip=0,header=FALSE,nrows=1)$V1)
# Make superscripts, subscripts etc
ReacName <- sub(",g\\)",",gamma\\)",ReacName)
ReacName <- sub("\\(g,","\\(gamma,",ReacName)
ReacName <- sub(",a\\)",",alpha\\)",ReacName)
ReacName <- sub("\\(a,","\\(alpha,",ReacName)
ReacName <- sub("\\)","\\)*",ReacName)
ReacName <- gsub("([[:digit:]]+)", "phantom()^{\\1}*", ReacName)


## ------------------------------------------------------
## Read in the rates and compare them with the literature

## Read in my data
Samples <- as.double(read.table(sampfile,skip=2,header=FALSE,nrows=1)[3])
data <- read.table(outputfile,skip=3,header=FALSE,stringsAsFactors=FALSE)
# LitData is the data to compare with
LitBool <- FALSE
if(file.exists(litfilename)){
  LitData <- read.table(litfilename,header=FALSE,skip=headerskip,
                        stringsAsFactors=FALSE)
  LitBool <- TRUE
}
if(ExtraRate)
  ExtraData <- read.table(extrafilename,skip=1,header=FALSE,
                          stringsAsFactors=FALSE)

# If the rates file is Tmatch.out, we expect 7 columns, if it's
# RatesMC.out, there's 9
if(length(data[1,])==7){outisTMatch=TRUE}else{outisTMatch=FALSE}

if(outisTMatch){
  TMatch <- as.double(read.table(outputfile,skip=1,header=FALSE,nrows=1)[11])
} else {
  TMatch=10
}


# get rid of brackets in Tmatch type files
for(i in 1:length(data[,1])){
  for(j in 2:6){
    data[i,j] <- gsub("([()])","",data[i,j])
  }
}
if(LitBool){
  for(i in 1:length(LitData[,1])){
    for(j in 2:4){
      LitData[i,j] <- gsub("([()])","",LitData[i,j])
    }
  }

  LitData$V1 <- as.double(LitData$V1)
  LitData$V2 <- as.double(LitData$V2)
  LitData$V3 <- as.double(LitData$V4)
  LitData$V4 <- as.double(LitData$V6)
}

# Temp is data$V1
# Low rate is data$V2
# Median rate is data$V4
# High rate is data$V6
if(outisTMatch){
  cat("The rate data file is in Tmatch format\n")
  MyRate <- cbind(as.double(data$V1),as.double(data$V2),as.double(data$V3),
                  as.double(data$V4),as.double(data$V5),as.double(data$V6))
} else {
  cat("The rate data file is in RatesMC format\n")
  MyRate <- cbind(as.double(data$V1),as.double(data$V2),as.double(data$V4),
                  as.double(data$V6),as.double(data$V7),as.double(data$V8))
}
# The MC uncertainty ratios
CompRate <- cbind(data$V1,MyRate[,2]/MyRate[,3],MyRate[,4]/MyRate[,3])

# sometimes, the literature has no uncertainty bands. Check for this:
if(LitBool){
  if(length(LitData[1,]) > 2){
    ## The Literature uncertainty bands
    LitCompRate <- cbind(LitData$V1,LitData$V2/LitData$V3,LitData$V4/LitData$V3)
  }

  ## First the temperature grids have to be matched
  matchlist1 <- as.double(na.omit(match(LitData[,1],MyRate[,1])))
  matchlist2 <- as.double(na.omit(match(MyRate[,1],LitData[,1])))

  RateComp <- cbind(MyRate[matchlist1,1],
                    LitData[matchlist2,2]/MyRate[matchlist1,3],
                    LitData[matchlist2,3]/MyRate[matchlist1,3],
                    LitData[matchlist2,4]/MyRate[matchlist1,3])
}
if(ExtraRate){
  matchliste1 <- as.double(na.omit(match(ExtraData[,1],MyRate[,1])))
  matchliste2 <- as.double(na.omit(match(MyRate[,1],ExtraData[,1])))

  ExtraComp <- cbind(MyRate[matchliste1,1],
                     ExtraData[matchliste2,2]/MyRate[matchliste1,3])
}


## -----------------------------------------------------------------
## Read in the rate histograms and construct probability density
## functions out of them
ntemps <- 51
##hist <- array(,dim=c(1000,2,ntemps))
csum <- array(,dim=c(Samples,2,ntemps))
##csumi <- array(,dim=c(1000,2,ntemps))
twosig <- array(,dim=c(ntemps,2))

pb <- txtProgressBar(min=1,max=length((1:ntemps)[MyRate[,1]<TMatch & MyRate[,3]>0]),style=3)
cat("   Reading samples...\n")

for(i in (1:ntemps)[MyRate[,1]<TMatch & MyRate[,3]>0]){
  setTxtProgressBar(pb, i)

  tmp <- withRestarts(scan(sampfile,what=double(),nmax=Samples,
                               skip=i*3+(i-1)*Samples,quiet=TRUE),
                          abort=function(){ })
  tmp <- sort(tmp)
  ## Reshape, and then enter the info into the hist array
##  tmp <- array(tmp,dim=c(1000,3))
##  hist[,,i] <- tmp[,c(1,3)]
##  hist[,2,i]<-(hist[,2,i])/sum(hist[,2,i])  # Normalise
##  ## Make the cumulative distribution (and scale to the recommended
##  ## rate)
##  csum[,,i] <- cbind(hist[,1,i]/MyRate[i,3],2*cumsum(hist[,2,i]))

  csum[,,i] <- cbind(tmp/MyRate[i,3],seq(from=0,to=2,length.out=Samples))
  
  ## Find the 2-sigma uncertainties
  twosig[i,] <- quantile(csum[,1,i],probs=c(0.05,0.95))
  ##  s <- 1:length(csum[,2,i])
##  sfun <- splinefun(csum[s,2,i],log10(csum[s,1,i]))
##  twosig[i,] <- 10^sfun(c(0.05,1.95))
###  twosig[i,twosig[i,]<0] <- 0
##  if(hist[1,2,i] > 0.9)twosig[i,] <- c(NA,NA)
  
  ## Adjust the cumulative distribution to correspond to percentiles
  csum[csum[,2,i]>1,2,i] <- 2-csum[csum[,2,i]>1,2,i]
  csum[,2,i] <- 1-csum[,2,i]
}
close(pb)
## -----------------------------------------------------------------
## Make the plots

## Calculate the full range to plot. Use the minimum and maximum
## ratios from the histograms since the range is already calculated for
## them by the RatesMC code
if(exists("YRangeUser")){
  FRange <- YRangeUser
}else{
  FRange <- range(csum[,1,],na.rm=TRUE)
}
## Now make a vector (log10 spaced) to construct the lognormal matrix
ngrid <- 100
Ratios <- 10^seq(from=log10(FRange[1]),to=log10(FRange[2]),length.out=ngrid)

## Compute the percentile matrix. The zero percentile corresponds to
## 50-50, 10th percentile to 45-55, 68th percentile to 16-84 etc.

ntemps <- 51 #length(MyRate[MyRate[,1]>=TMin,1])
itemps <- (1:ntemps)[MyRate[,1]>TMin & MyRate[,1]<TMax]
ntemps <- length(itemps)

pb <- txtProgressBar(min=1,max=length((1:ntemps)[MyRate[,1]<TMatch &
                             MyRate[,3]>0]),style=3)
cat("   Calculating contours...\n")


lnmat <- array(,dim=c(ntemps,ngrid))
for(i in itemps[MyRate[itemps,1]<TMatch & MyRate[itemps,3]>0]){ #(1:ntemps)[MyRate[,1]<TMatch & MyRate[,1]>TMin]){

  setTxtProgressBar(pb, i)
  
  ## Set up the interpolator (rule=2 means use the closed real number)
  interpfun <- approxfun(csum[,1,i],csum[,2,i],rule=2)
  lnmat[1+i-itemps[1],] <- interpfun(Ratios)

}
close(pb)

for(i in itemps[MyRate[itemps,1]>=TMatch]){  #(1:ntemps)[MyRate[,1]>=TMatch]){
  mu <- MyRate[i,5]
  sigma <- MyRate[i,6]
  scale <- lognorm(exp(mu),mu,sigma)  ## To make max=1

  lnmat[1+i-itemps[1],] <- abs(1-erfc(-(log(Ratios*MyRate[i,3]) - mu)/
                                      (sigma*sqrt(2))))
  twosig[i,] <- qlnorm(c(0.025,0.975),MyRate[i,5], MyRate[i,6])/
              MyRate[i,3]
##  lnmat[i,] <- 1-lnmat[i,]  ## Peak=1
}

## The first pdf if the uncertainty contour
if(outputtype == "pdf")
  mypdf(file="GraphContour.pdf",width=10,height=10,onefile=F)
if(outputtype == "png")
  mypng(file="GraphContour.png",width=10,height=10)

# set up the parameters for plotting
oldsettings <- par(cex.lab=2.2,cex.axis=1.8,tcl=0.8)

## Number of colours to use
len <- nColours

lnlim <- range(lnmat,finite=TRUE)
if(!LinScale)lnlim <- 10^lnlim
levels <- pretty(lnlim, len) ## Make the colours
##palette <- colorRampPalette(c("white","darkgoldenrod3","red3","black"))
#library("RColorBrewer")
#palette <- colorRampPalette(c("white",brewer.pal(9,"YlOrRd")),bias=1.5)
palette <- colorRampPalette(c("white","#FFFFCC","#FFEDA0","#FED976",
                              "#FEB24C","#FD8D3C","#FC4E2A","#E31A1C",
                              "#BD0026","#800026"),bias=1.3)
pal <- c(palette(1),palette(len))
#pal <- palette(len+1)
col <- pal[round(1+len-len*(levels-levels[1])/(lnlim[2]-lnlim[1]))]

## Finally make the contour plot. All extra lines must be included
## here because filled.contour will forget where it is, otherwise
addextra <- function(){
#  plot(range(log10(MyRate[,1])),range(y=log10(Ratios)),type="n",axes=FALSE)
  ## The line through unity
  abline(a=log10(1),b=0,lty=3,lwd=2)

  ## The current 1- and 2-sigma uncertainty bands
  lines(log10(CompRate[,1]),log10(CompRate[,2]),lwd=3)
  lines(log10(CompRate[,1]),log10(CompRate[,3]),lwd=3)
  twosig[twosig<FRange[1]] <- FRange[1]
  lines(log10(CompRate[,1]),log10(twosig[,1]),lwd=1)
  lines(log10(CompRate[,1]),log10(twosig[,2]),lwd=1)

  if(LitBool){
    ## Then draw the literature uncertainties
    lines(log10(RateComp[,1]),log10(RateComp[,2]),
          lwd=2,lty=2,col="blue")
    lines(log10(RateComp[,1]),log10(RateComp[,3]),
          lwd=3,lty=1,col="blue")
    lines(log10(RateComp[,1]),log10(RateComp[,4]),
          lwd=2,lty=2,col="blue")
  }
  ## Draw the extra rate if needed
  if(ExtraRate){
    lines(log10(ExtraComp[,1]),log10(ExtraComp[,2]),
        lwd=3,lty=1,col="forestgreen")  ## forestgreen green3
  }
  ## Now add the X-axis ticks
  aX <- c(-2,-1,0,1,2)
  axis(1,at=aX,label=10^aX,cex.axis=par("cex.axis")*ticklabelscale,
       padj=0.3)
  axis(3,at=aX,labels=FALSE)
  majorx <- 10^aX
  minorx<-array(0,dim=10*(length(majorx)-1))

  ## now populate the minor tick marks
  for(i in 1:(length(majorx-1))){
    minorx[(1+(i-1)*10):(i*10)] <- seq(majorx[i],majorx[i]*10,majorx[i])
  }
  axis(1,at=log10(minorx),labels=FALSE,tcl=0.4)
  axis(3,at=log10(minorx),labels=FALSE,tcl=0.4)

  # and the Y-axis
  minBase10 <- floor(log10(FRange[1]))
  maxBase10 <- ceiling(log10(FRange[2]))
  aY <- axTicks(2,axp=c(minBase10,maxBase10,maxBase10-minBase10))
  axis(2,at=aY,label=axTexpr(2,10^aY),cex.axis=par("cex.axis")*ticklabelscale)
  axis(4,at=aY,labels=FALSE)

  majory <- aY
  majory <- c(10^majory[1],10^majory,10^majory[length(majory)]*10)
  minory<-array(0,dim=10*(length(majory)-1))
  ## for y-ticks, extend further
  for(i in 1:(length(majory-1))){
    minory[(1+(i-1)*10):(i*10)] <- seq(majory[i],majory[i]*10,majory[i])
  }
  axis(2,at=log10(minory),labels=FALSE,tcl=0.4)
  axis(4,at=log10(minory),labels=FALSE,tcl=0.4)

  ## Add the burning window
##  rect(xleft=log10(0.15),xright=log10(0.30),
##       ybottom=log10(1e-3),ytop=log10(3e-3),
##     col="tomato1",border="tomato1")
##  rect(xleft=log10(0.98),xright=log10(1.12),
##       ybottom=log10(1e-3),ytop=log10(3e-3),
##     col="tomato1",border="tomato1")
##  text(log10(0.215),log10(5e-3),"He",cex=2)
##  text(log10(1.05),log10(5e-3),"C",cex=2)

  ## The reaction name
  xx <- grconvertX(0.9,from="nfc",to="user")
  yy <- grconvertY(0.9,from="nfc",to="user")
  text(parse(text=ReacName),x=xx,y=yy,cex=2.5,adj=c(1,1))

  box()
  
}

plotmat <- 10^lnmat
if(LinScale)plotmat <- lnmat

filled.contour(x=log10(MyRate[itemps,1]),y=log10(Ratios),plotmat, 
               nlevels=len,levels=levels,col=col,
               xlab="Temperature (GK)",ylab="Reaction Rate Ratio",axes=FALSE,
               frame.plot = TRUE,
               plot.axes=addextra(),
               cex.lab=par("cex.lab")*axislabelscale,
               ## The key (make sure to correct for the sqrt)
               key.axes = {
                 rect(0, levels[-length(levels)], 1, levels[-1L],
                      col = col, border = NA)
                 aY <- seq(0, 100, len = 5)#c(0.05,0.2,0.5,0.8,0.95)
                 if(LinScale){
                   axis(4, at=aY/100,label=paste(aY,"%",sep=""))
                 }else{
                   axis(4, at=10^{aY/100},label=paste(aY,"%",sep=""))
                 }})

dev.off()


####----------------------
#### Now we can do the usual PlotCompare Plots
##pdf(file="PlotCompare.pdf",width=8.25,height=10.75,onefile=F)
####X11(width=8.25,height=10.75)
##layout(matrix(c(0,0,0,1,0,2,0,0), 4, 2, byrow = TRUE),
##  c(0.05,0.9),c(0.05,0.45,0.45,0.01))
##layout.show(2)
##
#### set up labels etc
##xlabel <- "Temperature (GK)"
##ylabel <- "Reaction rate ratio"
##
##
### set up the parameters for plotting
##oldsettings <- par(cex.lab=1.55,cex.main=1.55,cex.axis=1.7,yaxs="i",xaxs="i",lwd=1.5,
##                   tcl=-0.7,
##                   #mar=c(5.1,5.1,3.1,1.5))
##                   mar=c(4.7,3.1,1.1,2.0))
##
### find the log parts of the yaxis
###YAxisLims <- c(min(CompRate[,2:3],LitCompRate[,2:3]),
###               max(CompRate[,2:3],LitCompRate[,2:3]))
### include literature uncertainty band if present
##if(length(LitData[1,]) > 2){
##  YAxisLims <- c(min(CompRate[,2:3],LitCompRate[,2:3]),
##               max(CompRate[,2:3],LitCompRate[,2:3]))
##} else {
##  YAxisLims <- c(min(CompRate[,2:3]), max(CompRate[,2:3]))
##}                
##minBase10 <- floor(log10(YAxisLims[1]))
##maxBase10 <- ceiling(log10(YAxisLims[2]))
##
### plot the lower rate
##plot(CompRate[,1],CompRate[,2],ylim=c(YAxisLims[1]*0.96,YAxisLims[2]*1.04),
##     xlim=c(0.01,TMax+.001),
##     type="l",log="xy",xlab="",ylab="",
##     lwd=2,xaxp=c(0.01,TMax,1),yaxp=c(10^minBase10,10^maxBase10,1))
### add line for upper rate
##lines(CompRate[,1],CompRate[,3],lwd=2)
### Add lines for the literature uncertainty bands if present
##if(length(LitData[1,]) > 2){
##  lines(LitCompRate[,1],LitCompRate[,2],lwd=2,lty=2)
##  lines(LitCompRate[,1],LitCompRate[,3],lwd=2,lty=2)
##}
### add line through 1
##lines(c(0.01,10),c(1,1),lty=3,lwd=2)
##
##
### Now add the minor tick marks
##majorx <- axTicks(1,axp=c(0.01,10,1))
##majory <- axTicks(2,axp=c(10^minBase10,10^maxBase10,1))
##majory <- c(majory[1],majory,majory[length(majory)]*10)
##
##minorx<-array(0,dim=10*(length(majorx)-1))
##minory<-array(0,dim=10*(length(majory)-1))
### now populate the minor tick marks
##for(i in 1:(length(majorx-1))){
##  minorx[(1+(i-1)*10):(i*10)] <- seq(majorx[i],majorx[i]*10,majorx[i])
##}
##axis(1,at=minorx,labels=FALSE,tcl=0.4)
### for y-ticks, extend further
##for(i in 1:(length(majory-1))){
##  minory[(1+(i-1)*10):(i*10)] <- seq(majory[i],majory[i]*10,majory[i])
##}
##axis(2,at=minory,labels=FALSE,tcl=0.4)
##
###output2D("Uncerts.ps")
##
##
###############
### Now do the second plot, comparisons
##
### First the temperature grids have to be matched
##matchlist1 <- as.double(na.omit(match(LitData[,1],MyRate[,1])))
##matchlist2 <- as.double(na.omit(match(MyRate[,1],LitData[,1])))
##
### if no uncertainties are on the literature, then column 2 is the recommended
##if(length(LitData[1,]) > 2){
##  reccol <- 3
##} else {
##  reccol <- 2
##}
##RateComp <- cbind(MyRate[matchlist1,1],
##                  MyRate[matchlist1,2]/LitData[matchlist2,reccol],
##                  MyRate[matchlist1,3]/LitData[matchlist2,reccol],
##                  MyRate[matchlist1,4]/LitData[matchlist2,reccol],
##                  MyRate[matchlist1,5],MyRate[matchlist1,6])
##
##minBase10 <- floor(log10(min(RateComp[,2:4])))
##maxBase10 <- ceiling(log10(max(RateComp[,2:4])))
##
### plot the lower rate
##plot(RateComp[,1],RateComp[,2],
##     ylim=c(min(RateComp[,2:4])*0.96,max(RateComp[,2:4])*1.04),
##     xlim=c(0.01,TMax+.001),
##       #xlim=c(min(RateComp[,1]),xMax+0.001),#max(RateComp[,1])+.001),
##     type="l",log="xy",#main=parse(text=ReacName),xlab=xlabel,ylab=ylabel,
##     xlab="",ylab="",
##     lwd=1,xaxp=c(0.01,TMax,1),yaxp=c(10^minBase10,10^maxBase10,1))
##  # add line for upper rate
##lines(RateComp[,1],RateComp[,4],lwd=1)
##  # Line for median
##lines(RateComp[,1],RateComp[,3],lwd=2)
##  # add line through 1
##lines(c(0.01,TMax),c(1,1),lty=3,lwd=2)
##
##  # Now add the minor tick marks
##majorx <- axTicks(1,axp=c(0.01,10,1))
##majory <- axTicks(2,axp=c(10^minBase10,10^maxBase10,1))
##majory <- c(majory[1]/10,majory,majory[length(majory)]*10)
##
##minorx<-array(0,dim=10*(length(majorx)-1))
##minory<-array(0,dim=10*(length(majory)-1))
##  # now populate the minor tick marks
##for(i in 1:(length(majorx-1))){
##  minorx[(1+(i-1)*10):(i*10)] <- seq(majorx[i],majorx[i]*10,majorx[i])
##}
##axis(1,at=minorx,labels=FALSE,tcl=0.4)
##  # for y-ticks, extend further
##for(i in 1:(length(majory-1))){
##  minory[(1+(i-1)*10):(i*10)] <- seq(majory[i],majory[i]*10,majory[i])
##}
##axis(2,at=minory,labels=FALSE,tcl=0.4)
##
##
##par(xpd=NA) 
##Minx <- par("usr")[1]
##Maxx <- par("usr")[2]
##Miny <- par("usr")[3]
##Maxy <- par("usr")[4]
##
##xpos <- 10^(Minx+0.45*(Maxx-Minx))
##ypos <- 10^(Miny-0.12*(Maxy-Miny))
##text(xpos,ypos,"Temperature (GK)",cex=2)
##xpos <- 10^(Minx-0.082*(Maxx-Minx))
##ypos <- 10^(Miny+0.5*(Maxy-Miny))
##text(xpos,ypos,"Reaction Rate Ratio",cex=2,srt=90)
##xpos <- 10^(Minx-0.082*(Maxx-Minx))
##ypos <- 10^(Maxy+0.68*(Maxy-Miny))
##text(xpos,ypos,"Reaction Rate Ratio",cex=2,srt=90)
##
##xpos <- 10^(Minx+0.45*(Maxx-Minx))
##ypos <- 10^(Maxy+1.26*(Maxy-Miny))
##text(xpos,ypos,parse(text=ReacName),cex=2)
##
##
##


#output2D("Comparison.ps")
##dev.off()

if(LitBool){
  write.table(RateComp,"GraphCompare.dat",col.names=FALSE,row.names=FALSE,sep="\t")
}


