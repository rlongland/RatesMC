## Plot the uncertainty bands
## RatesMC2
## Author: Richard L. Longland
## Date:   Dec. 2021
########################################

outputfile <- "../RatesMC.out"
litfilename <- ""#"RatesMC-comparison.out"
sampfile <- "../RatesMC.samp"

## Include another comparison rate?
## Include this as a separate [T, rate] file
ExtraRate <- FALSE
extrafilename <- "testing.dat" #"Tmatch.HF"

## Color scheme ("Blues" from ColorBrewer)
##pal <- c("#DEEBF7", "#9ECAE1", "#3182BD")
pal <- c("#BDD7E7", "#6BAED6", "#2171B5")

## Temeprature range to plot
TMin <- 0.01
TMax <- 10
## Y-axis range
YRangeUser <- NULL ##c(0.1e-10,1e10)

## Control to change the axis and tick label scales
axislabelscale <- 1.0
ticklabelscale <- 1.0

drawGrid <- TRUE

######################################################
# Function to print the displayed 2D curve to a file
lognorm <- function(x,mu,sigma){
  (1/(x*sigma*sqrt(2*3.142))) * exp(-((log(x)-mu)^2)/(2*sigma^2))
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
Samples <- as.double(read.table(sampfile,skip=3,header=FALSE,nrows=1)[3])
data <- read.table(outputfile,skip=4,header=FALSE,stringsAsFactors=FALSE)
## LitData is the data to compare with
## TODO
##LitBool <- FALSE
##if(file.exists(litfilename)){
##  LitData <- read.table(litfilename,header=FALSE,skip=headerskip,
##                        stringsAsFactors=FALSE)
##  LitBool <- TRUE
##}
##if(ExtraRate)
##  ExtraData <- read.table(extrafilename,skip=1,header=FALSE,
##                          stringsAsFactors=FALSE)

## If the rates file is Tmatch.out, we expect 7 columns, if it's
## RatesMC.out, there's 9
## TODO
## if(length(data[1,])==7){outisTMatch=TRUE}else{outisTMatch=FALSE}

##if(outisTMatch){
##  TMatch <- as.double(read.table(outputfile,skip=1,header=FALSE,nrows=1)[11])
##} else {
TMatch <- 11
##}
outisTMatch = FALSE

## TODO
## if(LitBool){
##   for(i in 1:length(LitData[,1])){
##     for(j in 2:4){
##       LitData[i,j] <- gsub("([()])","",LitData[i,j])
##     }
##   }
## 
##   LitData$V1 <- as.double(LitData$V1)
##   LitData$V2 <- as.double(LitData$V2)
##   LitData$V3 <- as.double(LitData$V4)
##   LitData$V4 <- as.double(LitData$V6)
## }

# Temp is data$V1
# Low rate is data$V2
# Median rate is data$V3
# High rate is data$V4
if(outisTMatch){
  cat("The rate data file is in Tmatch format\n")
  MyRate <- cbind(as.double(data$V1),as.double(data$V2),as.double(data$V3),
                  as.double(data$V4),as.double(data$V5),as.double(data$V6))
} else {
  cat("The rate data file is in RatesMC format\n")
  MyRate <- cbind(as.double(data$V1),as.double(data$V2),as.double(data$V3),
                  as.double(data$V4))
}
# The MC uncertainty ratios
CompRate <- cbind(data$V1,MyRate[,2]/MyRate[,3],MyRate[,4]/MyRate[,3])

## sometimes, the literature has no uncertainty bands. Check for this:
## TODO
## if(LitBool){
##   if(length(LitData[1,]) > 2){
##     ## The Literature uncertainty bands
##     LitCompRate <- cbind(LitData$V1,LitData$V2/LitData$V3,LitData$V4/LitData$V3)
##   }
## 
##   ## First the temperature grids have to be matched
##   matchlist1 <- as.double(na.omit(match(LitData[,1],MyRate[,1])))
##   matchlist2 <- as.double(na.omit(match(MyRate[,1],LitData[,1])))
## 
##   RateComp <- cbind(MyRate[matchlist1,1],
##                     LitData[matchlist2,2]/MyRate[matchlist1,3],
##                     LitData[matchlist2,3]/MyRate[matchlist1,3],
##                     LitData[matchlist2,4]/MyRate[matchlist1,3])
## }
## if(ExtraRate){
##   matchliste1 <- as.double(na.omit(match(ExtraData[,1],MyRate[,1])))
##   matchliste2 <- as.double(na.omit(match(MyRate[,1],ExtraData[,1])))
## 
##   ExtraComp <- cbind(MyRate[matchliste1,1],
##                      ExtraData[matchliste2,2]/MyRate[matchliste1,3])
## }


## -----------------------------------------------------------------
## Read in the rate samples and construct probability density
## functions out of them
ntemps <- dim(data)[1]
csum <- array(,dim=c(Samples,2,ntemps))
onesig <- array(,dim=c(ntemps,2))
twosig <- array(,dim=c(ntemps,2))
threesig <- array(,dim=c(ntemps,2))

pb <- txtProgressBar(min=1,max=length((1:ntemps)[MyRate[,1]<TMatch & MyRate[,3]>0]),style=3)
cat("   Reading samples...\n")

for(i in (1:ntemps)[MyRate[,1]<TMatch & MyRate[,3]>0]){
  setTxtProgressBar(pb, i)

  tmp <- withRestarts(scan(sampfile,what=double(),nmax=Samples,
                               skip=1+(i*3+(i-1)*Samples),quiet=TRUE),
                          abort=function(){ })
  tmp <- sort(tmp)
  tmp[tmp==0] <- NA
  
  csum[,,i] <- cbind(tmp/MyRate[i,3],seq(from=0,to=2,length.out=Samples))
  
  ## Find the 1- and 2-sigma uncertainties
  onesig[i,] <- quantile(csum[,1,i],probs=c(0.16,0.84),na.rm=TRUE)
  twosig[i,] <- quantile(csum[,1,i],probs=c(0.05,0.95),na.rm=TRUE)
  threesig[i,] <- quantile(csum[,1,i],probs=c(0.01,0.99),na.rm=TRUE)
  
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
if(!is.null(YRangeUser)){
    YRange <- YRangeUser
} else {
    YRange <- range(threesig,na.rm=TRUE)
}

## The first pdf of the uncertainty bands
mypdf(file="GraphUncertainties.pdf",width=6,height=6,onefile=F)
##while(dev.cur()>1)dev.off()
##myX11(width=8,height=8)

## set up the parameters for plotting
##oldsettings <- par(cex.lab=2.2,cex.axis=1.8,tcl=0.8)

## Finally make the plot.
plot(1,1,type='n',
     xlim=c(TMin,TMax),
     ylim=YRange,
     xaxs='i',
     log='xy',
     yaxt='n',xaxt='n',
     xlab="Temperature (GK)", ylab="Reaction Rate Ratio",
     cex.lab=par("cex.lab")*axislabelscale)

## x-axis ticks
aX <- c(-2,-1,0,1,2)

## y-axis ticks
minBase10 <- floor(log10(YRange[1]))
maxBase10 <- ceiling(log10(YRange[2]))
aY <- pretty(seq(minBase10,maxBase10,length.out=maxBase10-minBase10))


Temps <- MyRate[,1]
## 3-sigma uncertainties
xx <- c(Temps,rev(Temps),Temps[1])
yy <- c(threesig[,1],rev(threesig[,2]),threesig[1,1])
polygon(xx,yy,col=pal[1],border=NA)

## 2-sigma uncertainties
xx <- c(Temps,rev(Temps),Temps[1])
yy <- c(twosig[,1],rev(twosig[,2]),twosig[1,1])
polygon(xx,yy,col=pal[2],border=NA)

## 1-sigma uncertainties
xx <- c(Temps,rev(Temps),Temps[1])
yy <- c(onesig[,1],rev(onesig[,2]),onesig[1,1])
polygon(xx,yy,col=pal[3],border=NA)


if(drawGrid){
    abline(v=10^aX,lty="dotted",col="grey",lwd=0.5)
    abline(h=10^aY[aY!=0],lty="dotted",col="grey",lwd=0.5)
}

## Now add the X-axis ticks
axis(1,at=10^aX,label=10^aX,cex.axis=par("cex.axis")*ticklabelscale,
     padj=0.3)
axis(3,at=aX,labels=FALSE)
majorx <- 10^aX
minorx<-array(0,dim=10*(length(majorx)-1))

## now populate the minor tick marks
for(i in 1:(length(majorx-1))){
    minorx[(1+(i-1)*10):(i*10)] <- seq(majorx[i],majorx[i]*10,majorx[i])
}
axis(1,at=minorx,labels=FALSE,tcl=0.4)
axis(3,at=minorx,labels=FALSE,tcl=0.4)


## and the Y-axis

axis(2,at=10^aY,label=axTexpr(2,10^aY),cex.axis=par("cex.axis")*ticklabelscale)
axis(4,at=10^aY,labels=FALSE)

majory <- aY
majory <- c(10^majory[1],10^majory,10^majory[length(majory)]*10)
minory<-array(0,dim=10*(length(majory)-1))
## for y-ticks, extend further
for(i in 1:(length(majory-1))){
    minory[(1+(i-1)*10):(i*10)] <- seq(majory[i],majory[i]*10,majory[i])
}
axis(2,at=minory,labels=FALSE,tcl=0.4)
axis(4,at=minory,labels=FALSE,tcl=0.4)



abline(h=1,lty=2,lwd=2)

box()

## The reaction name
xx <- grconvertX(0.88,from="nfc",to="user")
yy <- grconvertY(0.88,from="nfc",to="user")
text(parse(text=ReacName),x=xx,y=yy,cex=1.7,adj=c(1,1))

dev.off()

##TODO
## if(LitBool){
##   write.table(RateComp,"GraphCompare.dat",col.names=FALSE,row.names=FALSE,sep="\t")
## }


