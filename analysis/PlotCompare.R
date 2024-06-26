# Plot the uncertainty bands
# Author: Richard L. Longland
########################################

outputfile <- "RatesMC.out"
fullfile    <- "RatesMC.full"
litfilename <- "RatesMC-comparison.out"
headerskip <- 4  ## 3 for RatesMC type files
                 ## 4 for RatesMC2 type files
headerskip.lit <- 4  ## 3 for RatesMC type files
                     ## 4 for RatesMC2 type files

extraSuper <- ""            ## Use to put a superscript at the end
                            ## e.g. "g" for ^{26}Al^g

## Y-axis range (set to NULL for automatic plotting)
YRangeUser <- c(1e-2,1e2)

RateColor <- "gray"  # colour for the error region shade
RateTrans <- 0.8     # transparancy of error region (1 is fully opaque)
CompColor <- "blue"  # colour for the comparison shade
CompTrans <- 0.3     # transparancy of comparison

## Show New/Old rather than Old/New
NewOverLiterature <- FALSE

######################################################
# Function to print the displayed 2D curve to a file
lognorm <- function(x,mu,sigma){
  (1/(x*sigma*sqrt(2*3.142))) * exp(-((log(x)-mu)^2)/(2*sigma^2))
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
        mgp=c(3,0.5,0),
	family="sans")     
  }
## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}
# Complementary error function
erfc    <- function(x) {  # 1 - erf(x)
    2 * pnorm(-sqrt(2) * x)
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


##library("gsl")

#cat("\nWhat Temperature would you like to stop (default = 10)\n")
#TMax <- scan(nmax=1,quiet=TRUE)
#if(length(TMax)==0)TMax=10
TMax <- 10


## Read and format the reaction name
ReacName = as.character(read.table(outputfile,skip=0,header=FALSE,nrows=1)$V1)
# Make superscripts, subscripts etc
ReacName <- sub(",g\\)",",gamma\\)",ReacName)
ReacName <- sub("\\(g,","\\(gamma,",ReacName)
ReacName <- sub(",a\\)",",alpha\\)",ReacName)
ReacName <- sub("\\(a,","\\(alpha,",ReacName)
ReacName <- sub(",","*','*",ReacName)
##ReacName <- sub("\\(a*","\\(alpha",ReacName)
ReacName <- sub("\\)","\\)*",ReacName)
ReacName <- gsub("([[:digit:]]+)", "phantom()^{\\1}*", ReacName)
## Add a superscript
if(nchar(extraSuper)>0)ReacName <- paste(ReacName,"^",extraSuper,sep="")


data <- read.table(outputfile,skip=headerskip,header=FALSE,stringsAsFactors=FALSE)
fulldata <- read.table(fullfile,skip=headerskip,header=FALSE,stringsAsFactors=FALSE)
# CompData is the data to compare with
LitData <- read.table(litfilename,header=FALSE,skip=headerskip.lit,
                      stringsAsFactors=FALSE)

# If the file is Tmatch.out, we expect 7 columns, if it's RatesMC.out,
# there's 9
if(length(data[1,])==7){outisTMatch=TRUE}else{outisTMatch=FALSE}
#print(outisTMatch)

# get rid of brackets in Tmatch file
if(length(LitData[1,])>5){
    for(i in 1:length(data[,1])){
	for(j in 2:4){
	    data[i,j] <- gsub("([()])","",data[i,j])
	}
    }
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
  cat("The data file is in Tmatch format\n")
  MyRate <- cbind(as.double(data$V1),as.double(data$V2),as.double(data$V3),
                  as.double(data$V4),as.double(data$V5),as.double(data$V6))
} else {
  cat("The data file is in RatesMC format\n")
  MyRate <- cbind(as.double(data$V1),as.double(data$V2),as.double(data$V3),
                  as.double(data$V4))#,as.double(data$V7),as.double(data$V8))
}
					# The MC uncertainty ratios
cut <- MyRate[,3]>0
MyRate <- MyRate[cut,]
fulldata <- fulldata[cut,]
CompRate <- cbind(MyRate[,1],MyRate[,2]/MyRate[,3],MyRate[,4]/MyRate[,3])

cut <- LitData[,3]>0
LitData <- LitData[cut,]
# sometimes, the literature has no uncertainty bands. Check for this:
if(length(LitData[1,]) > 2){
  # The Literature uncertainty bands
  LitCompRate <- cbind(LitData$V1,LitData$V2/LitData$V3,LitData$V4/LitData$V3)
}

## Make a contour plot separately
minper <- 0.01   ## The percentages to use
maxper <- 0.99
## Calculate the full range to plot (lqnorm calculates the inverse
## lognormal error function, ie. the position corresponding to the
## desired percentile)
cut <- MyRate[,3]>0
fulldata <- fulldata[cut,]
MyRate <- MyRate[cut,]
CompRate <- CompRate[cut,]
FRange <- c(min(qlnorm(minper,fulldata[,9],fulldata[,10])/MyRate[,3]),
            max(qlnorm(maxper,fulldata[,9],fulldata[,10])/MyRate[,3]))
## Now make a vector (log10 spaced) to construct the lognormal matrix
Ratios <- 10^seq(from=log10(FRange[1]),to=log10(FRange[2]),length.out=100)

## Compute the lognormal matrix. This is the computed lognormal
## distribution at each temperature.
lnmat <- array(,dim=c(length(MyRate[,1]),100))
for(i in 1:length(MyRate[,1])){
  mu <- fulldata[i,9]
  sigma <- fulldata[i,10]
  scale <- lognorm(exp(mu),mu,sigma)  ## To make max=1

  lnmat[i,] <- abs(1-erfc(-(log(Ratios*MyRate[i,3]) - mu)/(sigma*sqrt(2))))
  lnmat[i,] <- 1-lnmat[i,]  ## Peak=1
}


##----------------------
## Now we can do the usual PlotCompare Plots
mypdf(file="GraphCompare.pdf",width=6,height=5)
#X11(width=6,height=5)
#layout(matrix(c(0,0,0,1,0,2,0,0), 4, 2, byrow = TRUE),
#  c(0.05,0.9),c(0.05,0.45,0.45,0.01))
#layout.show(2)

## set up labels etc
xlabel <- "Temperature (GK)"
ylabel <- "Reaction rate ratio"


# set up the parameters for plotting
oldsettings <-  par(mar=c(4,5,2.5,2)+0.1,yaxs="i",xaxs="i")

# find the log parts of the yaxis
#YAxisLims <- c(min(CompRate[,2:3],LitCompRate[,2:3]),
#               max(CompRate[,2:3],LitCompRate[,2:3]))
# include literature uncertainty band if present
if(length(LitData[1,]) > 2){
  YAxisLims <- c(min(CompRate[,2:3]),#,LitCompRate[,2:3]),
               max(CompRate[,2:3])) #,LitCompRate[,2:3]))
} else {
  YAxisLims <- c(min(CompRate[,2:3],na.rm=TRUE), max(CompRate[,2:3],na.rm=TRUE))
}                

## #############
## # Now do the second plot, comparisons
## 
## # First the temperature grids have to be matched
matchlist1 <- as.double(na.omit(match(LitData[,1],MyRate[,1])))
matchlist2 <- as.double(na.omit(match(MyRate[,1],LitData[,1])))

## if no uncertainties are on the literature, then column 2 is the recommended
if(length(LitData[1,]) > 2){
    reccol <- 3
} else {
    reccol <- 2
}
##RateComp <- cbind(MyRate[matchlist1,1],
##                  MyRate[matchlist1,2]/LitData[matchlist2,reccol],
##                  MyRate[matchlist1,3]/LitData[matchlist2,reccol],
##                  MyRate[matchlist1,4]/LitData[matchlist2,reccol],
##                  MyRate[matchlist1,5],MyRate[matchlist1,6])
if(NewOverLiterature){
RateComp <- cbind(MyRate[matchlist1,1],
                  MyRate[matchlist1,2]/LitData[matchlist2,3],
                  MyRate[matchlist1,3]/LitData[matchlist2,3],
                  MyRate[matchlist1,4]/LitData[matchlist2,3],
                  fulldata[matchlist1,9],fulldata[matchlist1,10])
} else {
RateComp <- cbind(MyRate[matchlist1,1],
                  LitData[matchlist2,2]/MyRate[matchlist1,3],
                  LitData[matchlist2,3]/MyRate[matchlist1,3],
                  LitData[matchlist2,4]/MyRate[matchlist1,3],
                  fulldata[matchlist1,9],fulldata[matchlist1,10])
}
cYAxisLims <- range(RateComp[,c(2:4)])
YAxisLims <- range(c(cYAxisLims,YAxisLims))

if(!is.null(YRangeUser)){
    YAxisLims <- YRangeUser
}

##YAxisLims <- c(1e-1,1e1)
minBase10 <- floor(log10(YAxisLims[1]))
maxBase10 <- ceiling(log10(YAxisLims[2]))

poly <- function(T,low,high){

    shape <- rbind(cbind(T,high),
#                   c(max(T),high[length(high)]),
                   cbind(rev(T),rev(low)),
                   c(min(T),high[1]))

    return(shape)
}


# plot the axes
plot(1,1,ylim=c(YAxisLims[1]*0.96,YAxisLims[2]*1.04),
     xlim=c(0.01,TMax+.001),
     type="n",log="xy",xlab="",ylab="",
     lwd=2,xaxp=c(0.01,TMax,1),yaxt="n")#yaxp=c(10^minBase10,10^maxBase10,1))
aY <- 10^seq(from=minBase10,to=maxBase10,by=1)
axis(4,at=aY,labels=FALSE)

# Plot the low and high rate polygon
if(NewOverLiterature){
    polygon(poly(T=RateComp[,1],low=RateComp[,2],high=RateComp[,4]),
	    col=add.alpha(CompColor,CompTrans),lty=0)
} else {
    polygon(poly(T=CompRate[,1],low=CompRate[,2],high=CompRate[,3]),
	    col=add.alpha(RateColor,RateTrans),lty=0)
}
    
# add line for upper rate
lines(RateComp[,1],RateComp[,3],lwd=2)
# Add lines for the literature uncertainty bands if present
#if(length(LitData[1,]) > 2){
#  lines(LitCompRate[,1],LitCompRate[,2],lwd=2,lty=2)
#  lines(LitCompRate[,1],LitCompRate[,3],lwd=2,lty=2)
#}
# add line through 1
lines(c(0.01,10),c(1,1),lty=3,lwd=2)

## now, if the y-range is small, add some small labels
if(diff(log10(YAxisLims))<2){
    majory <- pretty(YAxisLims)
    axis(2,at=majory,tcl=0.4)
}else{
    axis(2,at=aY,labels=axTexpr(2,at=aY))
    
}

# Now add the minor tick marks
majorx <- axTicks(1,axp=c(0.01,10,1))
majory <- aY#axTicks(2,axp=c(10^minBase10,10^maxBase10,1))
#majory <- c(majory[1],majory,majory[length(majory)]*10)

minorx<-array(0,dim=10*(length(majorx)-1))
minory<-array(0,dim=10*(length(majory)-1))
# now populate the minor tick marks
for(i in 1:(length(majorx-1))){
  minorx[(1+(i-1)*10):(i*10)] <- seq(majorx[i],majorx[i]*10,majorx[i])
}
axis(1,at=minorx,labels=FALSE,tcl=0.4)
# for y-ticks, extend further
for(i in 1:(length(majory-1))){
  minory[(1+(i-1)*10):(i*10)] <- seq(majory[i],majory[i]*10,majory[i])
}
axis(2,at=minory,labels=FALSE,tcl=0.4)
axis(4,at=minory,labels=FALSE,tcl=0.4)


#output2D("Uncerts.ps")


## 
## minBase10 <- floor(log10(min(RateComp[,2:4])))
## maxBase10 <- ceiling(log10(max(RateComp[,2:4])))
## 
## # plot the lower rate
if(NewOverLiterature){
    polygon(poly(T=LitCompRate[,1],low=LitCompRate[,2],high=LitCompRate[,3]),
	    col=add.alpha(RateColor,RateTrans),lty=0)
} else {
    polygon(poly(T=RateComp[,1],low=RateComp[,2],high=RateComp[,4]),
	    col=add.alpha(CompColor,CompTrans),lty=0)
}##lines(RateComp[,1],RateComp[,3],lwd=1,lty=2)

## Add the axis labels
xpos <- grconvertX(0.5,from="npc",to="user")
ypos <- grconvertY(0.05,from="ndc",to="user")
text(xpos,ypos,"Temperature (GK)",cex=1.7,xpd=TRUE)
xpos <- grconvertX(0.06,from="ndc",to="user")
ypos <- grconvertY(0.5,from="npc",to="user")
text(xpos,ypos,"Reaction Rate Ratio",cex=1.7,srt=90,xpd=TRUE)
xpos <- grconvertX(0.5,from="ndc",to="user")
ypos <- grconvertY(0.95,from="ndc",to="user")
text(xpos,ypos,parse(text=ReacName),cex=1.7,xpd=TRUE)

box(which="plot")
#output2D("Comparison.ps")
dev.off()
write.table(RateComp,"GraphCompare.dat",col.names=FALSE,row.names=FALSE,sep="\t")


