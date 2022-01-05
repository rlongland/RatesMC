######################################################################
## PlotCorrelations.R
## Plot correlation between rate and input parameters
##
######################################################################

RatesMCFile   <- "RatesMC.in"
SampleFile    <- "RatesMC.samp"
ParameterFile <- "ParameterSamples.dat"

Temp <- 0.1  ## The temperature to look at

## Set limits to NA for auto plotting
xlim <- NA ##c(0,1)
ylim <- NA ##c(1e-50,1e-15)

######################################################################

##--------------------------------------------------
## Section to read reaction rate samples
##--------------------------------------------------
## Read all lines to look for temperature
lines <- readLines(SampleFile)
TLines <- grep("T9",lines)
## Which temperature entry is the one we want 
Tentry <- which(read.table(text=lines[TLines])[,3] == Temp)
## and corresponds to line:
goodTLine <- TLines[Tentry]
## The next temperature entry is
nextTLine <- TLines[Tentry+1]

## The list of lines with rate sample data
ilines <- (goodTLine+2):(nextTLine-2)

## Finally we can read the rate samples at this T9
data <- read.table(text=lines[ilines])[,1]

##--------------------------------------------------
## Section to read parameter samples
##--------------------------------------------------
pars <- read.table(ParameterFile,skip=2)

## Cut the data if it's a parameter that's exactly zero
cut <- !apply(pars,2,function(x)all(x==0))
pars.cut <- pars[cut]

names <- read.table("ParameterSamples.dat",nrow=1,skip=1)
names <- names[names != '|']
names.cut <- names[cut]

## Combine rates and parameters
dd <- cbind(log10(data),pars.cut)
dd.names <- c("Rate",names.cut)
##--------------------------------------------------
## Analyzing and Plotting
##--------------------------------------------------
## Close any open windows
##while(dev.cur()>1)dev.off()
##myX11()
mypdf("GraphCorrelation.pdf")

## Cut the data to only include 100 samples
NSamples <- 100
samp <- sample(1:dim(dd)[1],NSamples)

##library("PerformanceAnalytics")
##chart.Correlation(dd[samp,], histogram=FALSE,
##		  pch=18, col=add.alpha("red",0.4),
##		  labels=dd.names,yaxt='n')

# panel.cor puts correlation in upper panels, size proportional to correlation
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use="na.or.complete"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r, col="red")
}

pairs(dd[samp,],labels=dd.names, upper.panel=panel.cor, pch=21,
      col=NA, bg=add.alpha("red",0.4),gap=0,
      yaxt='n',xaxt='n')

dev.off()
