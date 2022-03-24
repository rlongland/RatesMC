######################################################################
## PlotCorrelations.R
## Plot correlation between rate and input parameters
##
######################################################################
library(grid)

RatesMCFile   <- "RatesMC.in"
SampleFile    <- "RatesMC.samp"
ParameterFile <- "ParameterSamples.dat"

Temp <- 0.04  ## The temperature to look at

## Set limits to NA for auto plotting
xlim <- NA ##c(0,1)
ylim <- NA ##c(1e-50,1e-15)

######################################################################
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
##--------------------------------------------------
## Analyzing and Plotting
##--------------------------------------------------
## Close any open windows
##while(dev.cur()>1)dev.off()
##myX11()
mypdf("GraphCorrelation.pdf",onefile=TRUE)

## Cut the data to only include 100 samples
NSamples <- 100
samp <- sample(1:dim(dd)[1],NSamples)

## Modify the names with the resonance
## nE <- length(names.cut[names.cut=="E"])
names.cut[names.cut=="E"] <- "E["
## nwg <- length(names.cut[names.cut=="wg"])
names.cut[names.cut=="wg"] <- "omega*gamma["
## nG1 <- length(names.cut[names.cut=="G1"])
names.cut[names.cut=="G1"] <- "Gamma[a"
## nG2 <- length(names.cut[names.cut=="G2"])
names.cut[names.cut=="G2"] <- "Gamma[b"
## nG3 <- length(names.cut[names.cut=="G3"])
names.cut[names.cut=="G3"] <- "Gamma[c"
ires <- 0
for(i in 5:dim(dd)[2]){
    if(length(grep("E",names.cut[i]))>0)ires  <- ires + 1
    names.cut[i] <- sprintf("%s%d]",names.cut[i],ires)
}
dd.names <- c("Rate",names.cut)

## First just plot 3x3 matrix of Rate vs. variable
layout(matrix(c(1:9), 3, 3, byrow = TRUE))
for(i in 2:dim(dd)[2]){

    lab <- names.cut[i-1]
    plot(x=dd[samp,i], y=dd[samp,1], xlab=parse(text=lab), ylab="log10(Rate)")

    c <- cor(dd[,1], dd[,i],use="complete.obs")
    col <- ifelse(abs(c)>0.2, "red", "black")
    mtext(side=3,paste("Corr:",format(c,digits=3)), line=0,
	  col=col)
}


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

if(dim(dd)[2] < 9){
    pairs(dd[samp,],labels=dd.names, lower.panel=panel.cor, pch=21,
	  col=NA, bg=add.alpha("red",0.4),gap=0,
	  yaxt='n',xaxt='n')
} else {
    seq <- seq(from=6,to=dim(dd)[2],by=7)
    for(i in 1:length(seq)){
	end <- min(seq[i]+6,dim(dd)[2])
	
	sub <- c(1:5,seq[i]:end)
	pairs(dd[samp,sub],labels=parse(text=dd.names[sub]), lower.panel=panel.cor, pch=21,
	  col=NA, bg=add.alpha("red",0.4),gap=0,
	  yaxt='n',xaxt='n')
    }
}

dev.off()
