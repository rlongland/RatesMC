######################################################################
## PlotIntegrand.R
## Plot the broad resonances from the integration routines
##
## Need to add vertical lines for narrow resonances
######################################################################

RatesMCFile   <- "RatesMC.in"
IntegrandFile <- "RatesMC.integ"

Temp <- 0.1  ## The temperature to plot the integrand at

## Set limits to NA for auto plotting
xlim <- c(0,1)
ylim <- c(1e-50,1e-15)

######################################################################

## Read all lines to look for temperature
lines <- readLines(IntegrandFile)
TLines <- grep("Temp",lines)
## Which temperature entry is the one we want 
Tentry <- which(read.table(text=lines[TLines])[,3] == Temp)
## and corresponds to line:
goodTLine <- TLines[Tentry]
## The next temperature entry is
nextTLine <- TLines[Tentry+1]

## The list of lines with integrand data
ilines <- (goodTLine+1):(nextTLine-2)

## Now each resonance is separated by spaces
isep <- c(1,grep("^$",lines[ilines]))
## Now we can read the integrand for each resonance we encounter
integ <- lapply(2:length(isep), function(x){
    data <- read.table(text=lines[ilines[isep[x - 1]:(isep[x]-1)]])
    data
})

## Close any open windows
while(dev.cur()>1)dev.off()
myX11()

## First make the plotting pane
ylim.default <- range(unlist(lapply(integ, function(x)range(x[,2]))))
if(is.na(ylim[1]))ylim <- ylim.default
xlim.default <- c(0,10)
if(is.na(xlim[1]))xlim <- xlim.default

plot(xlim,ylim,ylim=ylim,type='n',log='y',
     xlab="E (MeV)", ylab="Rate Integrand (arb. units)")

## The list of colours
lc <- length(integ)
hue <- c(seq(from=0,to=1-2/lc,length.out=floor(lc/2)),
         seq(from=1/lc,to=1-1/lc,length.out=ceiling(lc/2)))
sat <- 0.8
val <- 0.7
col <- hsv(hue,sat,val,alpha=1)
icol <- 1

## Plot each line
t <- lapply(integ, function(data){
    ## Sort data by energy
    ord <- order(data[,1])
    T9 <- data[ord,1]
    y <- data[ord,2]
    ## Finally plot
    lines(T9,y,col=col[icol])
    icol <<- icol+1
})

legend(x="topright",legend=paste("Res",1:length(integ)),col=col,lty=1)
