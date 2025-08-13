######################################################################
## PlotContribution.R
## Plot the contributions of individual resonances to the total rate
##
## Still not "production-ready" - labels need to be tweaked and tested
######################################################################

RatesMCFile      <- "RatesMC.in"
ContributionFile <- "RatesMC.cont"
TMatchFile       <- ""   ## set "" to ignore Tmatch file

extraSuper <- ""            ## Use to put a superscript at the end
                            ## e.g. "g" for ^{26}Al^g

TRange <- c(0.01,10)        ## Tempterature range in GK
Threshold <- 0.05           ## Only look at resonances that contribute
                            ## more than 'Threshold' over this
                            ## temperature range


randomColour <- FALSE        ## Randomize the colours


######################################################################
## You shouldn't need to touch anything below this line...
##library("RColorBrewer")

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

if(file.exists(TMatchFile)){
  title <- scan(RatesMCFile,what=character(),n=1)
  titleTMatch <- scan(TMatchFile,what=character(),n=1)
  if(title != titleTMatch){
    cat("\n",
        "Input file for reaction ",title,"\ndoesn't match Tmatch file ",
        "for ",titleTMatch,"!\n",sep="")
  } else {    
    TMatch <- read.table(TMatchFile,skip=1,nrow=1)[11]
    TRange[2] <- TMatch
  }
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
## Add a superscript
if(nchar(extraSuper)>0)ReacName <- paste(ReacName,"^",extraSuper,sep="")

## Find the beginning of the resonances
skip <- 28
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
}
## Now read the upper limit resonances in the same way
skip <- skip+5
while(
  substr(x <- try(unlist(read.table(RatesMCFile,skip=skip,nrows=1)[1]
                             )),1,2) != "**"){
  skip <- skip+1
  if(is.numeric(x))
    Energies <- c(Energies,x)
}
Energies <- format(Energies,digits=1,trim=TRUE)

## Look for resonances that contribute more than 5% in the specified
## temperature range
## Look only at the centroids
MedSeq <- seq(from=2,to=length(data[1,]))
rangeBool <- data[,1]>=TRange[1] & data[,1]<=TRange[2]
T <- data[rangeBool,1]
resCont <- data[rangeBool,MedSeq]

## toplot is the list of resonances to plot
toplot <- resCont##[rangeBool,]

## Plot tempteratures on a log scale
logT <- log10(T)

## Find the resonances that are contributors over this range
contribBool <- apply(resCont,2,function(x)any(x>Threshold))
## There is repetition here, so pare this down to resonance numbers
contributors <- unique(ceiling(which(contribBool)/3))
## Also find the resonances that don't contibute significantly
noncontributors <- unique(ceiling(which(!contribBool)/3))
noncontributors <- noncontributors[!sapply(
  noncontributors,function(x)any(x==contributors))]

## order the contributors according to their maximum contribution
maxCont <- apply(toplot,2,max)
contributors <- contributors[order(maxCont[contributors*3-1],decreasing=TRUE)]

## Write a list of the important resonance names for putting on the plot
## Remember that the first two are the analytical DC rate
resNames <- sapply(contributors,
                function(i){
                  if(i>2 & i<=length(Energies)+2)
                    t <- paste(Energies[i-2]," keV",sep="")
                  if(i>length(Energies)+2)
                    t <- paste("Intf ",i-(length(Energies)+2),sep="")
                  if(i==1)
                    t <- "DC"
                  if(i==2)
                    t <- paste("DC ",i,sep="")
                  t
                })


## The list of colours
lc <- length(contributors)
hue <- c(seq(from=0,to=1-2/lc,length.out=floor(lc/2)),
         seq(from=1/lc,to=1-1/lc,length.out=ceiling(lc/2)))
sat <- 0.8
val <- 0.7
col <- hsv(hue,sat,val,alpha=1)
if(randomColour){
    col <- sample(col)
}else if(length(col)>2){
    col <- col[unlist(lapply(1:3,function(i)seq(i,length(col),by=3)))]
}

density <- c(NA,rep(1.5,length(contributors)))*20
angle <- rep(seq(0,270,45),ceiling(length(contributors)/7))


## The contributions plot
mypdf("GraphContribution.pdf",width=6,height=5)
#while(dev.cur()>1)dev.off()
#myX11(width=7,height=5)

oldpar <- par(mar=c(4,5,2.5,1)+0.1)

maxy <- 1.38
plot(range(logT),c(0,maxy),type="n",
     ylim=c(0,maxy),yaxt='n',xaxt='n',xaxs='i',yaxs='i',
     xlab="",ylab="Fractional contribution")

abline(h=1,lty=2)


## Add the lines all in one go.
##matlines(logT,toplot,type="l",lty=1,lwd=2)
icol <- 1
plotpoly <- function(i){
  sel <- seq(from=i*3-2,to=i*3)
  x <- toplot[,sel]
  wid <- 0 #0.005
  shape <- rbind(cbind(logT,x[,3]+wid),
                 c(max(logT),x[length(x[,1]),1]),
                 cbind(rev(logT),rev(x[,1])-wid),
                 c(min(logT),x[1,3]))

  ## Calculate the overlap between this resonance and the ones that
  ## have already been plotted

  ol <- sapply(contributors[(1:length(contributors))<icol],function(j){
    ## contribution j
    csel <- seq(from=j*3-2,to=j*3,by=1)
    cx <- toplot[,csel]
    ## overlap amount
    ol <- apply(cbind(x[,3],cx[,3]),1,min)-apply(cbind(x[,1],cx[,1]),1,max)
    ## just count number of temperature points that overlap
    overlap <- sum(ol>0.05)
    overlap
  })

  ## Decide that significant overlap occurs with olthresh
  olthresh <- 3
  idensity <- length(ol[ol>olthresh])+1
  ## Add transparency if no overlap
  if(idensity==1)
    col <- add.alpha(col,0.7)

  ## Draw the actual contribution
  polygon(shape,col=col[icol],border = col[icol],
          angle=angle[icol],density=density[idensity],lwd=1.5)
  ##  cat(i,"     ",ol,"     ",idensity,"\n",col[icol],"\n")
  icol <<- icol+1
}
## plot a polygon for each contributing resonance
iret <- sapply(contributors,plotpoly)

## Add a line for everything else
if(length(noncontributors)>1){
  others <- rowSums(toplot[,(noncontributors*3-1)])
} else {
  others <- toplot[,(noncontributors*3-1)]
}
if(length(others)>0)lines(x=logT,y=others,lty=3)

## Need to plot the labels
xxmin <- grconvertX(0.1,from="npc",to="user")
xxmax <- grconvertX(0.9,from="npc",to="user")
yymin <- 1.03
yymax <- maxy-0.03

## Make a matrix of used label coordinates
stepx <- grconvertX(0.04,from="npc",to="user")-
    grconvertX(0,from="npc",to="user")
lx <- seq(from=xxmin,to=xxmax,by=stepx)
stepy <- 0.06
ly <- seq(from=yymax,to=yymin,by=-stepy)
usedCoord <- matrix(FALSE,ncol=length(lx),nrow=length(ly))

## Complicated function to figure out where to place the text
plotCoord <- sapply(contributors,function(i){
  sel <- i*3-1
  ## Find the temperature that corresponds to the max contribution for
  ## this resonance
  ixx <- which.min(abs(logT[which.max(toplot[,sel])]-lx))
  xx <- lx[ixx]

  for(offset in c(0,1,-1)){
      ixx.offset <- ixx+offset
      xx <- lx[ixx.offset]
      ## Clear out a range around this position to that text doesn't overlap
      clearrange <- 3
      iixx <- seq(from=max(1,ixx.offset-clearrange),
		  to=min(length(lx),ixx.offset+clearrange),by=1)

      ## Find the row to put this text that doesn't already have something
      iyy <- which(!usedCoord[,ixx.offset])[1]
      yy <- ly[iyy]
      if(!is.na(yy)) break

  }
  ## Mark usedCoord to make sure future text isn't written here
  usedCoord[iyy,iixx] <<- TRUE

  if(is.na(xx) | is.na(yy)){
    cat("\033[91m\n WARNING!:\n",
        " Resonance label (",resNames[which(contributors==i)],") couldn't be plotted. There's no room!\n",
	"\033[0m",sep="")
    usedCoord[iyy,iixx] <<- FALSE
    ##print(xx)
    ##print(yy)

  }
  ## return the chosen text coordinate
  c(xx,yy)
})
## Plot the text in this position
text(plotCoord[1,],plotCoord[2,],resNames,
     col=add.alpha(col[1:length(resNames)],1),pos=1,offset=0)

## The x-axis. Allow for small ranges...
if(max(T)/min(T) < 100){   ## Small Ranges
  pt <- pretty(T)
  aXM <- log10(pt)
  nminor <- 4
  ptt <- seq(from=min(pt),to=max(pt),
             length.out=length(pt)*nminor-(nminor-1))
  aXm <- log10(ptt)
} else {  ## Large ranges
  aXM <- c(0.01,0.1,1,10)
  aXm <- array(,dim=10*length(aXM))
  for(i in 1:(length(aXM))){
    aXm[(1+(i-1)*10):(i*10)] <- seq(aXM[i],aXM[i]*10,aXM[i])
  }
  aXM <- log10(aXM)
  aXm <- log10(aXm)
}
axis(1,at=aXM,labels=10^aXM)
mtext("Temperature (GK)",side=1,line=2.3,cex=1.5)
axis(1,at=aXm,labels=FALSE,tcl=0.3)
##axis(3,at=aXM,labels=FALSE)
##axis(3,at=aXm,labels=FALSE,tcl=0.3)

## y-axis
axis(2,at=c(0,0.2,0.4,0.6,0.8,1.0))

## Plot the title
xpos <- grconvertX(0.5,from="npc",to="user")
ypos <- grconvertY(1,from="ndc",to="user")
text(xpos,ypos,parse(text=ReacName),cex=1.7,xpd=TRUE,pos=1,offset=0.7)


## Stick a box around the plot to clean up
box()

dev.off()

