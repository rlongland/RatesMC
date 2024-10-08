name="RatesMC.full"
samplename = "RatesMC.samp"
ntemps = 61

lognorm_norm <- function(){
    sum(lognorm(hist[,1]))
    max(hist[,3])
}

lognorm <- function(x){
  (1/(sigma*sqrt(2*pi))) * exp(-((log(x)-mu)^2)/(2*sigma^2))
}


######################################################
# MAIN PROGRAM

#pdf(file="RatesMC.pdf",paper="letter")
#postscript(file="RatesMC.pdf",paper="letter")
pdf(file="GraphPanelall.pdf",width=8.25,height=10.75)

layout(matrix(c(1:9), 3, 3, byrow = TRUE))
#layout.show(9)

#settings <- par(mar=c(2.1,5.1,2.1,1.1))

Samples = as.double(read.table(name,skip=2,header=FALSE,nrows=1)[3])
data <- read.table(name,skip=4,header=FALSE,stringsAsFactors=FALSE)
NumberPlot<-1

# set up the options for error handling
#eofError <- simpleError("Reached the end of file")
#old.op <- options(error=function() eofError)

pb <- txtProgressBar(min=1,max=ntemps,style=3)
cat("   Reading samples...\n")


for( i in 1:ntemps ) {
  setTxtProgressBar(pb, i)
  
##  samples <- withRestarts(read.table(samplename,skip=i*3+(i-1)*Samples,
##                                  nrows=Samples,header=FALSE),
##                       abort=function(){ })$V1
  samples <- withRestarts(scan(samplename,what=double(),nmax=Samples,
                               skip=1+(i*3+(i-1)*Samples),quiet=TRUE),
                          abort=function(){ })
  # check to make sure hist was read
  if(is.null(samples) | length(samples)==0){
    cat("\nReached the end of file\n")
    break
  }

  histt <- hist(log10(samples),breaks=1000,plot=FALSE)
  hist <- cbind(histt$mids,histt$mids,histt$density/100)
  ##  hist <- read.table(histname,skip=templist[i]*3+(templist[i]-1)*1000,
  ##                     nrows=1000,header=FALSE)
  ##hist[,3]<-(hist[,3])/sum(hist[,3])

  # get rid of brackets
  data[i,9] <- gsub("([()])","",data[i,9])
  data[i,10] <- gsub("([()])","",data[i,10])
  
  mu <- as.double(data[i,9])
  sigma <- as.double(data[i,10])

  mean <- exp(mu+(sigma^2)/2)
  var <- (exp(sigma^2) -1)*exp((2*mu + sigma^2))

  # make histogram
  ##plot(hist[,1],hist[,3]/lognorm_norm(),type="l",col="tomato1",
  ##     ylim=c(0,max(hist[,3])/lognorm_norm()),
  ##     xlab="",ylab="",cex.axis=1.3,yaxs="i")
  plot(hist[,1],hist[,3],type="l",col="tomato1",
       xlab="",ylab="",cex.axis=1.3,yaxs="i")

  # Make the lognormal fit
  xlims <- c(min(hist[,1]),max(hist[,1]))
  x <- seq(from=log10(min(samples[samples>0])),
	   to=log10(max(samples[samples>0])),
	   length.out=length(hist[,3]))
  y <- lognorm(10^x)
  ##y <- y/max(y)
  y <- sum(hist[,3])*y/sum(y)
  
  ## The normalisation
  ##n <- mean((hist[,3]/y)[(hist[,3])>0.5*max(hist[,3])])
  n <- 1
  lines(x,y*n,lwd=1.5)

#  lines(x,lognorm(x)/lognorm_norm(),lwd=1.5)

#  xpos<-par()$usr[1]+(8/10)*(par()$usr[2]-par()$usr[1])
#  ypos<- (8/10)*max(hist[,3])
  xpos <- grconvertX(0.75,from="npc",to="user")
  ypos <- grconvertY(0.85,from="npc",to="user")
  
  AD_text <- paste("T9 =",data[i,1], "\nA-D =",format(data[i,11],digits=3))

  text(xpos,ypos,AD_text,cex=1.5)


  #  open a new window when the last was full
  if(i%%9 == 0){
#    X11()
#    dev.copy(device=x11)
#    dev.set(which=dev.prev())
    layout(matrix(c(1:9), 3, 3, byrow = TRUE))
    NumberPlot <- NumberPlot+1
  }

#  output2D(name=paste("Plot",NumberPlot,".eps",sep=""))

}

close(pb)
dev.off()

#restore error handling
#options(old.op)

#dev.copy(X11)


fit <- function(){
  fn <- function(p) sum(y - (1/(x*p[2]*sqrt(2*3.142))) * exp(-((log(x)-p[1])^2)/(2*p[2]^2)))
  p <- c(mu,sigma)
  y <- hist[,3]
  x <- hist[,1]

  out <- nlm(fn,p)

  mu <- out$estimate[1]
  sigma <- out$estimate[2]

  lines(x,lognorm(x)*max(hist[,3])/max(lognorm(x)),col="blue")
}
  
