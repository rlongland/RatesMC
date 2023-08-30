name="RatesMC.full"
samplename = "RatesMC.samp"
ntemps = 6
templist = c(11,14,18,28,33,37)

extraSuper <- ""            ## Use to put a superscript at the end
                            ## e.g. "g" for ^{26}Al^g


lognorm_norm <- function(){
  sum(lognorm(hist[,1]))
  max(hist[,3])
}

lognorm <- function(x){
  (1/(x*sigma*sqrt(2*3.142))) * exp(-((log(x)-mu)^2)/(2*sigma^2))
}


######################################################
# MAIN PROGRAM

#X11(width=8.25,height=10.75)
#open the postscript file
#postscript(file="Panel.ps",horizontal=F, onefile=F)
pdf(file="GraphPanel6.pdf",width=8.25,height=10.75,onefile=F)

#layout(matrix(c(1,1,1,9,3,4,9,5,6,9,7,8,2,2,2), 5, 3, byrow = TRUE),
#  c(0.1,2,2),c(0.5,3,3,3,0.5))
#layout.show(9)
layout(matrix(c(0,0,0,0,1,2,0,3,4,0,5,6,0,0,0), 5, 3, byrow = TRUE),
  c(0.1,2,2),c(0.5,3,3,3,0.5))
layout.show(6)

# set the plot parameters to better arrange the plots
# mar=c(bottom,left,top,right)
settings <- par(mar=c(2.1,5.1,2.1,1.1))

## Read and format the reaction name
ReacName = as.character(read.table(name,skip=0,header=FALSE,nrows=1)$V1)
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

FirstLetter <- regexpr("([a-z])",ReacName)[1]
word1<-substr(ReacName,1, FirstLetter-1)
word2<-substr(ReacName,FirstLetter, 100000)
GamPos <- regexpr(",g",word2)[1]
if(GamPos > 0) {
  word2 <- substr(word2,1,GamPos)
  word3 <- "gamma"
}

Samples = as.double(read.table(name,skip=2,header=FALSE,nrows=1)[3])
data <- read.table(name,skip=4,header=FALSE,stringsAsFactors=FALSE)
NumberPlot<-1


for( i in 1:ntemps ) {

  samples <- withRestarts(read.table(samplename,skip=1+(templist[i]*3+
                                     (templist[i]-1)*Samples),
                                  nrows=Samples,header=FALSE),
                       abort=function(){ })$V1

  histt <- hist(samples,breaks=1000,plot=FALSE)
  hist <- cbind(histt$mids,histt$mids,histt$counts)
#  hist <- read.table(histname,skip=templist[i]*3+(templist[i]-1)*1000,
#                     nrows=1000,header=FALSE)
  hist[,3]<-(hist[,3])/sum(hist[,3])

  # get rid of brackets
  data[templist[i],9] <- gsub("([()])","",data[templist[i],9])
  data[templist[i],10] <- gsub("([()])","",data[templist[i],10])
  
  mu <- as.double(data[templist[i],9])
  sigma <- as.double(data[templist[i],10])

  mean <- exp(mu+(sigma^2)/2)
  var <- (exp(sigma^2) -1)*exp((2*mu + sigma^2))

  # make histogram
  #plot(hist[,1],hist[,3],type="l",col="red",log="y",ylim=c(0.1,max(hist[,3])))
  plot(hist[,1],hist[,3]/lognorm_norm(),type="l",col="tomato1",
       ylim=c(0,max(hist[,3]/lognorm_norm())),
       xlab="",ylab="",cex.axis=1.3,yaxs="i")
#  print(lognorm_norm())
  
  # Make the lognormal fit
  xlims <- c(min(hist[,1]),max(hist[,1]))
  x<-seq(from=xlims[1]/2,to=xlims[2]*2,length.out=2000)
#  lines(x,lognorm(x)/lognorm_norm(),lwd=1.5)
  h <- lognorm(hist[,1])/max(lognorm(hist[,1]))
  ## The normalisation
  n <- median((h/(hist[,3]/lognorm_norm()))[h>0.3])
  lines(x,lognorm(x)/(n*max(lognorm(x))))

  xpos <- grconvertX(0.75,from="npc",to="user")
  ypos <- grconvertY(0.85,from="npc",to="user")

  AD_text <- paste("T9 =",data[templist[i],1],"\nA-D =",data[templist[i],9])

  text(xpos,ypos,AD_text,cex=1.5)

}

#mtext(expression(paste("Reaction rate (cm"^{3}," mol"^{-1}," s"^{-1},")")),side=1,adj=-40,padj=1.95,cex=1.4)
#mtext("Probability (arb. units)",side=2,adj=11,padj=-26.3,cex=1.4)
#mtext(ReacName,side=3,padj=-39,adj=-0.7,cex=1.4)

#reset_to_01_coordinates()
par(xpd=NA) 
Minx <- par("usr")[1]
Maxx <- par("usr")[2]
Miny <- par("usr")[3]
Maxy <- par("usr")[4]

xpos <- Minx-0.25*(Maxx-Minx)
ypos <- Miny-0.2*(Maxy-Miny)
text(xpos,ypos,expression(paste("Reaction rate (cm"^{3}," mol"^{-1}," s"^{-1},")")),cex=2)
xpos <- Minx-1.48*(Maxx-Minx)
ypos <- Maxy+0.65*(Maxy-Miny)
text(xpos,ypos,"Probability (arb. units)",cex=2,srt=90)

xpos <- Minx-0.25*(Maxx-Minx)
ypos <- Maxy+2.555*(Maxy-Miny)
text(xpos,ypos,parse(text=ReacName),cex=2)


#output2D(name="Panel.eps")
dev.off()



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
  
