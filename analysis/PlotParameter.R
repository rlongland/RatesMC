######################################################################
## PlotCorrelations.R
## Plot correlation between rate and input parameters
##
######################################################################

ParameterFile <- "ParameterSamples.dat"

######################################################################

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

##--------------------------------------------------
## Analyzing and Plotting
##--------------------------------------------------
## Close any open windows
##while(dev.cur()>1)dev.off()
##myX11()
pdf("GraphParameters.pdf",width=8.25,height=10.75)

layout(matrix(c(1:9), 3, 3, byrow = TRUE))

for(i in 1:length(names.cut)){

    hist(pars.cut[,i],500,main=names.cut[i],xlab="",ylab="")

    var.par <-var(pars.cut[,i])
    mean.par <- mean(pars.cut[,i])

    cat("Par",names.cut[i],": Mean =",mean.par, " SD =",sqrt(var.par),"\n")
    xx <- grconvertX(0.9,from="npc",to="user")
    yy <- grconvertY(0.9,from="npc",to="user")

    mtext(bquote(.(format(mean.par,digits=4)) ~ "+/-" ~ .(format(sqrt(var.par),digits=4))),
	  side=1,line=2,cex=0.8)
    
    if(i%%9 == 0){
	layout(matrix(c(1:9), 3, 3, byrow = TRUE))
    }
}

dev.off()
