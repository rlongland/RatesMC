######################################################################
## makeLaTeX.R
## Read in the LaTeX output, reformat the table
##
## Author: R. Longland
## Year  : 2023
######################################################################
filename <- "RatesMC.latex"
filecontent <- readChar(filename, file.info(filename)$size)
filecontent <- gsub("\n", "", filecontent)
filecontent <- gsub("\\\\\\\\", "\n", filecontent)

##writeChar(filecontent,"test.dat")
con <- textConnection(filecontent)
dat <- read.table(con,sep="&",nrows = 60, colClasses="character")

cut1 <- 1:30
cut2 <- 31:60
dat2c <- cbind(dat[cut1,],dat[cut2,])

##dat2c.format <- format(dat2c,nsmall=3)

write.table(file="RatesMC.latex2",dat2c,quote=F,sep=" & ",
	    eol="  \\\\ \n",col.names=F,row.names=F)

##----------------------------------------------------------------------
filename <- "TMatch.latex"
if(file.exists(filename)){
    filecontent <- readChar(filename, file.info(filename)$size)
    filecontent <- gsub("\n", "", filecontent)
    filecontent <- gsub("\\\\\\\\", "\n", filecontent)

    ##writeChar(filecontent,"test.dat")
    con <- textConnection(filecontent)
    dat <- read.table(con,sep="&",nrows = 60, colClasses="character")

    cut1 <- 1:30
    cut2 <- 31:60
    dat2c <- cbind(dat[cut1,],dat[cut2,])

    ##dat2c.format <- format(dat2c,nsmall=3)

    write.table(file="TMatch.latex2",dat2c,quote=F,sep=" & ",
		eol="  \\\\ \n",col.names=F,row.names=F)
}
