# Plotting.R produces the plots for ABangle.
# It takes in the data to be plotted through a pipe.

# Read from the pipe.
data <- read.table("stdin", header=TRUE)

# Find the colour with which to plot the set
Colours <- data$cols[data$cols != "nan"]

# Find the number of sets that have been passed:
nsets <- length( Colours[ !is.nan(Colours) ])

# Find which angles are being plotted
anglenames <- data$anglenames[data$anglenames != "nan"]

# Find which type of plot each set needs to be plotted with
plottype <- data$plottype[data$plottype != "nan"]

# Find the line type for each of the plots.
Ltype <- data$Ltype[ !is.nan(data$Ltype) ]

# Open the file name of the image
Ylim=list()
for ( angle in anglenames ){
		X <- data[ paste( angle, "Full",sep="") ][,]
		X <- X[!is.nan(X)]
		FROM=min(X)
		TO=max(X)
		Ylim[[angle]] <- c(0, max( hist(X, breaks=50, plot=FALSE)$density ))
		}


repeat{
	finish=TRUE
	png( paste( data$fname[data$fname != "nan"],sep="" ), width=720, height=480 )
	if (length( anglenames ) != 1){ par(mfcol=c(2,3) )}
	par(cex.lab=1.5, cex.main=2, cex.axis=1.2)
	# for each of the angles plot the data
	for ( angle in anglenames ){
		# Draw the background distribution
		X <- data[ paste( angle, "Full",sep="") ][,]
		X <- X[!is.nan(X)]
		FROM=min(X)
		TO=max(X)
		ymax<-list()
		print(Ylim[[angle]])
		if (angle == "dc")
			{
				h<-hist( X ,breaks=50,main=paste( angle ),xlab="Distance (Angstroms)",freq=FALSE,ylim=Ylim[[angle]])
				ymax[[angle]] <- c( max( h$density ))
			}
		else
			{
				h<-hist( X ,breaks=50,main=paste( angle ),xlab="Angle (Degrees)",freq=FALSE,ylim=Ylim[[angle]] )
				ymax[[angle]] <- c( max( h$density ))
			}

		# For each of the angle sets plot the angles using the appropriate plotter and color
		for (i in 1:nsets){
			if (plottype[i] == "density"){
				x <- data[ paste( angle, i, sep="")][,]
				x <- x[!is.nan(x) ]
				d<-density( x, from=FROM, to=TO )
				lines(d , col=paste(Colours[i]), lwd=3, lty=Ltype[i] )
				ymax[[angle]] <- c( ymax[[angle]], max( d$y) )
			}


		
			if (plottype[i] == "lines"){
				x <- data[ paste( angle, i, sep="")][,]
				x <- x[!is.nan(x) ]
				for (j in 1:length(x) ){
					lines( c( x[j], x[j]), c(0,1000), col=paste(Colours[i]), lwd=3, lty=Ltype[i] )
					}
			}

	
			}
			if (Ylim[[angle]][2] != max(ymax[[angle]]) )
			{Ylim[[angle]][2] = max(ymax[[angle]])
			finish=FALSE}
		}

	dev.off()

	if (finish)
	{break}


}


if ( nsets != 0){
	Legendname <- paste( sub( ".png", "legend.png", paste( data$fname[data$fname != "nan"],sep="" )), sep="", collapse="")
	png( Legendname, width=600)
		plot( c(0,1), c(0,1), type="n",tck=0,xaxt='n',yaxt='n',xlab="",ylab="",axes=FALSE)

		Labels = c()
		for (i in 1:nsets){
			Labels <- c( Labels, paste(sub( "~","'", sub("/", "\n", sub("@", " ", paste( data[ paste( "Label", i, sep="" ) ][ data[ paste( "Label", i, sep="" ) ] != "nan"], sep="") ))), sep="", collapse="") )
			}

		legend( 0,0.9,legend=Labels, lty=Ltype,lwd=3,col=paste(Colours),cex=1,text.col=paste(Colours),bty="n")

		dev.off()
	}












 





