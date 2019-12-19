##################################################################################################

# The following script is to reproduce Figure 3 in the manuscript. 
# The two panels are plotted in sequence, and saved in the default working directory.

##################################################################################################

# First, begin by installing the require packages.

# 'deSolve' is required to solve the differential equation
# 'shape' is required to include the specialised arros that denote the phase plane.

install.packages(c("deSolve", "shape"), dependencies = TRUE)

# load packages
library(deSolve)
library(shape)


# Set all the default parameters for the model
a <- 1.5; i <- 0.01; T <- 10; d <- 0.03; m <- 1


# To begin, it is first necessary to define all the derivative equations that need to be solved

#-----------------------------
# The derivative functions
#-----------------------------

# The counterfactual is only defined according to the dynamics of a logit function
eco_counter <- function(t, H, parms) {
  dH <- (((log(H/(100-H)))/100) * (1/a) * H )  
  list(dH, dy = dH,i)
}

# The mitigation scenario is include delayed (lagged) compensation as soon as the time of the simulation is longer than parameter T 
eco_mitigate <- function(t, H, parms) {
  if (t < T)
    lag <- 0
  else
    lag <- lagvalue(t - T)
  
  dH <- (((log(H/(100-H)))/100) * (1/a) * H ) - (i * H) + (m * i * lag *(1+d)^T)
  list(dH, dy = dH,i)
}


# Set the incremental values for the mutltiplier, m.
multiplier <- seq(1,10,by=0.1)

# Run a loop for each value for the multiplier
for (j in 1:length(multiplier)){

  # Select the valeu of 'm' from the j-th value in the set of multipliers
  m <- multiplier[j]

  # Define all the starting values for habitat cover
  # Note: These are integers (excluding stable equilibria) in order to define the plot surface
  startH <- seq(1,99,by=1)

  # Define the time sequence
  # Note: These are integers in order to define the plot surface
  times <- seq(0, 100, by = 1)

  # Create a blank output matrix that records whether 'zero loss' is achieved
  zeroloss <- matrix(NA,ncol=length(times),nrow=length(startH))

  # Run a loop that solves the differential equations for every starting level of habitat cover
  for (k in 1:length(startH)) {
    # Solve equations for counterfactual and mitigation scenarios  
    tcf <- dede(y = startH[k], times = times, func = eco_counter, parms = NULL, atol = 1e-10)
    tmit <- dede(y = startH[k], times = times, func = eco_mitigate, parms = NULL, atol = 1e-10)

    # Cap maximum habitat cover to 100% for counterfactual
    CF <- tcf[,2] 
    CF[which(CF=="NaN")] <- 100
    CF[which(CF>100)] <- 100
    
    # Cap maximum habitat cover to 100% for mitigation scenario
    Mit <- tmit[,2] 
    Mit[which(Mit=="NaN")] <- 100
    Mit[which(Mit>100)] <- 100

    # Record the mitgation debt for the time series to the blank output matrix
    zeroloss[k,] <- CF - Mit 
  }

  # This command creates a multidimensional array
  # If, for the specific value of 'm', zero loss is acheived for that starting value of H during that year, 
  # make the value in the array equal the value of 'm'. If not, use a dummy value of 99.
  nnl <- ifelse(zeroloss<=0,m,99)

  # If this is the first run of the loop, the matrix above is saved as a list called 'mult' (for multidimensional array),
  # if it is not the first loop, the matrix is turned into a list and then appended to the existing array called 'mult'.
  if (j == 1) {
	  mult <- list(nnl)
  } else{
	 mult <- append(mult,list(nnl))
  }
}


# Define the name of the plot (plus dimensions); alterantive directories can be identified here
png(filename="Fig3.png",width=10,height=20,units="cm",res=300)

# Set plot margins and layouts
par(mai=c(0.8,0.8,0.18,0.1))
par(mfrow=c(2,1))



###############################################
#                                             #
#                 Panel A                     #
#                                             #
###############################################

# This command takes the multidimensional array ('mult') and determins the minimum value of 'm' that would allow for 
# 'zero loss' by each year for every starting level of H
min.m <- apply(simplify2array(mult), 1:2, min)

# Because 'zero loss' is always possible in the first year (before impact), we set this to the first year values to the dummy value of 99.
min.m[,1] <- 99

# The panel is actually a scatterplot, so define the x- and y-values for each combination of time and starting habitat
yval <- rep(startH,each=length(times))
xval <- rep(times,length(startH))

# Convert the output matrix to a vector (z-values) that will define the colour of the minimum multiplier
zval <- as.vector(t(min.m))

# Define a colour ramp 
ramp <- colorRampPalette(c("grey90","blue","darkblue","black"),bias=1.9,interpolate="linear")(length(multiplier))


# This makes the plot for each combination of time and starting habitat, colouring the points according to the mitigation debt
plot(xval+2020,yval,pch=15,cex=0.7, las=1, cex.axis=1, cex.lab= 1.2, mgp=c(2.5,0.6,0), xlab="Year",
     ylab="Starting habitat cover (%)",col=ramp[(zval*10)-9])

# Add line for delayed compensation
abline(v=T+2020,lty=2,lwd=2)

# Create a vector of colours for use int he legend
leg.cols <- ramp[(c(1,3,5,7,9)*10)-9]

# Add legend
legend(2090,40,pch=15,col=leg.cols,legend=c("m = 1","m = 3","m = 5","m = 7", "m = 10"),cex=0.8, bg="white")

# Add points and arrows to denote the phase plane
points(x = (par("usr")[2]), y = 50, pch = 21, cex=1.5, col = "black",las = 1, xpd = TRUE,bg="white")
points(x = (par("usr")[2]), y = 100, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
points(x = (par("usr")[2]), y = 0, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 60,y1=75,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 40,y1=25,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)

# Label the panel
text(1995,105,"a",font=2,cex=1.8,las = 1, xpd = TRUE)

###############################################
#                                             #
#                 Panel B                     #
#                                             #
###############################################

# Creates a plot of the multplier value needed to achieve 'zero loss' by 2050 (i.e. 30 years after the start of simulation)
plot(min.m[,30],startH,pch=16,cex=0.7, las=1, cex.axis=1, cex.lab= 1.2, mgp=c(2.5,0.6,0), xlab="Multiplier",
     ylab="Starting habitat cover (%)",col="black",type="l", xlim=c(1,4),lwd=2)

# Add a line for multiplier needed to achieve 'zero loss' by 2035 (i.e. 15 years after the start of simulation)
lines(min.m[,15],startH,col="darkblue")
# Add a line for multiplier needed to achieve 'zero loss' by 2040 (i.e. 20 years after the start of simulation)
lines(min.m[,20],startH,col="darkblue",lty=3)
# Add a line for multiplier needed to achieve 'zero loss' by 2045 (i.e. 25 years after the start of simulation)
lines(min.m[,25],startH,col="blue")

# Add a legend to the plot
legend(3,30,col=c("darkblue","darkblue","blue","black"),legend=c("2035", "2040","2045","2050"),lty=c(1,3,1,1),lwd=c(1,1,1,2), bg="white", cex=0.8)

# Add points and arrows to denote the phase plane
points(x = (par("usr")[2]), y = 50, pch = 21, cex=1.5, col = "black",las = 1, xpd = TRUE,bg="white")
points(x = (par("usr")[2]), y = 100, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
points(x = (par("usr")[2]), y = 0, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 60,y1=75,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 40,y1=25,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)

# Label the panel
text(0.25,105,"b",font=2,cex=1.8,las = 1, xpd = TRUE)

# Close plotting device and save the figure to file
dev.off()
