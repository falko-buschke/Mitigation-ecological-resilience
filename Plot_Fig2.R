##################################################################################################

# The following script is to reproduce Figure 2 in the manuscript. 
# The four panels are plotted in sequence, and saved in the default working directory.

##################################################################################################

# First, begin by installing the require packages.

# 'deSolve' is required to solve the differential equation
# 'shape' is required to include the specialised arros that denote the phase plane.
install.packages(c("deSolve", "shape"), dependencies = TRUE)

# load packages
library(deSolve)
library(shape)


# Set all the default parameters for the model
a <- 1.5; i <- 0.01; T <- 10; d <- 0.03

# Define the name of the plot (plus dimensions); alrernative directories can be identified here
png(filename="Fig2.png",width=32,height=8,units="cm",res=300)

# Set plot margins and layouts
par(mai=c(0.5,0.6,0.15,0.1))
par(mfrow=c(1,4))

###############################################
#                                             #
#                 Panel A                     #
#                                             #
###############################################

H <- 0:100      # Create a vector of habitat cover, from 0% to 100%

# Simulate the expected habitat response (dH/dt) for each level of H for the counterfactual and impact scenarios
H_counter <- (((log(H/(100-H)))/100) * (1/a))*H
H_impact <- ((((log(H/(100-H)))/100) * (1/a))*H) - (i*H)

# Make base plots
plot(H~H_counter,type="l",las=1, cex.axis=1.2, cex.lab= 1.5, mgp=c(2.5,0.6,0), xlab= expression(Habitat~response~(dH/dt)),
     ylab="Habitat cover (%)",col="black",ylim=c(0,100),xlim=c(-1,3))
lines(H~H_impact,col="blue",lty=1)
# Add grey line to denote zero response (i.e. equilibria)
abline(v=0,lwd=2,col="grey")

# Add dashed lines that link the equilibria pointss to the phase plane denoted on the secondary vertical axis
lines(c(0,3,c(H[which(H_impact>=0)[1]],H[which(H_impact>=0)[1]]),lty=3,col="blue")
lines(c(0,(par("usr")[2])),c(50,50),lty=3)

# Add vertical line for second phase plane in the counerfactual
abline(v=3,col="darkblue")

# Add arrows to the phase plane that denote habitat responses relative to equilbrium points (requires 'shape' package)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 60,y1=75,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 40,y1=25,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)
Arrows(x0 = 3, x1 = 3,y0 = 80,y1=H[which(H_impact>=0)[1]]+10,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15, col="darkblue")
Arrows(x0 = 3, x1 = 3,y0 = 75,y1=65,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15, col="darkblue")

# Add solid and open points to denote stable and unstable equilibria, respectively, to the phase plane 
points(x = (par("usr")[2]), y = 50, pch = 21, cex=1.5, col = "black",las = 1, xpd = TRUE,bg="white")
points(x = (par("usr")[2]), y = 100, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
points(x = (par("usr")[2]), y = 0, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
points(x = 3, y = H[which(H_impact>=0)[1]], pch = 21, cex=1.5, col = "blue",bg="white",las = 1, xpd = TRUE)
points(x = 3, y = 100, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
points(x = 3, y = 0, pch = 16, cex=1.5, col = "darkblue",las = 1, xpd = TRUE)

# Add a legend to the plot
legend(0.75,20,lty=1,legend=c("Impact (i) = 0","Impact (i) = 0.01"),col=c("black","blue"),cex=1)


# Add a legend to the plot
legend(0.75,20,lty=1,legend=c("Impact (i) = 0","Impact (i) = 0.01"),col=c("black","blue"),cex=1)

# Add the panel label
text(-2,105,"a",font=2,cex=1.8,las = 1, xpd = TRUE)



###############################################
#                                             #
#                 Panel B                     #
#                                             #
###############################################

# To begin this panel, it is first necessary to define all the derivative equations that need to be solved

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
  
  dH <- (((log(H/(100-H)))/100) * (1/a) * H ) - (i * H) + (i * lag *(1+d)^T)
  list(dH, dy = dH,i)
}


########################
########################

# Define the time series and increment ('by' defines the time increment, lower values result in smoother curves)
times <- seq(0, 100, by = .1)

startH <- 20    # Define the Starting habitat cover

# Solve the equations numerically to get a time-series for the counterfactual (tcf) and mitigation (tmit) scenarios
tcf <- dede(y = startH, times = times, func = eco_counter, parms = NULL, atol = 1e-10)
tmit <- dede(y = startH, times = times, func = eco_mitigate, parms = NULL, atol = 1e-10)

# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
CF <- tcf[,2] 
CF[which(CF=="NaN")] <- 100
CF[which(CF>100)] <- 100

# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
Mit <- tmit[,2] 
Mit[which(Mit=="NaN")] <- 100
Mit[which(Mit>100)] <- 100

# Plot the trajecotries for the scenarios at starting habitat cover of 20%
plot(times+2020,CF,type="l",las=1, cex.axis=1.2, cex.lab= 1.5, mgp=c(2.5,0.6,0), xlab="Year",
     ylab="Habitat cover (%)",col="red",ylim=c(0,100))
lines(times+2020,Mit,lty=2,col="red")
# Add point where curves intersect
points(times[which((CF-Mit)<=0)[2]]+2020,CF[which((CF-Mit)<=0)[2]], cex=1.5,pch=18,col="red")

######################################################

# Repeat sequence for different starting habitat

######################################################

startH <- 40    # Define the Starting habitat cover

# Solve the equations numerically to get a time-series for the counterfactual (tcf) and mitigation (tmit) scenarios
tcf <- dede(y = startH, times = times, func = eco_counter, parms = NULL, atol = 1e-10)
tmit <- dede(y = startH, times = times, func = eco_mitigate, parms = NULL, atol = 1e-10)

# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
CF <- tcf[,2] 
CF[which(CF=="NaN")] <- 100
CF[which(CF>100)] <- 100

# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
Mit <- tmit[,2] 
Mit[which(Mit=="NaN")] <- 100
Mit[which(Mit>100)] <- 100

# Plot the trajectories for the scenarios at starting habitat cover of 40%
lines(times+2020,CF,lty=1,col="orange")
lines(times+2020,Mit,lty=2,col="orange")
# Add point where curves intersect
points(times[which((CF-Mit)<=0)[2]]+2020,CF[which((CF-Mit)<=0)[2]], cex=1.5,pch=18,col="orange")


######################################################

# Repeat sequence for different starting habitat

######################################################

startH <- 60    # Define the Starting habitat cover

# Solve the equations numerically to get a time-series for the counterfactual (tcf) and mitigation (tmit) scenarios
tcf <- dede(y = startH, times = times, func = eco_counter, parms = NULL, atol = 1e-10)
tmit <- dede(y = startH, times = times, func = eco_mitigate, parms = NULL, atol = 1e-10)

# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
CF <- tcf[,2] 
CF[which(CF=="NaN")] <- 100
CF[which(CF>100)] <- 100

# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
Mit <- tmit[,2] 
Mit[which(Mit=="NaN")] <- 100
Mit[which(Mit>100)] <- 100

# Plot the trajectories for the scenarios at starting habitat cover of 60%
lines(times+2020,CF,lty=1,col="gold")
lines(times+2020,Mit,lty=2,col="gold")
# Add point where curves intersect
points(times[which((CF-Mit)<=0)[2]]+2020,CF[which((CF-Mit)<=0)[2]], cex=1.5,pch=18,col="gold")


######################################################

# Repeat sequence for different starting habitat

######################################################

startH <- 80    # Define the Starting habitat cover

# Solve the equations numerically to get a time-series for the counterfactual (tcf) and mitigation (tmit) scenarios
tcf <- dede(y = startH, times = times, func = eco_counter, parms = NULL, atol = 1e-10)
tmit <- dede(y = startH, times = times, func = eco_mitigate, parms = NULL, atol = 1e-10)

# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
CF <- tcf[,2] 
CF[which(CF=="NaN")] <- 100
CF[which(CF>100)] <- 100

# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
Mit <- tmit[,2] 
Mit[which(Mit=="NaN")] <- 100
Mit[which(Mit>100)] <- 100

# Plot the trajectories for the scenarios at starting habitat cover of 80%
lines(times+2020,CF,lty=1,col="darkgreen")
lines(times+2020,Mit,lty=2,col="darkgreen")
# Add point where curves intersect
points(times[which((CF-Mit)<=0)[2]]+2020,CF[which((CF-Mit)<=0)[2]], cex=1.5,pch=18,col="darkgreen")

##################################################

# Add dashed line to indicate mitigation delay
abline(v=T+2020,lty=2,lwd=2)

# Add points and arrows to denote the phase plane
points(x = (par("usr")[2]), y = 50, pch = 21, cex=1.5, col = "black",las = 1, xpd = TRUE,bg="white")
points(x = (par("usr")[2]), y = 100, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
points(x = (par("usr")[2]), y = 0, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)

Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 60,y1=75,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 40,y1=25,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)


# Add a legend
legend(2070,70,lty=c(1,2),legend=c("Counterfactual","Mitigation"),cex=1)

# Label the plot panel
text(1995,105,"b",font=2,cex=1.8,las = 1, xpd = TRUE)


###############################################
#                                             #
#                 Panel C                     #
#                                             #
###############################################


# This panel inherits the differntial equations from above

# Define all the starting values for habitat cover
# Note: These are integers (excluding stable equilibria) in order to define the plot surface
startH <- seq(1,99,by=1)

# Define the time sequence
# Note: These are integers in order to define the plot surface
times <- seq(0, 100, by = 1)

# Creat a blank output matrix that records the mitigation debt
output <- matrix(NA,ncol=length(times),nrow=length(startH))

# Run a loop that solves the differential equations for every starting level of habitat cover
for (k in 1:length(startH)) {

# Solve equations for counterfactual and mitigation scenarios  
  tcf <- dede(y = startH[k], times = times, func = eco_counter, parms = NULL, atol = 1e-10)
  tmit <- dede(y = startH[k], times = times, func = eco_mitigate, parms = NULL, atol = 1e-10)
  
 # Cap maximum habitat cover to 100% 
  CF <- tcf[,2] 
  CF[which(CF=="NaN")] <- 100
  CF[which(CF>100)] <- 100
  
  Mit <- tmit[,2] 
  Mit[which(Mit=="NaN")] <- 100
  Mit[which(Mit>100)] <- 100

# Record the mitgation debt for the time series to the blanck output matrix
  output[k,] <- CF - Mit 
}

# The panel is actually a scatterplot, so define the x- and y-values for each combination of time and starting habitat
yval <- rep(1:99,rep(101,99))
xval <- rep(2020:2120,99)

# Convert the output matrix to a vextor (z-values) that will define the colout of mitigation debt
zval <- as.vector(t(output))

# Define a colour ramp 
ramp <- colorRampPalette(c("darkblue","blue","white","red","darkred"),interpolate="linear")(41)


# The following function identifies the time where NNL is achieved
# While seemingly simple, it is based on the trasition from net loss to net gain (becasue exact NNL, where the difference between scenarios is zero,
# ... is quite rare). Therefore, it first replaced perfect NNL with a small net gain, then multiplies sequential
# ... time value. Thus, the time where the product is negative denotes the transition from net loss to net gain.
nnl <- rep(NA,99)
for (l in 1:99){
  vr1 <- output[l,]
  vr1[which(vr1==0)] <- -0.00000001
  vr1[1] <- 0
  vr2 <- rep(NA,99)
  
  for (j in 1:99) {
    vr2[j] <- vr1[j]*vr1[j+1]
  }
  if (min(vr2)<0) {
    nnl[l]<-which(vr2<0)
  } else {nnl[l] <- 0}
}

# This makes the plot for each combination ot time and starting habitat, colouring the points according to the mitigation debt
plot(xval,yval,pch=15,cex=0.48, las=1, cex.axis=1.2, cex.lab= 1.5, mgp=c(2.5,0.6,0), xlab="Year",
     ylab="Starting habitat cover (%)",col=ramp[ceiling(zval+21)])

# Add the grey line to denote NNL (i.e. the transition from net loss to net gain)
lines(nnl+2020,c(1:99),lty=3,lwd=2,col="darkgrey")

# Add line for delayed compensation
abline(v=T+2020,lty=2,lwd=2)

# Add legend
legend(2095,100,pch=15,col=c("darkblue","blue","grey","red","darkred"),legend=c("+20","+10","0","-10","-20"),cex=1)

# Add points and arrows to denote the phase plane
points(x = (par("usr")[2]), y = 50, pch = 21, cex=1.5, col = "black",las = 1, xpd = TRUE,bg="white")
points(x = (par("usr")[2]), y = 100, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
points(x = (par("usr")[2]), y = 0, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 60,y1=75,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 40,y1=25,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)

# Incorporate panel label 
text(1995,105,"c",font=2,cex=1.8,las = 1, xpd = TRUE)



###############################################
#                                             #
#                 Panel D                    #
#                                             #
###############################################

# This panel requires the simulated data used to identify the mitigation space
# ... (i.e. the combination of Mitigation delay (T) and Impact (i) that allow for NNL)

# Read the data from file
space<- read.table("NNLSpaceHabitat.txt",header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

# Make the plot, beginning with starting habitat cover of 70%
plot(space$Delay_70,space$Impact_70,type="l",ylim=c(0,.15),xlim=c(0,30),xlab="Mitigation delay (T)",
     cex.axis=1.2, cex.lab= 1.5, mgp=c(2.5,0.6,0),las=1,ylab="Impact (i)",col="lightblue")

# Add a polygon to denote Net los by 2050 (added first before the other lines, which will be ploted above this polygon)
polygon(c(space$Delay_95,30),c(space$Impact_95,0.15),col="lightgrey",border=F)

# Sequentially add lines for different starting habitat cover
lines(space$Delay_75,space$Impact_75,col="lightblue",lty=2)
lines(space$Delay_80,space$Impact_80,col="blue")
lines(space$Delay_85,space$Impact_85,col="blue",lty=2)
lines(space$Delay_90,space$Impact_90,col="darkblue")
lines(space$Delay_95,space$Impact_95,col="darkblue",lty=2)

# Add a label to show that shaded polygon implies net loss
text(19,0.07,"Net loss by 2050",cex=1.5)

# Add a legend
legend(16,0.147,lty=c(2,1), 
	col=c("darkblue", "darkblue","blue","blue","lightblue","lightblue"),
	legend=c("95 (+2%)","90 (+1.5%)","85 (+1.2%)","80 (+0.9%)","75 (+0.7%)","70 (+0.6%)"),cex=1)

# Label the panel
text(-7.5,0.1585,"d",font=2,cex=1.8,las = 1, xpd = TRUE)

# Close the plotting device to save plot to file
dev.off()


