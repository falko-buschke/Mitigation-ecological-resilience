##################################################################################################

# The following script is to reproduce Figure S1 in the supplementary material. 
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
r <- 0.03 ; K <-  100; A <- 25; i <- 0.01; T <- 10; d <- 0.03

# Define the name of the plot (plus dimensions); alrernative directories can be identified here
png(filename="FigS1.png",width=32,height=8,units="cm",res=300)

# Set plot margins and layouts
par(mai=c(0.5,0.6,0.15,0.1))
par(mfrow=c(1,4))

###############################################
#                                             #
#                 Panel A                     #
#                                             #
###############################################

P <- 0:100      # Create a vector of population size, from 0 to 100

# Simulate the expected population response (dP/dt) for each level of P for the counterfactual and impact scenarios
P_counter <- (((log(P/(100-P)))/100) * (1/a))*P
P_impact <- ((((log(P/(100-P)))/100) * (1/a))*P) - (i*P)

P_counter <- (r * P * ((P/A)-1) * (1 - P/K)) 
P_impact <- (r * P * ((P/A)-1) * (1 - P/K)) - (i*P)

# Make base plots
plot(P~P_counter,type="l",las=1, cex.axis=1.2, cex.lab= 1.5, mgp=c(2.5,0.6,0), xlab= expression(Population~response~(dP/dt)),
     ylab="Population size",col="black",ylim=c(0,100),xlim=c(-1,3))
lines(P~P_impact,col="blue",lty=1)
# Add grey line to denote zero response (i.e. equilibria)
abline(v=0,lwd=2,col="grey")

# Add dashed lines that link the equilibria points to the phase plane denoted on the secondary vertical axis
lines(c(0,(par("usr")[2])),c(P[which(P_impact>=0)[2]],P[which(P_impact>=0)[2]]),lty=3,col="blue")
lines(c(0,(par("usr")[2])),c(A,A),lty=3)
lines(c(0,(par("usr")[2])),c(86.5,86.5),lty=3,col="blue")

# Add arrows to the phase plane that denote population responses relative to equilbrium points (requires 'shape' package)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 50,y1=60,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 20,y1=12,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 100,y1=93,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15, col="darkblue")
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 38,y1=32,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15, col="darkblue")

# Add solid and open points to denote stable and unstable equilibria, respectively, to the pahse plane 

points(x = (par("usr")[2]), y = 25, pch = 21, cex=1.5, col = "black",las = 1, xpd = TRUE,bg="white")
points(x = (par("usr")[2]), y = 100, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
points(x = (par("usr")[2]), y = 0, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
points(x = (par("usr")[2]), y = P[which(P_impact>=0)[2]], pch = 21, cex=1.5, col = "blue",bg="white",las = 1, xpd = TRUE)
points(x = (par("usr")[2]), y = 86.5, pch = 16, cex=1.5, col = "darkblue",las = 1, xpd = TRUE,bg="white")

# Add a legend to the plot
legend(0.3,20,lty=1,legend=c("Impact (i) = 0","Impact (i) = 0.01"),col=c("black","blue"),cex=1)

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

# The counterfactual is only defined according to the dynamics of a logistic function with Allee effect
pop_counter <- function(t, P, parms) {
  dP <- (r * P * ((P/A)-1) * (1 - P/K))  
  list(dP, dy = dP,i)
}

# The mitigation scenario is include delayed (lagged) compensation as soon as the time of the simulation is longer than parameter T 
pop_mitigate <- function(t, P, parms) {
  if (t < T)
    lag <- 0
  else
    lag <- lagvalue(t - T)
  
  dP <- ((r * P * ((P/A)-1) * (1 - P/K))) - (i * P) + (i * lag *(1+d)^T)
  list(dP, dy = dP,i)
}

########################
########################

# Define the time series and increment ('by' defines the time increment, lower values result in smoother curves)
times <- seq(0, 100, by = .1)

startP <- 20    # Define the Starting population

# Solve the equations numerically to get a time-series for the counterfactual (tcf) and mitigation (tmit) scenarios
tcf <- dede(y = startP, times = times, func = pop_counter, parms = NULL, atol = 1e-10)
tmit <- dede(y = startP, times = times, func = pop_mitigate, parms = NULL, atol = 1e-10)

# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
CF <- tcf[,2] 
CF[which(CF=="NaN")] <- 100
CF[which(CF>100)] <- 100

# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
Mit <- tmit[,2] 
Mit[which(Mit=="NaN")] <- 100
Mit[which(Mit>100)] <- 100

# Plot the trajecotries for the scenarios at starting population
plot(times+2020,CF,type="l",las=1, cex.axis=1.2, cex.lab= 1.5, mgp=c(2.5,0.6,0), xlab="Year",
     ylab="Population size",col="red",ylim=c(0,100))
lines(times+2020,Mit,lty=2,col="red")
# Add point where curves intersect
points(times[which((CF-Mit)<=0)[2]]+2020,CF[which((CF-Mit)<=0)[2]], cex=1.5,pch=18,col="red")

######################################################

# Repeat sequence for different starting populations

######################################################

startP <- 40    # Define the Starting population

# Solve the equations numerically to get a time-series for the counterfactual (tcf) and mitigation (tmit) scenarios
tcf <- dede(y = startP, times = times, func = pop_counter, parms = NULL, atol = 1e-10)
tmit <- dede(y = startP, times = times, func = pop_mitigate, parms = NULL, atol = 1e-10)

# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
CF <- tcf[,2] 
CF[which(CF=="NaN")] <- 100
CF[which(CF>100)] <- 100

# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
Mit <- tmit[,2] 
Mit[which(Mit=="NaN")] <- 100
Mit[which(Mit>100)] <- 100

# Plot the trajectories for the scenarios at starting population
lines(times+2020,CF,lty=1,col="orange")
lines(times+2020,Mit,lty=2,col="orange")
# Add point where curves intersect
points(times[which((CF-Mit)<=0)[2]]+2020,CF[which((CF-Mit)<=0)[2]], cex=1.5,pch=18,col="orange")


######################################################

# Repeat sequence for different starting population

######################################################

startP <- 60    # Define the Starting population

# Solve the equations numerically to get a time-series for the counterfactual (tcf) and mitigation (tmit) scenarios
tcf <- dede(y = startP, times = times, func = pop_counter, parms = NULL, atol = 1e-10)
tmit <- dede(y = startP, times = times, func = pop_mitigate, parms = NULL, atol = 1e-10)

# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
CF <- tcf[,2] 
CF[which(CF=="NaN")] <- 100
CF[which(CF>100)] <- 100

# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
Mit <- tmit[,2] 
Mit[which(Mit=="NaN")] <- 100
Mit[which(Mit>100)] <- 100

# Plot the trajectories for the scenarios at starting population
lines(times+2020,CF,lty=1,col="gold")
lines(times+2020,Mit,lty=2,col="gold")
# Add point where curves intersect
points(times[which((CF-Mit)<=0)[2]]+2020,CF[which((CF-Mit)<=0)[2]], cex=1.5,pch=18,col="gold")


######################################################

# Repeat sequence for different starting population

######################################################

startP <- 80    # Define the Starting population

# Solve the equations numerically to get a time-series for the counterfactual (tcf) and mitigation (tmit) scenarios
tcf <- dede(y = startP, times = times, func = pop_counter, parms = NULL, atol = 1e-10)
tmit <- dede(y = startP, times = times, func = pop_mitigate, parms = NULL, atol = 1e-10)

# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
CF <- tcf[,2] 
CF[which(CF=="NaN")] <- 100
CF[which(CF>100)] <- 100

# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
Mit <- tmit[,2] 
Mit[which(Mit=="NaN")] <- 100
Mit[which(Mit>100)] <- 100

# Plot the trajectories for the scenarios at starting population
lines(times+2020,CF,lty=1,col="darkgreen")
lines(times+2020,Mit,lty=2,col="darkgreen")
# Add point where curves intersect
points(times[which((CF-Mit)<=0)[2]]+2020,CF[which((CF-Mit)<=0)[2]], cex=1.5,pch=18,col="darkgreen")

##################################################

# Add dashed line to indicate mitigation delay
abline(v=T+2020,lty=2,lwd=2)

# Add points and arrows to denote the phase plane
points(x = (par("usr")[2]), y = A, pch = 21, cex=1.5, col = "black",las = 1, xpd = TRUE,bg="white")
points(x = (par("usr")[2]), y = K, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
points(x = (par("usr")[2]), y = 0, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 50,y1=60,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 20,y1=12,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)

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

# Define all the starting values for population size
# Note: These are integers (excluding stable equilibria) in order to define the plot surface
startP <- seq(1,99,by=1)

# Define the time sequence
# Note: These are integers in order to define the plot surface
times <- seq(0, 100, by = 1)

# Creat a blank output matrix that records the mitigation debt
output <- matrix(NA,ncol=length(times),nrow=length(startP))

# Run a loop that solves the differential equations for every starting population
for (k in 1:length(startP)) {

# Solve equations for counterfactual and mitigation scenarios  
  tcf <- dede(y = startP[k], times = times, func = pop_counter, parms = NULL, atol = 1e-10)
  tmit <- dede(y = startP[k], times = times, func = pop_mitigate, parms = NULL, atol = 1e-10)
  
 # Cap maximum population 100 
  CF <- tcf[,2] 
  CF[which(CF=="NaN")] <- 100
  CF[which(CF>100)] <- 100
  
  Mit <- tmit[,2] 
  Mit[which(Mit=="NaN")] <- 100
  Mit[which(Mit>100)] <- 100

# Record the mitgation debt for the time series to the blanck output matrix
  output[k,] <- CF - Mit 
}

# The panel is actually a scatterplot, so define the x- and y-values for each combination of time and starting population
yval <- rep(1:99,rep(101,99))
xval <- rep(2020:2120,99)

# Convert the output matrix to a vextor (z-values) that will define the colout of mitigation debt
zval <- as.vector(t(output))

# Define a colour ramp 
ramp <- colorRampPalette(c("darkblue","blue","white","red","darkred"),interpolate="linear")(21)


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
  } else {nnl[l] <- NA}
}


# This makes the plot for each combination ot time and starting population, colouring the points according to the mitigation debt
plot(xval,yval,pch=15,cex=0.48, las=1, cex.axis=1.2, cex.lab= 1.5, mgp=c(2.5,0.6,0), xlab="Year",
     ylab="Starting population",col=ramp[ceiling(zval+11)])

# Add the grey line to denote NNL (i.e. the transition from net loss to net gain)
lines(nnl+2020,c(1:99),lty=3,lwd=2,col="darkgrey")

# Add line for delayed compensation
abline(v=T+2020,lty=2,lwd=2)

# Add legend
legend(2095,100,pch=15,col=c("darkblue","blue","grey","red","darkred"),legend=c("+10","+5","0","-5","-10"),cex=1)

# Add points and arrows to denote the phase plane
points(x = (par("usr")[2]), y = A, pch = 21, cex=1.5, col = "black",las = 1, xpd = TRUE,bg="white")
points(x = (par("usr")[2]), y = K, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
points(x = (par("usr")[2]), y = 0, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 50,y1=60,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 20,y1=12,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)

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
space<- read.table("NNLSpaceSpecies.txt",header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

# Make the plot, beginning with starting population of 70
plot(space$Delay_70,space$Impact_70,type="l",ylim=c(0,.15),xlim=c(0,30),xlab="Mitigation delay (T)",
     cex.axis=1.2, cex.lab= 1.5, mgp=c(2.5,0.6,0),las=1,ylab="Impact (i)",col="lightblue")

# Add a polygon to denote Net los by 2050 (added first before the other lines, which will be ploted above this polygon)
polygon(c(space$Delay_95,30,30),c(space$Impact_95,0.15,0),col="lightgrey",border=F)

# Sequentially add lines for different starting population
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
  legend=c("95 (+0.4%)","90 (+0.8%)","85 (+1.1%)","80 (+1.3%)","75 (+1.5%)","70 (+1.6%)"),cex=1)

# Label the panel
text(-7.5,0.1585,"d",font=2,cex=1.8,las = 1, xpd = TRUE)

# Close the plotting device to save plot to file
dev.off()


