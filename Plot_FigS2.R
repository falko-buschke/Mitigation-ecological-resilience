##################################################################################################

# The following script is to reproduce Figure S2 in the manuscript. 
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
a <- 1.5; i <- 0.01; T <- 10; d <- 0.03; m <- 1

times <- seq(0, 100, by = 0.1) # Define the time sequence for the simulation

startH <- seq(1,99,by=1)    # Define the Starting habitat cover

# Define the name of the plot (plus dimensions); alternative directories can be identified here
png(filename="FigS2.png",width=32,height=8,units="cm",res=300)

# Set plot margins and layouts
par(mai=c(0.5,0.6,0.15,0.1))
par(mfrow=c(1,4))



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
  
  dH <- (((log(H/(100-H)))/100) * (1/a) * H ) - (i * H) + (i * m * lag *(1+d)^T)
  list(dH, dy = dH,i)
}




###############################################
#                                             #
#               Panel A 	              #
#		(logit shape, a)	      #
#                                             #
###############################################

# Define the set of parameters 'a' to be used in the senstivity analysis
a.set <- c(0.5,1,1.5,2,2.5)

# Set the colours and line types for each curve
cols <- c("lightblue","blue", "blue", "darkblue","darkblue"); ltype <- c(1,2,1,2,1)


# Run a loop that incrementally adjusts the values for parameter 'a'
for (j in 1:length(a.set)) {

	# Set the value of parameter 'a' from the j-th value in the set of parameter values
	a <- a.set[j]

	# Create a blank vector to store the year at which 'zero loss' is achieved
	tnnl <- rep(NA,length(startH))

	# Run a loop to simulate dynamics for each starting value of H
	for (k in 1:length(startH)) {

		# Solve the equations numerically to get a time-series for the counterfactual (tcf) and mitigation (tmit) scenarios
		tcf <- dede(y = startH[k], times = times, func = eco_counter, parms = NULL, atol = 1e-10)
		tmit <- dede(y = startH[k], times = times, func = eco_mitigate, parms = NULL, atol = 1e-10)

		# Cap the maximum results to 100 in the counterfactual (i.e. maximum percentage vegetation cover) 
		CF <- tcf[,2] 
		CF[which(CF=="NaN")] <- 100
		CF[which(CF>100)] <- 100

		# Cap the maximum results to 100 in the mitigation scenario (i.e. maximum percentage vegetation cover) 
		Mit <- tmit[,2] 
		Mit[which(Mit=="NaN")] <- 100
		Mit[which(Mit>100)] <- 100

		# Record the year in whch 'zero loss' is achieved
		tnnl[k] <- times[which((CF-Mit)<=0)[2]]
	}

	# If 'zero loss' is not achieved, set the value to 100, the maximum time in this simulation
	tnnl[is.na(tnnl)] <- 100

	# Make the plot for panel (a). Use an if statement to make first plot, then add lines incrementally for each valeu of 'a'
	if (j ==1) {
		plot(tnnl+2020, startH, type="l",las=1,ylab="Starting habitat cover (%)", cex.axis=1.2, cex.lab= 1.5, mgp=c(2.5,0.6,0), xlab="Years", col=cols[j],lty=ltype[j], xlim=c(2020,2120))
	} else {lines(tnnl+2020,startH, col=cols[j], lty=ltype[j]) }

}

# Add line for delayed compensation
abline(v=T+2020,lty=2,lwd=2)

# Add legend
legend(2080,40, c("a = 0.5", "a = 1.0","a = 1.5","a = 2.0", "a = 2.5"), lty=ltype, col=cols)

# Add points and arrows to denote the phase plane
points(x = (par("usr")[2]), y = 50, pch = 21, cex=1.5, col = "black",las = 1, xpd = TRUE,bg="white")
points(x = (par("usr")[2]), y = 100, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
points(x = (par("usr")[2]), y = 0, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 60,y1=75,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 40,y1=25,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)

# Incorporate panel label 
text(1995,105,"a",font=2,cex=1.8,las = 1, xpd = TRUE)

# Reset the value of parameter 'a' to the default value
a <- 1.5


###############################################
#                                             #
#               Panel B                       #
#		(Development impact, i)	      #
#                                             #
###############################################

# Define the set of parameters 'i' to be used in the senstivity analysis
i.set <- c(0.005,0.01,0.015,0.02,0.025)

# Set the colours and line types for each curve
cols <- c("lightblue","blue", "blue", "darkblue","darkblue"); ltype <- c(1,2,1,2,1)

# Run a loop that incrementally adjusts the values for parameter 'i'
for (j in 1:length(i.set)) {

	# Set the value of parameter 'i' from the j-th value in the set of parameter values
	i <- i.set[j]

	# Create a blank vector to store the year at which 'zero loss' is achieved
	tnnl <- rep(NA,length(startH))

	# Run a loop to simulate dynamics for each starting value of H
	for (k in 1:length(startH)) {
		# Solve the equations numerically to get a time-series for the counterfactual (tcf) and mitigation (tmit) scenarios
		tcf <- dede(y = startH[k], times = times, func = eco_counter, parms = NULL, atol = 1e-10)
		tmit <- dede(y = startH[k], times = times, func = eco_mitigate, parms = NULL, atol = 1e-10)

		# Cap the maximum results to 100 for the counterfactual(i.e. maximum percentage vegetation cover) 
		CF <- tcf[,2] 
		CF[which(CF=="NaN")] <- 100
		CF[which(CF>100)] <- 100

		# Cap the maximum results to 100 for the mitigation scenario (i.e. maximum percentage vegetation cover) 
		Mit <- tmit[,2] 
		Mit[which(Mit=="NaN")] <- 100
		Mit[which(Mit>100)] <- 100

		# Record the year in whch 'zero loss' is achieved
		tnnl[k] <- times[which((CF-Mit)<=0)[2]]
	}

	# If 'zero loss' is not achieved, set the value to 100, the maximum time in this simulation
	tnnl[is.na(tnnl)] <- 100

	# Make the plot for panel (b). Use an if statement to make first plot, then add lines incrementally for each valeu of 'a'
	if (j ==1) {
		plot(tnnl+2020, startH, type="l",las=1,ylab="Starting habitat cover (%)", cex.axis=1.2, cex.lab= 1.5, mgp=c(2.5,0.6,0), xlab="Years", col=cols[j],lty=ltype[j], xlim=c(2020,2120))
	} else {lines(tnnl+2020,startH, col=cols[j], lty=ltype[j]) }
}

# Add line for delayed compensation
abline(v=T+2020,lty=2,lwd=2)

# Add legend
legend(2080,40, c("i = 0.005","i = 0.01","i = 0.015", "i = 0.02", "i = 0.25"), lty=ltype, col=cols)

# Add points and arrows to denote the phase plane
points(x = (par("usr")[2]), y = 50, pch = 21, cex=1.5, col = "black",las = 1, xpd = TRUE,bg="white")
points(x = (par("usr")[2]), y = 100, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
points(x = (par("usr")[2]), y = 0, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 60,y1=75,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 40,y1=25,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)

# Incorporate panel label 
text(1995,105,"b",font=2,cex=1.8,las = 1, xpd = TRUE)

# Reset the value of parameter 'a' to the default value
i <- 0.01


###############################################
#                                             #
#               Panel C                       #
#		(Discount rate, d)	      #
#                                             #
###############################################

# Define the set of parameters 'd' to be used in the senstivity analysis
d.set <- c(0.01,0.02,0.03,0.04,0.05)

# Set the colours and line types for each curve
cols <- c("lightblue","blue", "blue", "darkblue","darkblue"); ltype <- c(1,2,1,2,1)

# Run a loop that incrementally adjusts the values for parameter 'd'
for (j in 1:length(d.set)) {

	# Set the value of parameter 'd' from the j-th value in the set of parameter values
	d <- d.set[j]

	# Create a blank vector to store the year at which 'zero loss' is achieved
	tnnl <- rep(NA,length(startH))

	# Run a loop to simulate dynamics for each starting value of H
	for (k in 1:length(startH)) {
		# Solve the equations numerically to get a time-series for the counterfactual (tcf) and mitigation (tmit) scenarios
		tcf <- dede(y = startH[k], times = times, func = eco_counter, parms = NULL, atol = 1e-10)
		tmit <- dede(y = startH[k], times = times, func = eco_mitigate, parms = NULL, atol = 1e-10)

		# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
		CF <- tcf[,2] 
		CF[which(CF=="NaN")] <- 100
		CF[which(CF>100)] <- 100

		# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
		Mit <- tmit[,2] 
		Mit[which(Mit=="NaN")] <- 100
		Mit[which(Mit>100)] <- 100
	
		# Record the year in whch 'zero loss' is achieved
		tnnl[k] <- times[which((CF-Mit)<=0)[2]]
	}
	
	# If 'zero loss' is not achieved, set the value to 100, the maximum time in this simulation
	tnnl[is.na(tnnl)] <- 100

	# Make the plot for panel (c). Use an if statement to make first plot, then add lines incrementally for each valeu of 'a'
	if (j ==1) {
		plot(tnnl+2020, startH, type="l",las=1,ylab="Starting habitat cover (%)", cex.axis=1.2, cex.lab= 1.5, mgp=c(2.5,0.6,0), xlab="Years", col=cols[j],lty=ltype[j], xlim=c(2020,2120))
	} else {lines(tnnl+2020,startH, col=cols[j], lty=ltype[j]) }
}

# Add line for delayed compensation
abline(v=T+2020,lty=2,lwd=2)

# Add legend
legend(2085,100, c("d = 0.01", "d = 0.02","d = 0.03","d = 0.04", "d = 0.05"), lty=ltype, col=cols)

# Add points and arrows to denote the phase plane
points(x = (par("usr")[2]), y = 50, pch = 21, cex=1.5, col = "black",las = 1, xpd = TRUE,bg="white")
points(x = (par("usr")[2]), y = 100, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
points(x = (par("usr")[2]), y = 0, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 60,y1=75,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 40,y1=25,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)

# Incorporate panel label 
text(1995,105,"c",font=2,cex=1.8,las = 1, xpd = TRUE)

d <- 0.03

###############################################
#                                             #
#               Panel D                       #
#		(Compensation delay, T)	      #
#                                             #
###############################################

# Define the set of parameters 'T' to be used in the senstivity analysis
T.set <- c(2,6,10,14,18)

# Set the colours and line types for each curve
cols <- c("lightblue","blue", "blue", "darkblue","darkblue"); ltype <- c(1,2,1,2,1)

# Run a loop that incrementally adjusts the values for parameter 'T'
for (j in 1:length(T.set)) {

	# Set the value of parameter 'T' from the j-th value in the set of parameter values
	T <- T.set[j]

	# Create a blank vector to store the year at which 'zero loss' is achieved
	tnnl <- rep(NA,length(startH))
	
	# Run a loop to simulate dynamics for each starting value of H
	for (k in 1:length(startH)) {
		# Solve the equations numerically to get a time-series for the counterfactual (tcf) and mitigation (tmit) scenarios
		tcf <- dede(y = startH[k], times = times, func = eco_counter, parms = NULL, atol = 1e-10)
		tmit <- dede(y = startH[k], times = times, func = eco_mitigate, parms = NULL, atol = 1e-10)

		# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
		CF <- tcf[,2] 
		CF[which(CF=="NaN")] <- 100
		CF[which(CF>100)] <- 100

		# Cap the maximum results to 100 (i.e. maximum percentage vegetation cover) 
		Mit <- tmit[,2] 
		Mit[which(Mit=="NaN")] <- 100
		Mit[which(Mit>100)] <- 100
		
		# Record the year in whch 'zero loss' is achieved
		tnnl[k] <- times[which((CF-Mit)<=0)[2]]

	}
	
	# If 'zero loss' is not achieved, set the value to 100, the maximum time in this simulation
	tnnl[is.na(tnnl)] <- 100

	# Make the plot for panel (d). Use an if statement to make first plot, then add lines incrementally for each valeu of 'a'
	if (j ==1) {
		plot(tnnl+2020, startH, type="l",las=1,ylab="Starting habitat cover (%)", cex.axis=1.2, cex.lab= 1.5, mgp=c(2.5,0.6,0), xlab="Years", col=cols[j],lty=ltype[j], xlim=c(2020,2120))
	} else {lines(tnnl+2020,startH, col=cols[j], lty=ltype[j]) }
}

# Add legend
legend(2080,40, c("T = 2", "T = 6","T = 10","T = 14", "T = 18"), lty=ltype, col=cols)

# Add points and arrows to denote the phase plane
points(x = (par("usr")[2]), y = 50, pch = 21, cex=1.5, col = "black",las = 1, xpd = TRUE,bg="white")
points(x = (par("usr")[2]), y = 100, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
points(x = (par("usr")[2]), y = 0, pch = 16, cex=1.5, col = "black",las = 1, xpd = TRUE)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 60,y1=75,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)
Arrows(x0 = (par("usr")[2]), x1 = (par("usr")[2]),y0 = 40,y1=25,las = 1, xpd = TRUE,arr.type="triangle", arr.width=0.15)

# Incorporate panel label 
text(1995,105,"d",font=2,cex=1.8,las = 1, xpd = TRUE)

# Close plotting device and save figure to file
dev.off()
