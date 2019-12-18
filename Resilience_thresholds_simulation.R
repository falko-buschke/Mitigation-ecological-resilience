##################################################################################

# The following code tests sequentially whether NNL is attainbale by 2050.
# It explores combinations of impact magnitude (i) and mitigation delay (T).
# User must manually select the statrting biodiversity (habitat cover, H, or population size, P),
# ...and the simulation will save to file the combinations of i and T that allow NNL by 2050.

# Please not that this is quite a time intensive simulation, so to test the code:
#			- Use smaller increments for T and or i
#			- Use shorter time increments	

##################################################################################

# First, begin by installing the require packages.

# 'deSolve' is required to solve the differential equation
install.packages(c("deSolve", "shape"), dependencies = TRUE)

# load packages
library(deSolve)


# Define all the model parameters

r <- 0.03 ; K <-  100; A <- 25; a <- 1.5; d <- 0.03

#Define the starting value for habitat cover and population size.
Start_value <- 70


# It is necessary to define the derivative equations that need to be solved
#-----------------------------
# The derivative functions
#-----------------------------

#-----------------------------
# Ecosystems
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

#-----------------------------
# For population size
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


# This is a free variable to select whether the simulation is for habitat cover or population size.
# For habitat, make:
cf_func <- eco_counter
imp_func <- eco_mitigate

# For population, make:
#cf_func <- pop_counter
#imp_func <- pop_mitigate

####################################################################################
####################################################################################

# Create a vector for all the values of T (i.e. the delay vector)
dl.vec <- seq(0.1,30,by=0.1)

# Create a vector for all the values of i (i.e. the impact vector)
imp <- seq(0.001,0.15,by=0.001)

# Create a blank matrx that will be used to record the time until NNL for each cobination of T and i
NNL.t <- matrix(NA,nrow=length(imp),ncol=length(dl.vec))


# Use a nested loop to make calculation for each combination of i and T. 

for (k in 1:length(dl.vec)){
	for (j in 1:length(imp)){
		T <- dl.vec[k]
		i <- imp[j]

# The following section of code was added to speed up simulations by not calculating NNL for
# combination of large T and large i (becasue NNL will not be achevied int he time frame)
		if(T> 10 & i > 0.1) {
			NNL.t[j,k] <- 0
		} else {

#-----------------------------
# Time sequence
#-----------------------------
# The time sequence use was for the next 50 years (2050 is 30 years from now).
# The simualtion can be shortened by reducing the time-series, or widening the time interval.
times <- seq(0, 50, by = 0.1)

# Selects the start value as defined in line 28 above
startB <- Start_value

# Solve the differntial equations
tcf <- dede(y = startB, times = times, func = cf_func, parms = NULL, atol = 1e-10)
toff <- dede(y = startB, times = times, func = imp_func, parms = NULL, atol = 1e-10)

# Cap maximum values to 100
Nat <- tcf[,2] 
Nat[which(Nat=="NaN")] <- 100

Off <- toff[,2] 
Off[which(Off=="NaN")] <- 100

# Calculate the mitigation debt
NNL <- Off - Nat

# Report if there is no net loss or net gain by 2050 (i.e. 30 years in the future)
NNL.t[j,k] <- ifelse(NNL[which(times==30)]>=0,1,0)
}


# The following function is just a loop counter to keep track of the simulation. 
} 
        if (k/10 == floor(k/10))
        	{print(k)
        	flush.console()}
}


##########################################################################
##########################################################################

# The following loop is used to identify the longest value of T that corresponds to the highest value of i
# ... but which still leads to NNL outcomes by 2050

# Create blank vector to record maxumim values for delay T.
dl.val <- rep(NA,length(imp))

for (n in 1:length(imp)){
	if (sum(NNL.t[n,])>0){
		dl.val[n] <- max(which(NNL.t[n,]==1))
	} else {
		dl.val[n] <- 1	
	}
}

# Make a plot of the combinations of T and i that allow NNL by 2050
plot(dl.vec[dl.val],imp,type="l",ylim=c(0,.15),xlim=c(0,30),xlab="Mitigation delay (T)",las=1,ylab="Impact (i)",col="red")

# Write outputs to file
write.table(cbind(dl.vec[dl.val],imp),file= "OffsetSpace.txt",
	quote=TRUE,sep="\t",row.names=FALSE,col.names=c("Delay","Impact"))


