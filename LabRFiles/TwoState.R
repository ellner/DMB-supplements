# Dynamic Models in Biology, Stephen Ellner and John Guckenheimer
# TwoState.R simulates a group of particles all moving 
# back and forth between two compartments, using a coin-tossing
# process. In time interval dt, each particle in the left box has
# probability R*dt of moving to the right box, and each one in the right
# box has probability L*dt of moving into the left box. For dt small
# and a large number of particles, the process should approximate
# a two-compartment differential equation
#
#   dQ1/dt= -R*Q1 + L*Q2      dQ2/dt= -L*Q2 + R*Q11
#
# Simulations of the Markov Chain are plotted along with solutions
# of the differential equation. 

rm(list=ls(all=TRUE)); # clear past actions from R's memory 
require(deSolve); 

# parameters
R=0.3; L=0.2; dt=0.1; tmax=25; 

# initial condition
nL0=250; nR0=0; ntot=nL0+nR0; 

# movement probabilities per time step 
Rdt=R*dt; Ldt=L*dt;

# initialize variables for current population size
nL=nL0; nR=nR0; 

# initialize vectors storing the results over time
ntL=nL; ntR=nR;  t=0; tvals=0;

while(t<tmax) {
   # how many move from L to R?
    moveR=ifelse(nL>0, sum(runif(nL)<Rdt), 0);
   # how many move from R to L?
    moveL=ifelse(nR>0, sum(runif(nR)<Ldt), 0);
   # make the changes
     nL=nL + moveL - moveR;
     nR=nR + moveR - moveL;
   # store the results for plotting
   	t=t+dt; tvals=c(tvals,t);
	ntL=c(ntL,nL); ntR=c(ntR,nR);
} 

# solve the ODEs for the continuous time chain (dt -> 0)
diffus2=function(t,y,parms) {
  dy= c( -R*y[1]+L*y[2], R*y[1]-L*y[2] );
  return(list(dy))
}	 

out=lsoda(y=c(nL0,nR0),times=seq(0,tmax,by=dt),func=diffus2,parms=0);

# plot stochastic simulation and deterministic solutions
par(cex.axis=1.25,cex.lab=1.25); 
matplot(tvals,cbind(ntL,ntR),type="l",lty=1,lwd=2,ylim=c(0,ntot), xlab="time",ylab="Particles")
matpoints(out[,1],out[,2:3],type="l",lty=1,lwd=2); 	 
title(main="Two-compartment diffusion")
legend(0.5*tmax,ntot,legend=c("Left","Right"),col=c("Black","Red"),lty=1,lwd=2)

