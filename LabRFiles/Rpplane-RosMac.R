rm(list=ls(all=TRUE))
setwd("e:\\Projects\\pplane package"); # edit this to work on your computer 
source("DMBpplane2019.R"); 

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Function to compute state variable derivatives. 
# Note that the state variables are sent in as two separate arguments, 
# not as a vector, and that time t is not an argument.
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

RosMac=function(n1,n2,parms) {
    r1=parms[1]; d2=parms[2]; 
    a1=parms[3]; a2=parms[4];
    B=parms[5]; K=parms[6]; 
    dx1=r1*n1*(1-n1/K) - (a1*n1*n2)/(B+n1);
    dx2=(a2*n1*n2)/(B+n1)-d2*n2;	
    return(c(dx1,dx2)) ;
}
y_label <- "Prey"; x_label <- "Predators"; main_title <- "Rosenzweig-MacArthur model" 

parms=c(r1=1,d2=1,a1=2,a2=2,B=200,K=100) 
Rpplane(fun=RosMac, xlim=c(-0.02,220),ylim=c(-0.01,250),parms=parms); 

