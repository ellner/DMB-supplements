rm(list=ls(all=TRUE)); # clear past actions from R's memory 

require(deSolve); 

repress=function(t,y,parms){
   dy = rep(0,6);
   alpha=parms[1]; alpha0=parms[2]; b=parms[3]; n=parms[4]; 
   dy[1] = -y[1] + alpha/(1+y[6]^n)+ alpha0;
   dy[2] = -y[2] + alpha/(1+y[4]^n)+ alpha0;
   dy[3] = -y[3] + alpha/(1+y[5]^n)+ alpha0;
   dy[4] = -b*(y[4]-y[1]);
   dy[5] = -b*(y[5]-y[2]);
   dy[6] = -b*(y[6]-y[3]);
   return(list(dy))
}

