#   Solve the spatial FHN model by method of lines 
#   For simplicity, we assume a square array of grid points

rm(list=ls(all=TRUE)); # clear past actions from R's memory 
library(deSolve);
setwd("e:/classes/dynmod/LabRfiles"); 

sfn2=function(t,y,parms) {
	e=parms[1]; b=parms[2]; dx2=parms[3];
	u=y[1:ny2]; v=y[(ny2+1):ny]; 
	umat=matrix(u,nxy,nxy); 
    
    uec = cbind(umat[,1],umat,umat[,nxy]); # u with extra columns
	uer = rbind(umat[1,],umat,umat[nxy,]); # u with extra rows 
	ul = uer[3:(nxy+2),]+uer[1:nxy,]+uec[,1:nxy]+uec[,3:(nxy+2)]-4*umat;
	uf = (u-0.3333*u*u*u-v)/e + dx2*matrix(ul,ncol=1);
	vf = e*(u + b-0.5*v);
  
   	return(list(dy=c(uf,vf)))
}	

# Parameters and initial conditions. 
nxy = 60; ny2=nxy*nxy; ny=2*ny2;  
dx=3; dx2 = 1/(dx*dx); e = 0.1; b = 0.67;
u=5*outer(rep(1,nxy),((1:nxy)-0.4*nxy)/nxy); 
v = t(u); 

pdf(file="sfn.pdf",onefile=TRUE);

y0=c(as.vector(u),as.vector(v)); h  = 0.05; 
for(i in 0:500) { 
   	u=matrix(y0[1:ny2],nxy,nxy); 
    image(1:dim(u)[2],1:dim(u)[1],t(u),col=rainbow(n=1000,start=0,end=0.7))
    title(main=paste("Range=",round(min(u),digits=2),"to", round(max(u),digits=2))); 
    out=rk4(y0,times=c(0:10)*h,func=sfn2,parms=c(e,b,dx2)); # plots every 10*h time units!! 
	y0=out[nrow(out),-1]; 
    cat(i,"\n"); 
}
dev.off();



