rm(list=ls(all=TRUE)); 
setwd("e:/classes/dynmod/LabRfiles"); 
require(deSolve); 

K992D = function(t,y,parms) {
    a=parms[1]; v= parms[2]; m=parms[3];   
    
    ###### Water first, vegetation second 
    w=matrix(y[1:nx2],nx,nx); n=matrix(y[-(1:nx2)],nx,nx); 

    # spatial Laplacian of n 
    nec = cbind(n[,nx],n,n[,1]); # n with extra columns, periodic BC
	ner = rbind(n[nx,],n,n[1,]); # n with extra rows, periodic BC 
    nl = (ner[3:(nx+2),]+ner[1:nx,]+nec[,1:nx]+nec[,3:(nx+2)]-4*n)/(dx^2);
    
    # derivative of w with respect to x 
    wec = cbind(w[,nx],w,w[,1]); # w with extra columns, periodic BC 
    dwdx = (wec[,3:(nx+2)]-wec[,1:nx])/(2*dx);  
    
    dw = matrix(a - w - w*(n^2) + v*dwdx, ncol=1); 
    dn = matrix(w*(n^2) - m*n + nl, ncol=1); 
    
   	return(list(dy=c(dw,dn)))
}	

# parameters 
nx = 60; dx = 1; nx2 = nx^2; # nx by nx grid, at distance of dx apart. 
parms=c(a=2,v=80,m=0.5); 

# initial conditions 
n0=seq(0,1,length=nx); n0=(4*n0*(1-n0))^2; 
n0 = 0.5+0.25*outer(n0,n0); 
w0 = matrix(parms["a"],nx,nx);  

graphics.off(); 
image(1:nx,1:nx,t(n0)/max(n0),main="Vegetation",col=rainbow(n=1000,start=0,end=0.7)); 

h=0.01; y0=c(matrix(w0,ncol=1),matrix(n0,ncol=1)); 

pdf(file="K992D.pdf",onefile=TRUE);
for(i in 0:1000) { 
    par(mfrow=c(2,1)) 
    n=matrix(y0[-(1:nx2)],nx,nx);  w=matrix(y0[1:nx2],nx,nx); 
    image(1:nx,1:nx,t(n)/max(n),col=rainbow(n=1000,start=0,end=0.7))
    title(main=paste("Vegetation. Range=",round(min(n),digits=2),"to", round(max(n),digits=2))); 
    image(1:nx,1:nx,t(w)/max(w),col=rainbow(n=1000,start=0,end=0.7))
    title(main=paste("Water Range=",round(min(w),digits=2),"to", round(max(w),digits=2))); 
    cat(i,"\n"); 
    out=rk4(y0,times=c(0:10)*h,func=K992D,parms=parms); 
	y0=out[nrow(out),-1]; 
}
dev.off();

