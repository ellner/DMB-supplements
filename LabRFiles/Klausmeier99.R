rm(list=ls(all=TRUE)); 
setwd("e:/classes/dynmod"); 
require(deSolve); 

K99 = function(t,y,parms) {
    a=parms[1]; v= parms[2]; m=parms[3];   
    w=y[1:nx]; n=y[-(1:nx)]; 
    wplus=c(w[nx],w,w[1]);
    nplus=c(n[nx],n,n[1]); 
    dwdx = (wplus[3:(nx+2)]-wplus[1:nx])/(2*dx);  
    d2ndx2 = (nplus[3:(nx+2)]+ nplus[1:nx] - 2*nplus[2:(nx+1)])/(dx^2);  
    dw = a - w - w*(n^2) + v*dwdx; 
    dn = w*(n^2) - m*n + d2ndx2; 
    return(list(dy=c(dw,dn)))
} 


############ Solve in 1D 

# parameters 
nx = 100; dx = 1; nx2 = nx^2; # nx by nx grid, at distance of dx apart. 
parms=c(a=2,v=80,m=0.5); 

# initial conditions 
n0=seq(0,1,length=nx); n0=(4*n0*(1-n0))^2; 
w0 = rep(parms["a"],nx);  
y=c(w0,n0); matplot(1:nx,cbind(w0,n0)); 

# solve 
out=ode(y,times=seq(0,100,by=.2),K99,parms=parms,method="vode"); 

# plot 
ymax=max(out[,-1]); 
pdf(file="K99.pdf",onefile=TRUE);
for(j in 1:nrow(out)){
     matplot(1:nx,matrix(out[j,-1],ncol=2,byrow=FALSE),col=c("blue","green3"),xlab="Location", ylab="Plants and water",lty=1,
        type="l",lwd=2,ylim=c(0,ymax),main=round(out[j,1],digits=3)); 
}
dev.off(); 



