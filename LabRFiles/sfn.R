rm(list=ls(all=TRUE)); # clear past actions from R's memory 

setwd("e:/classes/dynmod/LabRfiles"); 

sfn=function(u,v) {
	ny = dim(u)[1]; nx = dim(u)[2];
	uec = cbind(u[,1],u,u[,nx]); # u with extra columns
	uer = rbind(u[1,],u,u[ny,]); # u with extra rows 
	ul = uer[3:(ny+2),]+uer[1:ny,]+uec[,1:nx]+uec[,3:(nx+2)]-4*u;
	uf = (u-0.3333*u*u*u-v)/e + dx2*ul;
	vf = e*(u + b-0.5*v);
	return(list(uf=uf,vf=vf)); 
} 

sfn.run=function(nsteps,init) {
   u=init$u; v=init$v; 
   for (i in 1:nsteps){
	out=sfn(u,v);
	u = u+h*out$uf;	v = v+h*out$vf;
	if (i%%5 == 1) image(1:dim(u)[2],1:dim(u)[1],t(u),col=rainbow(n=1000,start=0,end=0.7))
    cat(i,"\n"); 
   }
   return(list(u=u,v=v)); 
}

# Parameters and initial conditions: set 1
nxy = 100; dx2 = 1/4; h  = 0.04; e = 0.1; b = 0.67;
u=5*outer(rep(1,nxy),((1:nxy)-0.4*nxy)/nxy); 
v = t(u); 
init1=list(u=u,v=v); 

pdf(file="sfn.pdf",onefile=TRUE);
nsteps=1000; # or however many time steps you want 
out=sfn.run(nsteps,init1);
dev.off();

# Parameters and initial conditions: set 2
nxy = 60; dx2 = 1/9; h  = 0.04; e = 0.1; b = 0.67;
v=5*outer( (1:nxy)-0.4*nxy, rep(1,2*nxy))/nxy; 
u = 5*matrix(1,nxy,1)%*%((1:nxy)-0.4*nxy)/nxy;
u = cbind(u,u[,nxy:1])
init2=list(u=u,v=v); 

pdf(file="sfn.pdf",onefile=TRUE);
nsteps=1000; # or however many time steps you want 
out=sfn.run(nsteps,init2);
dev.off(); 

