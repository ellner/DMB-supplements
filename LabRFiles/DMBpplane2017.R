# Last modified May 17, 2016 by S.P. Ellner
# 
# Defines the function Rpplane to do phase-plane analysis of
# two-dimensional systems of ODEs. 
#
# Based on the script pplane.r written by Daniel Kaplan,
# Dept. of Mathematics, Macalester College, kaplan@macalester.edu
#
# Modifications by S. Ellner and J. Guckenheimer for use with
# Dynamic Models in Biology by S.P. Ellner and J. Guckenheimer,
# Princeton University Press (2006)

# Revision history (SPE unless otherwise notes) 
# 2008-2010: removed Windows-specific functions to make it macOS compatible. 

# Feb 2013: added 'save plot to PDF' option

# January 2016: changed from lsoda() to ode() for calculating trajectories,
# and added the ability to choose any ODE solution method available in ode().
# By default, Rpplane uses vode, the Livermore ODEPACK stiff solver.
# If that is too slow, try setting solver = "lsoda" in the call to Rpplane.

# Oct 2016 (Evan G. Cooch): adding aesthetic features, like changing colour of titles

# May 2017: make phasearrows() handle vector fields not defined at some points. 

cat("\n");
cat("DMBpplane script now loaded...\n");
cat("(version date - May 17, 2016)  \n\n");

require(deSolve); require(zoom); par(xaxs="i",yaxs="i");

x_label <<- "Variable 1"; y_label <<- "Variable 2"; 
main_title <<-c("Variable 1", "vs", "Variable 2"); 

############################################################################################### 
#  Plotting functions, largely the same as in the original pplane.r by Daniel Kaplan
###############################################################################################
phasearrows <- function(fun,xlims,ylims,resol=25, col='black',add=F,parms=NULL,jitter=FALSE) {
  if (add==F) {
     plot(1,xlim=xlims, ylim=ylims, type='n',xlab="",ylab="");
     title(xlab=x_label, col.lab="blue");
     title(ylab=y_label, col.lab="forestgreen");
     technicolorTitle(main_title, c("blue", "black","forestgreen"))
  }
  x <- matrix(seq(xlims[1],xlims[2], length=resol), byrow=T, resol,resol);
  y <- matrix(seq(ylims[1],ylims[2], length=resol),byrow=F, resol, resol);
  npts <- resol*resol;
  #### the lines below 'jitter' the location of arrows.
  if(jitter) {
    xspace <- abs(diff(xlims))/(resol*10);
    yspace <- abs(diff(ylims))/(resol*10);
    x <- x + matrix(runif(npts, -xspace, xspace),resol,resol);
    y <- y + matrix(runif(npts, -yspace, yspace),resol,resol);
  }
  z <- fun(x,y,parms);
  z1 <- matrix(z[1:npts], resol, resol);
  z2 <- matrix(z[(npts+1):(2*npts)], resol, resol);
  
  # an obscure but effective scaling to make arrows at a good size for plotting 
  max_x <- max(abs(z1),na.rm=TRUE); 
  max_y <- max(abs(z2),na.rm=TRUE);
  dt <- min(abs(diff(xlims))/max_x, abs(diff(ylims))/max_y)/resol;
  lens <- sqrt(z1^2 + z2^2);
  lens2 <- lens/max(lens,na.rm=TRUE);
  arrows(x, y, x + dt*z1/(lens2+.1), y + dt*z2/(lens2+.1), length=.04, col=col);
}

showcontours <- function(fun,xlims, ylims,resol=250,add=F, colors=c('red', 'blue'),parms=NULL) {
  dx=0.01*diff(xlims); dy=0.01*diff(ylims);
  x <- matrix(seq(xlims[1]-dx,xlims[2]+dx, length=resol), byrow=F, resol,resol);
  y <- matrix(seq(ylims[1]-dy,ylims[2]+dy, length=resol),byrow=T, resol, resol);
  npts = resol*resol;
  z <- fun(x,y,parms);
  z1 <- matrix(z[1:npts], resol, resol);
  z2 <- matrix(z[(npts+1):(2*npts)], resol, resol);
  contour(x[,1],y[1,],z1, add=add, col=colors[1]);
  contour(x[,1],y[1,],z2, add=T, col=colors[2]);
  technicolorTitle(main_title, c("blue", "black","forestgreen"))
}

nullclines <- function(fun,xlims, ylims, resol=250, add=F,parms=NULL) {
  dx=0.01*diff(xlims); dy=0.01*diff(ylims);
  x <- matrix(seq(xlims[1]-dx,xlims[2]+dx, length=resol), byrow=F, resol,resol);
  y <- matrix(seq(ylims[1]-dy,ylims[2]+dy, length=resol),byrow=T, resol, resol);
  npts = resol*resol;
  z <- fun(x,y,parms);
  z1 <- matrix(z[1:npts], resol, resol);
  z2 <- matrix(z[(npts+1):(2*npts)], resol, resol);
  contour(x[,1],y[1,],z1,levels=c(0), drawlabels=F,add=add, col="blue",lwd=2);
  contour(x[,1],y[1,],z2,levels=c(0), drawlabels=F,add=T, col="forestgreen",lwd=2);
  technicolorTitle(main_title, c("blue", "black","forestgreen"))

}

grid=function(fun,xlim,ylim,parms,ngrid,maxtime=50,add=F,color="purple",solver="vode") {
	 if (add==F) {
	     plot(1,xlim=xlim, ylim=ylim, type='n',,xlab="",ylab="");
     title(xlab=x_label, col.lab="blue")
     title(ylab=y_label, col.lab="forestgreen");
	}
	xvals=seq(xlim[1],xlim[2],length=ngrid);
	yvals=seq(ylim[1],ylim[2],length=ngrid);
	for(i in 1:ngrid) {
	for(j in 1:ngrid) {
	out=ode(times=seq(0,maxtime,.02),y=c(xvals[i],yvals[j]),func=fun,parms=parms,method=solver);
	points(out[,2],out[,3],type="l",lwd=1,col=color);
	out=ode(times=-seq(0,maxtime,.02),y=c(xvals[i],yvals[j]),func=fun,parms=parms,method=solver);
	points(out[,2],out[,3],type="l",lwd=1,col=color);

	}}
}

# Use multiple colors in main title, matched to the colors used to plot the nullclines
# Provided by Evan G. Cooch 
technicolorTitle <- function(words, colors, cex=1) {
    widths <- strwidth(words,cex=cex)
    spaces <- rep(strwidth(" ",cex=cex), length(widths)-1)
    middle <- mean(par("usr")[1:2])
    total <- sum(widths) + sum(spaces)
    start <- c(0,cumsum(widths[-length(widths)] + spaces))
    start <- start + middle - total/2
    mtext(words, 3, 1, at=start, adj=0, col=colors,cex=cex)
    }

############################################################################################### 
#  Utility functions: Newton's method, finding Jacobians, drawing Stable/Unstable manifolds
###############################################################################################

# Newton's method to find equilibria of vector field.
# func() must have the same input arguments and returns as for lsoda/rk4.
# Inputs:
#   x0 = intial guess at equilibrium. If x0 is not supplied in the call,
#        the user chooses it from the current graphics device via locator()
#         and the equilibrium is plotted to the same device. Plotting
#         symbol is closed/open=stable/unstable, circle/triangle=eigenvalues imaginary/real.
#   tol= Convergence tolerance
#   niter = Maximum number of iterations
#   inc = finite-difference increment for derivative estimates
# Coded 5/25/06 by SPE based on Matlab version by JG

newton=function(func,x0=NULL,parms=NULL,tol=1e-16,niter=40,inc=1e-6,plotit=FALSE) {
  x=x0; #initial x
  if (is.null(x0)){x = locator(n=1); x=c(x$x,x$y)};
  nx = length(x); # length of state vector
  ######### Newton iteration loop: start
  for(i in 1:niter){
   y = func(0,x,parms)[[1]]
   df = matrix(0,nx,nx); # Compute df
   for(j in 1:nx) {
	#Increment vector for estimating derivative wrt jth coordinate
	v=rep(0,nx);
	v[j] = inc;
        df[,j]=  (func(t,x+v,parms)[[1]] - func(t,x-v,parms)[[1]])/(2*inc)
    }
    if (sum(y^2) < tol){  #check for convergence
        if(plotit){
	     ev=eigen(df)$values;
		 pch1=1+as.numeric(Im(ev[1])!=0); # 1 if real, 2 if complex
		 pch2=1+sum(Re(ev)>0); #number of positive eigenvalues
	     pchs=matrix( c(17,5,2,16,1,1),2,3,byrow=T);
 	     par(xpd=TRUE); 	
	     points(x[1],x[2],type="p",pch=pchs[pch1,pch2],cex=1.5,lwd=2)
	     par(xpd=FALSE);
 	}
	cat("========================================","\n")
	cat("Fixed point at (x,y) = ",x,"\n");
	cat("Jacobian Df=","\n"); print(df); cat("Eigenvalues","\n"); print(eigen(df)$values);
	cat("Eigenvectors", "\n"); print(eigen(df)$vectors);
	cat("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^","\n","\n");
	x0=as.character(round(x[1],digits=3)); y0=as.character(round(x[2],digits=3))
	evals=round(eigen(df)$values,digits=3);
	e1=as.character(evals[1]); e2=as.character(evals[2]);
	return(list(x=x,df=df))
    }	
    x = x - solve(df,y) # one more step if needed
    cat(i, x, "\n") #print out the next iterate
  }
  ######### Newton iteration loop: end
 cat("Convergence failed");

}

# Compute Jacobian of a planar vector field at a point (x,y),
# either input or chosen with locator().
jacobianAtXY <- function(fun,x=NULL, y=NULL,inc=1e-7){
  if (is.null(x)|is.null(y)) {
    x0 <- locator(n=1); x <- x0$x; y <- x0$y;
  }
  foo <- fun(x,y); h = inc;
  foox <- fun(x+h,y); fooy <- fun(x,y+h);
  A <- (foox[1] - foo[1])/h;
  B <- (fooy[1] - foo[1])/h;
  C <- (foox[2] - foo[2])/h;
  D <- (fooy[2] - foo[2])/h;
  return(matrix( c(A,B,C,D ),2,2,byrow=T))
}

# Draw stable manifolds for a saddle at or near the point x0 (2-vector). 
# fun.lsoda() must have the same input arguments and returns as ode/lsoda/rk4.
# Exact location of the saddle, and the Jacobian there, are found by a call to newton(). 
DrawManifolds=function(fun.lsoda,parms,x0=NULL,maxtime=1000,dt=0.02,solver="vode") {
	xbar=newton(fun.lsoda,x0=x0,parms=parms,plotit=FALSE);
	x=xbar$x; df=xbar$df; V=eigen(df)$vectors; ev=eigen(df)$values;
	if (ev[1]*ev[2] > 0) {
	  cat("Fixed point is not a saddle \n");
	}else{
      i1=which(ev>0); i2=which(ev<0);
	  v1=V[,i1]; v2=V[,i2]; eps=0.002; xmax=.5+max(abs(x));
	  maxtime1=log(2500*xmax)/abs(ev[i1]); maxtime1=max(maxtime,maxtime1)
	  out1=ode(times=seq(0,maxtime1,dt),y=x+eps*v1,func=fun.lsoda,parms=parms,method=solver);
      points(out1[,2],out1[,3],type="l",lwd=2,col="red");
	  out2=ode(times=seq(0,maxtime1,dt),y=x-eps*v1,func=fun.lsoda,parms=parms,method=solver);
      points(out2[,2],out2[,3],type="l",lwd=2,col="red");

	  maxtime2=log(2500*xmax)/abs(ev[i2]); maxtime2=max(maxtime,maxtime2)
	  out3=ode(times=-seq(0,maxtime2,dt),y=x+eps*v2,func=fun.lsoda,parms=parms,method=solver);
      points(out3[,2],out3[,3],type="l",lwd=2,col="black");
	  out4=ode(times=-seq(0,maxtime2,dt),y=x-eps*v2,func=fun.lsoda,parms=parms,method=solver);
      points(out4[,2],out4[,3],type="l",lwd=2,col="black");
	  legend("topright",c("Stable","Unstable"),lty=1,lwd=2,col=c("black","red"),bty="n")

	}
}

# GUI-challenged version of Polking's pplane.m for planar vector fields.
# The function fun() defining the vector field must be in the format 
# fun=function(x,y,parms)   
#    x and y are vectors of values of the state variables
#    parms is a vector of parameter values 
#    The returned value must be a vector, c(dxdt,dydt).
# See the examples at the end of this script. 
# Control parameters: 
#   maxtime, dt: total duration and time-increment for drawing trajectories.
#   ngrid: grid of trajectories is drawn for (ngrid*ngrid) initial conditions. 
#   resol: phase arrows are drawn at (resol*resol) points
#   nullFAC: nullclines are computed based on vector field at (nullFac*resol) by (nullFac*resol) grid of points
Rpplane=function(fun,xlim,ylim,parms=NULL,add=FALSE,ngrid=5,maxtime=100,dt=0.02,resol=25,nullFac=40,solver="vode") {
   fun.lsoda=function(t,y,p) {dx=fun(y[1],y[2],parms=p); return(list(dx))}
   menu.go=1; while(menu.go>0) {
    jl= select.list(c("1:  Phase arrows",
	"2:  Nullclines",
	"3:  Find fixed point (click on plot)",
    "4:  Start Forward trajectory [purple] (click on plot)",
	"5:  Start Backward trajectory [orange] (click on plot)",
	"6:  Extend current trajectory",
	"7:  Local S/U manifolds for saddle (click on plot)",
    "8:  Draw a grid of trajectories",
	"9:  Zoom",
	"10: Exit",
    "11: Save as PDF in working directory"),
	title = "R-pplane: Select action");
	j=substr(jl,1,2);
    if(j=="1:") {phasearrows(fun=fun,xlims=xlim,ylims=ylim,parms=parms,resol=resol,add=add); add=TRUE}
	if(j=="2:") {nullclines(fun=fun,xlims=xlim,ylims=ylim,resol=nullFac*resol,parms=parms,add=add); add=TRUE};
    if(j=="3:") {xbar=newton(fun.lsoda,parms=parms,plotit=TRUE); }
	if(j=="4:") {
		x=locator(n=1);
		out=ode(times=seq(0,maxtime,dt),y=c(x$x,x$y),func=fun.lsoda,parms=parms,method=solver);
		points(out[,2],out[,3],type="l",lwd=2,col="purple");
	}
	if(j=="5:") {
		x=locator(n=1);
		out=ode(times=seq(0,-maxtime,-dt),y=c(x$x,x$y),func=fun.lsoda,parms=parms,method=solver);
		points(out[,2],out[,3],type="l",lwd=2,col="orange3");
	}
	if(j=="6:") {
		nt=dim(out)[1]; x=out[nt,2:3]; times=out[,1];
	 	dt=times[2]-times[1]; pcol=ifelse(dt>0,"purple","grey20");
		out=ode(times=out[,1],y=out[nt,2:3],func=fun.lsoda,parms=parms,method=solver);
		  points(out[,2],out[,3],type="l",lwd=2,col=pcol);
	}
    if(j=="7:") {DrawManifolds(fun.lsoda,parms=parms,maxtime=maxtime,dt=dt,solver=solver) }
	if(j=="8:") {grid(fun.lsoda,xlim=xlim,ylim=ylim,parms=parms,ngrid=ngrid,add=add,solver=solver); add=TRUE}
    if(j=="9:") {zm(type="n")}
	if(j=="10") menu.go=0
	if(j=="11") dev.copy2pdf(file="PPlane.pdf")
  }
}


###########################################################
# Examples below 
###########################################################

# toggle switch function for PPlane 
toggle=function(u,v,parms) {
    du = -u + parms[1]/(1+v^parms[2])
    dv = -v + parms[1]/(1+u^parms[3])
	return(c(du,dv))
}

# Compare with: Toggle switch function for ode or lsoda!  
Toggle=function(t,y,parms) {
	u=y[1]; v=y[2];
	du= -u + parms[1]/(1+v^parms[2]);
	dv= -v + parms[1]/(1+u^parms[3]);
	dX=c(du,dv);
	return(list(dX));
}


# A simple patch dynamics problem
# Species 1 better competitor than species 2
PatchTest = function(psi1,psi2,parms) {
    m1 <- parms[1]
    m2 <- parms[2]
    e <- parms[3]
    dpsi1 = m1*psi1*(1-psi1)-e*psi1
    dpsi2 = m2*psi2*(1-psi1-psi2)-e*psi2-m1*psi1*psi2
    return(c(dpsi1,dpsi2))
}

#parms <- c(0.2,0.9,0.1)
#y_label <<- "Species 2"; x_label <<- "Species 1"; # labels for x and y axes 
#main_title <- c("Prop(1)", "vs", "Prop(2)") # set words you want in main title for phase-plane
#Rpplane(fun=PatchTest,xlim=c(0,1),ylim=c(0,1),parms=parms,resol=25,dt=0.02,solver="vode")
