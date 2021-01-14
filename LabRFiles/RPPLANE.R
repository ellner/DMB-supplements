#' Draw Phase Arrows
#'
#' Draws phase arrows for a planar system of ODEs dx/dt= f(x,y), dy/dt = g(x,y)
#' @param fun Function defining the vector field, fun(x,y,parms). See 'toggle' for an example. 
#' @param xlims Equivalent to 'xlim' in plot(); vector of length 2. 
#' @param ylims Equivalent to 'ylim' in plot(); vector of length 2.
#' @param resol Plotting resolution. Phase arrows are plotted at a regular (resol x resol) grid in the plotting region.  
#' @param col Color of arrows.
#' @param add Logical; equivalent to 'add' in plot()
#' @param parms Vector of parameters for fun().
#' @param jitter Logical: if TRUE, arrow locations are jittered slightly.  
#' @keywords phasearrows
#' @importFrom stats runif
#' @export
#' @examples
#' phasearrows()

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

#' Draw nullclines
#'
#' Draws nullclines for a planar system of ODEs dx/dt= f(x,y), dy/dt = g(x,y), using contour() to find where f and g are zero. 
#' @param fun Function defining the vector field, fun(x,y,parms). See 'toggle' for an example. 
#' @param xlims Equivalent to 'xlim' in plot(); vector of length 2. 
#' @param ylims Equivalent to 'ylim' in plot(); vector of length 2.
#' @param resol Controls grid size. Contouring is based on values computed at a regular (resol x resol) grid in the plotting region.  
#' @param add Logical; equivalent to 'add' in plot()
#' @param parms Vector of parameters for fun().
#' @keywords nullclines
#' @export
#' @examples
#' nullclines()
nullclines <- function(fun,xlims, ylims, resol=250, add=F,parms=NULL) {
  dx=0.01*diff(xlims); dy=0.01*diff(ylims);
  x <- matrix(seq(xlims[1]-dx,xlims[2]+dx, length=resol), byrow=F, resol,resol);
  y <- matrix(seq(ylims[1]-dy,ylims[2]+dy, length=resol), byrow=T, resol, resol);
  npts = resol*resol;
  z <- fun(x,y,parms);
  z1 <- matrix(z[1:npts], resol, resol);
  z2 <- matrix(z[(npts+1):(2*npts)], resol, resol);
  contour(x[,1],y[1,],z1,levels=c(0), drawlabels=F,add=add, col="blue",lwd=2);
  contour(x[,1],y[1,],z2,levels=c(0), drawlabels=F,add=T, col="forestgreen",lwd=2);
  title(xlab=x_label, col.lab="blue");
  title(ylab=y_label, col.lab="forestgreen");
  technicolorTitle(main_title, c("blue", "black","forestgreen"))

}

#' Draw trajectories from a grid of initial conditions
#'
#' Draws trajectories for a planar system of ODEs dx/dt= f(x,y), dy/dt = g(x,y), from a grid of initial conditions. 
#' @param fun Function defining the vector field, fun(x,y,parms). See 'toggle' for an example. 
#' @param xlim Equivalent to 'xlim' in plot(); vector of length 2. 
#' @param ylim Equivalent to 'ylim' in plot(); vector of length 2.
#' @param parms Vector of parameters for fun().
#' @param ngrid Controls grid size. Trajectories are started at a regular (rgrid x ngrid) grid in the plotting region.  
#' @param maxtime Solution time interval for trajectoris
#' @param color Color to use for plotting trajectories
#' @param add Logical; equivalent to 'add' in plot()
#' @param solver Choice of ODE solver method, any valid method for ode() from deSolve.
#' @keywords trajectories, grid
#' @export
#' @importFrom deSolve ode
#' @examples
#' gridStart()
gridStart=function(fun,xlim,ylim,parms,ngrid,maxtime=50,add=F,color="purple",solver="vode") {
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

#' Use multiple colors in the main title 
#'
#' Use multiple colors in main title. Called by Rpplane, not intended for access by users.
#' @param words Text for the title
#' @param colors Colors for the title
#' @param cex Equivalent to 'cex' in plot().
#' @keywords titles, color
#' @examples
#' technicolorTitle()
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

#' Newton's method to find equilibria for a system of autonomous ODEs, written in the format required by ode() in the deSolve package. 
#'
#' Finds equilibrium of a system of autonomous ODEs dX/dt= f(X) by Newton's method, starting from a specified initial point. 
#' @param func Function defining the vector field, in the format required by ode() in deSolve package. 
#' @param x0 Starting point for Newton iteration. 
#' @param parms Vector of parameters for func().
#' @param tol Tolerance for convergence check. Iteration terminates when the sum(f(X)^2) is <tol.  
#' @param niter Maximum number of iterations before return.
#' @param inc Increment for finite-difference estimates of function gradient. 
#' @param plotit plotit Logical; if TRUE add the equilibrium point to an existing plot
#' @keywords equilibrium, Newton's method
#' @export
#' @examples
#' newton()
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
        df[,j]=  (func(0,x+v,parms)[[1]] - func(0,x-v,parms)[[1]])/(2*inc)
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

#' Newton's method to find equilibria for a planar system of autonomous ODEs. 
#'
#' Finds equilibrium of a planar system of autonomous ODEs dX/dt= f(X) by Newton's method, starting from a specified initial point. 
#' @param fun Function defining the vector field, written in the format required by Rpplane() 
#' @param x0 Starting point for Newton iteration. 
#' @param parms Vector of parameters for func().
#' @param tol Tolerance for convergence check. Iteration terminates when the sum(f(X)^2) is <tol.  
#' @param niter Maximum number of iterations before return.
#' @param inc Increment for finite-difference estimates of function gradient. 
#' @param plotit plotit Logical; if TRUE add the equilibrium point to an existing plot
#' @keywords equilibrium, Newton's method
#' @export
#' @examples
#' newton2()
#'require(Rpplane)
#'fhn=function(x,y,parms) { # Fitzhugh-Nagumo function for Rpplane 
#'	dx = parms["c"]*(x-(1/3)*x^3-y+parms["j"]);
#'	dy = (1/parms["c"])*(x+parms["a"]-parms["b"]*y)
#'	return(c(dx,dy))
#'}
#' parms = c(c=8,a=0.7,b=0.8,j=0.3);
#' newton2(fun=fhn,x0=c(-1,-0.3),parms=parms)

newton2=function(fun,x0=NULL,parms=NULL,tol=1e-16,niter=40,inc=1e-6,plotit=FALSE) {
  x=x0; #initial x
  if (is.null(x0)){x = locator(n=1); x=c(x$x,x$y)};
  nx = 2; # length of state vector
  ######### Newton iteration loop: start
  for(i in 1:niter){
   y = fun(x[1],x[2],parms)
   df = matrix(0,nx,nx); # Compute df
   for(j in 1:nx) {
	#Increment vector for estimating derivative wrt jth coordinate
	v=rep(0,nx);
	v[j] = inc;
        df[,j]=  (fun(x[1]+v[1],x[2]+v[2],parms) - fun(x[1]-v[1],x[2]-v[2],parms))/(2*inc)
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

#' Draw stable and unstable manifolds of a saddle point
#'
#' Draws stable/unstable manifolds of a saddle point of planar system of ODEs, at or near a specified point. 
#' @param func Function defining the vector field, in the format required by lsoda/ode/rk4 in deSolve package.  
#' @param x0 Starting point for Newton iteration to find exact saddle location. 
#' @param parms Vector of parameters for func().
#' @param maxtime Time limit for forward/backward integration of the system to draw the manifolds.
#' @param dt Time increment for output of solution values to plot.
#' @param solver Solution method to be used by ode() for computing manifolds.
#' @keywords saddle point, stable manifold, unstable manifold
#' @importFrom deSolve ode
#' @importFrom grDevices dev.copy2pdf
#' @importFrom graphics arrows contour legend locator mtext par plot points strwidth title matplot
#' @importFrom utils select.list
#' @export
#' @examples
#' drawManifolds()
drawManifolds=function(func,parms,x0=NULL,maxtime=1000,dt=0.02,solver="vode") {
	xbar=newton(func,x0=x0,parms=parms,plotit=FALSE);
	x=xbar$x; df=xbar$df; V=eigen(df)$vectors; ev=eigen(df)$values;
	if (ev[1]*ev[2] > 0) {
	  cat("Fixed point is not a saddle \n");
	}else{
      i1=which(ev>0); i2=which(ev<0);
	  v1=V[,i1]; v2=V[,i2]; eps=0.002; xmax=.5+max(abs(x));
	  maxtime1=log(2500*xmax)/abs(ev[i1]); maxtime1=max(maxtime,maxtime1)
	  out1=ode(times=seq(0,maxtime1,dt),y=x+eps*v1,func=func,parms=parms,method=solver);
      points(out1[,2],out1[,3],type="l",lwd=2,col="red");
	  out2=ode(times=seq(0,maxtime1,dt),y=x-eps*v1,func=func,parms=parms,method=solver);
      points(out2[,2],out2[,3],type="l",lwd=2,col="red");

	  maxtime2=log(2500*xmax)/abs(ev[i2]); maxtime2=max(maxtime,maxtime2)
	  out3=ode(times=-seq(0,maxtime2,dt),y=x+eps*v2,func=func,parms=parms,method=solver);
      points(out3[,2],out3[,3],type="l",lwd=2,col="black");
	  out4=ode(times=-seq(0,maxtime2,dt),y=x-eps*v2,func=func,parms=parms,method=solver);
      points(out4[,2],out4[,3],type="l",lwd=2,col="black");
	  legend("topright",c("Stable","Unstable"),lty=1,lwd=2,col=c("black","red"),bty="n")

	}
}

#' Phase plane analysis for a planar system of ODEs dx/dt=f(x,y), dy/dt=g(x,y)
#'
#' Draws stable/unstable manifolds of a saddle point of planar system of ODEs, at or near a specified point. 
#' @param fun Function defining the vector field, fun(x,y,parms). See 'toggle' for an example. 
#' @param xlim Equivalent to 'xlim' in plot(); vector of length 2. 
#' @param ylim Equivalent to 'ylim' in plot(); vector of length 2.
#' @param parms Vector of parameters for fun().
#' @param add Logical; equivalent to 'add' in plot()
#' @param ngrid Grid resolution for drawing a grid of trajectories. Trajectories are started at a regular (ngrid x ngrid) grid in the plotting region. 
#' @param maxtime Time limit for forward/backward integration of trajectories.
#' @param dt Time increment for output of solution values to plot.
#' @param resol Grid resolution plotting phase arrows.
#' @param nullFac Grid resolution used for finding nullclines by contouring is resol*nullFac. 
#' @param solver Solution method to be used by ode() for computing manifolds.
#' @param x_lab Character: text for x-axis.
#' @param y_lab Character: text for y-axis.
#' @keywords Phase plane analysis
#' @importFrom deSolve ode
#' @importFrom zoom zm
#' @export
#' @examples
#' 
#' toggle=function(u,v,parms) { # toggle switch function for PPlane 
#'     du = -u + parms[1]/(1+v^parms[2])
#'     dv = -v + parms[1]/(1+u^parms[3])
#' 	return(c(du,dv))
#' }
#' 
#' fhn=function(x,y,parms) { # Fitzhugh-Nagumo function for PPlane 
#' 	dx = parms["c"]*(x-(1/3)*x^3-y+parms["j"]);
#' 	dy = (1/parms["c"])*(x+parms["a"]-parms["b"]*y)
#' 	return(c(dx,dy))
#' }
#' 
#' parms = c(c=10,a=0.7,b=0.8,j=0.3);
#' require(deSolve); require(zoom); 
#' Rpplane(fun=fhn,c(-3,3),c(-2,2),parms=parms)
#'
Rpplane=function(fun,xlim,ylim,parms=NULL,add=FALSE,ngrid=5,maxtime=100,dt=0.02,resol=25,nullFac=40,solver="vode",
    x_lab="Variable 1",y_lab="Variable 2",main_title=NULL) { 
    par(xaxs="i",yaxs="i");
    if(is.null(main_title)) main_title <<-c(x_lab,"vs",y_lab);
    x_label <<- x_lab; y_label <<- y_lab; 
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
    if(j=="3:") {xbar=newton2(fun,parms=parms,plotit=TRUE); }
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
    if(j=="7:") {drawManifolds(fun.lsoda,parms=parms,maxtime=maxtime,dt=dt,solver=solver) }
	if(j=="8:") {gridStart(fun.lsoda,xlim=xlim,ylim=ylim,parms=parms,ngrid=ngrid,add=add,solver=solver); add=TRUE}
    if(j=="9:") {zoom::zm(type="n")}
	if(j=="10") menu.go=0
	if(j=="11") dev.copy2pdf(file="PPlane.pdf")
  }
}
