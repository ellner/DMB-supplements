%Dynamic Models in Biology, Stephen Ellner and John Guckenheimer
%First example m file for the Matlab Introductory materials

%Define row vectors
Light=[20 20 20 20 21 24 44 60 90 94 101]
rmax=[1.73 1.65 2.02 1.89 2.61 1.36 2.37 2.08 2.69 2.32 3.67]

%Compute least squares linear fit to rmax as a function of Light
C=polyfit(Light,rmax,1)

%Prepare to plot linear fit by defining mesh and values
xvals=(0:120); 
yhat=polyval(C,xvals);

%Plot it
plot(Light,rmax,'+',xvals,yhat);


