%Dynamic Models in Biology, Stephen Ellner and John Guckenheimer
%Secoond example m file for the Matlab Introductory materials

%Load data from file
X=load('ChlorellaGrowth.txt');

%Define vectors by slicing data
Light=X(:,1); rmax=X(:,2); 

%Repeat computations from Intro1
C=polyfit(Light,rmax,1)
xvals=0:120; 
yhat=polyval(C,xvals);
plot(Light,rmax,'+',xvals,yhat);

%Annotate plot with labels, title and text
xlabel('Light intensity (uE/m2/s)','Fontsize',14); 
ylabel('Maximum growth rate (1/d)','Fontsize',14);
title('Data from Fussmann et al. (2000)','FontSize',14);
text(30,3.5,'y= 1.581+ 0.0136*x','Fontsize',14);


