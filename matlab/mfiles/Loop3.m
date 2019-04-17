%Dynamic Models in Biology, Stephen Ellner and John Guckenheimer
%Third example of Matlab loops: solve geometric growth model until population goes over 1000

%Set initial size
initsize=10;

%Create matrix to hold results sizes and store initial size 
popsize=initsize; 

%Create variable to hold the current population size
popnow=initsize; 

%Calculate population sizes and append to popsize
%Note that the length of the vector popsize grows
while(popnow<1000); 
	popnow=popnow*2;
    popsize=[popsize;popnow];
end;

%To plot we need to know how many values are in popsize
[r,c]=size(popsize);
xvals=(1:r)-1;
plot(xvals, popsize,'-o','Color','red');

xlabel('Generation','Fontsize',16); 
ylabel('Population size','Fontsize',16);
title('Geometric growth model','Fontsize',18);
