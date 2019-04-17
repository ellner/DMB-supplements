% solve geometric growth model until population goes over 1000

% set initial size
initsize=10;

% create matrix to hold results sizes and store initial size 
popsize=initsize; 

% create variable to hold the current population size
popnow=initsize; 

% calculate population sizes and append to popsize
while(popnow<1000); 
	popnow=popnow*2;
    popsize=[popsize;popnow];
end;

% to plot we need to know how many values are in popsize
[r,c]=size(popsize);
xvals=(1:r)-1;
plot(xvals, popsize,'-o','Color','red');

xlabel('Generation','Fontsize',16); 
ylabel('Population size','Fontsize',16);
title('Geometric growth model','Fontsize',18);
