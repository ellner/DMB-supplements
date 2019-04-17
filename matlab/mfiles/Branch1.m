%Dynamic Models in Biology, Stephen Ellner and John Guckenheimer
%This m-file illustrates the use of conditional if-else statements 
%The program simulates a geometric growth model and plots the results

%Set the initial population size
initsize=10;

%Create matrix to hold results sizes and store initial size 
popsize=initsize; 

%Create variable to hold the current population size
popnow=initsize;

%Calculate population sizes and append to popsize
while popnow<1000;
    %Determine the growth rate from the population size
    if(popnow<250);
    	popnow=popnow*2;
    else;
        popnow=popnow*1.5;
    end;
    popsize=[popsize;popnow];
end;

%To plot we need to know how many values are in popsize
[r,c]=size(popsize);
xvals=(1:r)-1;
plot(xvals, log10(popsize),'-o');

xlabel('Generation','Fontsize',16); 
ylabel('Population size','Fontsize',16);
title('Geometric growth model','Fontsize',18);
