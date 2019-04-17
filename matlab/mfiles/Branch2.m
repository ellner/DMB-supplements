%Dynamic Models in Biology, Stephen Ellner and John Guckenheimer
%This m-file illustrates the use of conditional if-elseif-else statements 
%The program simulates a geometric growth model and plots the results

%Set the initial population size and number of generations
initsize=10;
r = 50;

%Create matrix to hold results sizes and store initial size
popsize=zeros(1,50);
popsize(1)=initsize; 

%Create variable to hold the current population size
popnow=initsize;

%Calculate population sizes and append to popsize
for i=2:r;
    %Determine the growth rate from the population size
    if(popnow<250);
    	popnow=popnow*2;
    elseif (popnow<500);
        popnow=popnow*1.5;
    else;
        popnow=popnow*0.95;
    end;
    popsize(i) = popnow;
end;

%Plot
xvals=[1:r]-1;
plot(xvals, log10(popsize),'-o');

xlabel('Generation','Fontsize',16); 
ylabel('Population size','Fontsize',16);
title('Limited growth model','Fontsize',18);
