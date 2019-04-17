# solve geometric growth model and plot the results

rm(list=ls(all=TRUE)); # clear past actions from R's memory 

# set the initial population size
initsize=10

# create matrix to hold results sizes and store initial size 
popsize=initsize; 

# create variable to hold the current population size
popnow=initsize;

# calculate population sizes and append to popsize
while(popnow<1000) {
    if(popnow<250){
        popnow=popnow*2;
    }else{
        popnow=popnow*1.5;
    }
    popsize=c(popsize,popnow);
}

tvals=1:length(popsize); 
plot(tvals, log(popsize),type="o",col="red",
xlab="Generation",ylab="Population size",pch=16,cex=1.25);
title(main="Geometric growth model"); 


