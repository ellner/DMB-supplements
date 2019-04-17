rm(list=ls(all=TRUE)); # clear past actions from R's memory 

# Iterate geometric growth model until population goes over 1000

# set initial size
initsize=10;

# initialize vector to hold results 
popsize=initsize; 

# create variable to hold the current population size
popnow=initsize; 

# calculate population sizes and append to popsize
while(popnow<1000) { 
        popnow=popnow*2;
        popsize=c(popsize,popnow);
}

tvals=1:length(popsize); 
plot(tvals, popsize,type="o",col="red",
xlab="Generation",ylab="Population size",pch=16,cex=1.25);

title(main="Geometric growth model"); 


