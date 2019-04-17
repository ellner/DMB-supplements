rm(list=ls(all=TRUE)); # clear past actions from R's memory 

# initial population size
initsize=4; 

# create vector to hold results and store initial size 
popsize=rep(0,10); popsize[1]=initsize;

# calculate population size at times 2 through 10, write to Command Window
for (j in 2:10 ) { 
	popsize[j]=2*popsize[j-1];
	x=log(popsize[j]);
	cat(j,x,"\n");
}
