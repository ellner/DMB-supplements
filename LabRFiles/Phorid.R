s0 = 0.0035; s1 = 0.035; m0 = 0.116; m1 = 0.036;

N = 500; # number of rows = number of columns in grid representing the forest
Trees = matrix(0,N,N);
nSites = (N-2)^2;

Trees[2:(N-1),2:(N-1)] = (runif(nSites) < 0.5); 

here = 2:(N-1); less=1:(N-2); more=3:N;


for(j in 1:2500) {
MooreSum = Trees[less,here] +Trees[more,here] +
		Trees[less,less] +Trees[here,less]+Trees[more,less] +
		Trees[less,more] + Trees[here,more] + Trees[more,more];
pdie = m0 + m1*MooreSum;
pgro = s0 + s1*MooreSum; 
die = matrix( (runif(nSites)< pdie), N-2,N-2) 
gro = matrix( (runif(nSites)< pgro), N-2,N-2) 

newTrees = oldTrees = Trees[here,here];
newTrees[(oldTrees==1)&(die==1)]<-0; 
newTrees[(oldTrees==0)&(gro==1)]<-1;
Trees[here,here] <- newTrees;
cat(j,"\n")
}


