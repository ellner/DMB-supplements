rm(list=ls(all=TRUE)); # clear past actions from R's memory 

### Read in the data from a text file
### Note: edit the path!   
X=read.table("e:/classes/DynMod/LabRfiles/ChlorellaStart.txt"); 
tvals=X[,1]; cvals=X[,2]; 

## use lm() to estimate model parameters 
fit=lm(log(cvals[1:5])~tvals[1:5]); 
x0=exp(fit$coef[1]); R=exp(fit$coef[2]); 
b=(1-1/R)/35; 

## Solve the model using a for-loop.
## Don't panic -- we'll learn about for-loops!  
xvals=numeric(34); xvals[1]=x0; 
for(j in 2:34) {xvals[j]=R*xvals[j-1]*(1-b*xvals[j-1])}

## plot data and model solutions
plot(tvals,cvals,type="p",pch=16,cex=1.5, col="green",xlab="Time (days)", ylab="Algal density"); 
points(0:33,xvals,type="o",pch=16,cex=1,col="black",lwd=2); 
