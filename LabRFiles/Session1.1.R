rm(list=ls(all=TRUE)); # clear past actions from R's memory 

### Enter the data
### Note: a line starting with # is a comment, and R ignores these.  
cvals=c(1.1,1.4,4.1,5.5,5.4,10.7,6.0,22.0,19.8,26.7,31.4,30.9,27.1,40.2,36.1,36.8,31.6);
tvals=c(0,1,2,3,4,5,6,7,8,9,10,12,14,21,25,32,33);  

### Plot the data 
plot(tvals,log(cvals),xlab="Time (days)",ylab="log Chlorella density"); 

### Fit the exponential growth phase using lm()
y=log(cvals[1:5]); x=tvals[1:5]; 
fit = lm(y~x); summary(fit); 

### add the fit to the data plot
abline(fit); 

### parameter estimates from the fitted model 
x0=exp(fit$coef[1]); R=exp(fit$coef[2]); 
b=(1-1/R)/35; b; 

