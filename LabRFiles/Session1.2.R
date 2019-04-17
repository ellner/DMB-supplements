
rm(list=ls(all=TRUE)); # clear past actions from R's memory 

### Read in the data from a text file
### Note: edit the path!   
X=read.table("e:/classes/DynMod/LabRfiles/ChlorellaStart.txt"); 
tvals=X[,1]; cvals=X[,2]; 

### Plot the data 
par(cex=1.5,cex.main=0.9); 
plot(tvals,log(cvals),xlab="Time (days)", ylab="Ln Algal density",pch=16); 
title(main="Data from Fussmann et al. (2000) system");
# Note: 
# xlab and ylab are x and y axis labels, pch is "plotting character"
# cex is 'character expansion' - cex=1.5 increases symbol & label sizes by 50%
# cex.main sets the character expansion for the main title of the plot 

### Fit the exponential growth phase using lm(), and
### add the fit to the plot 
y=log(cvals[1:5]); x=tvals[1:5]; 
fit = lm(y~x); abline(fit); 

# Have the regression equation 'display itself' on the graph
c1=round(fit$coef[1],digits=3); 	# intercept	
c2=round(fit$coef[2],digits=3); 	# slope
text(15,2,paste("ln Algae=",c1,"+",c2,"t")); 
# You can use ?round, ?text and ?paste to read about these commands
# for working with plots  

### parameter estimates from the fitted model 
x0=exp(fit$coef[1]); R=exp(fit$coef[2]); 
b=(1-1/R)/35; 

### Write parameter values to the console, see ?cat for details.
cat(x0,R,b,"\n"); 


