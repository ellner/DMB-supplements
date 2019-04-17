rm(list=ls(all=TRUE)); # clear past actions from R's memory 
graphics.off(); 

nt = 5000000;
p = 0.6;      # probability of heads
rdata = as.numeric(runif(nt) < p); # heads = 1, tails=0 

######### Compute the run lengths for heads and tails 
z = which(rdata[1:(nt-1)]!=rdata[(2:nt)]); # where is the next value different from the current value? 
lengths=diff(c(0,z));   # run length = distance between places where value changes
values=rdata[z];        # what is the value just before the change? 

rh = lengths[values==1];  # run lengths of heads
rt = lengths[values==0];  # run lengths of tails 

############# Plot results 
par(mfrow=c(2,1),cex.axis=1.35,cex.lab=1.35); 
hhist = hist(rh,(0:max(rh))+0.5,xlim=c(1,max(rh)));     # heads histogram
#>> for exercise 11.4: compute theoretical predictions, and overlay them onto this plot 

thist = hist(rt,(0:max(rt)+0.5),xlim=c(1,max(rt)));     # tails histogram
#>> for exercise 11.4: compute theoretical predictions, and overlay them onto this plot 

### Next, plot the results on semi-log scale 
### If run-lengths have a geometric distribution, this should give
### a straight-line plot. 
### Do ?hist and read the Value section to see what's being plotted

dev.new(); 
par(mfrow=c(2,1),cex.axis=1.35,cex.lab=1.35);  
plot(hhist$mids,log(hhist$density),xlab="Run length",ylab="Frequency"); 
fit=lm(log(hhist$density[1:8])~hhist$mids[1:8]); abline(fit,col="blue",lty=2); 
title(main="Heads run-length frequencies on semilog scale")
#>> for exercise 11.4: compute theoretical predictions, and overlay them onto this plot 

plot(thist$mids,log(thist$density),xlab="Run length",ylab="Frequency"); 
fit=lm(log(thist$density[1:8])~thist$mids[1:8]); abline(fit,col="blue",lty=2); 
title(main="Tails run-length frequencies on semilog scale")
#>> for exercise 11.4: compute theoretical predictions, and overlay them onto this plot 

