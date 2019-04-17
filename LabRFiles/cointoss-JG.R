rm(list=ls(all=TRUE)); # clear past actions from R's memory 
graphics.off(); 

nt = 5000000;
p = 0.6;      # probability of heads
r = runif(nt+1); 
rd = (r < p); # tails = 0, heads = 1
rdiff = rd[1:nt] + 2*rd[2:(nt+1)];         # binary code for transitions: tt = 0, ht = 1, th = 2, hh = 3
ht = which(rdiff == 1);                    # heads to tails indices
th = which(rdiff == 2);                    # tails to heads indices
rn  = min(length(ht),length(th));          # numbers of transitions

###Compute residence times from ht and th transition times 
if (rd[1]==0) {                             # if tails first
   rh = ht[1:rn] - th[1:rn];                  # heads residence times
   rt = c(th[1],th[2:rn] - ht[1:rn-1]);       # tails residence times
}else{                                      # if heads first
   rh = c(ht[1],ht[2:rn] - th[1:rn-1]);       # heads residence times
   rt = th[1:rn] - ht[1:rn];                  # tails residence times
}

### Plot results 
par(mfrow=c(2,1),cex.axis=1.35,cex.lab=1.35); 
hhist = hist(rh,(0:max(rh))+0.5,xlim=c(1,max(rh)));     # heads histogram
thist = hist(rt,(0:max(rt)+0.5),xlim=c(1,max(rt)));     # tails histogram

### Plot results on log scale 
### If run-lengths have a geometric distribution, this should give
### a straight-line plot. 
### Do ?hist  and read the Value section to see what's being plotted

dev.new(); 
par(mfrow=c(2,1),cex.axis=1.35,cex.lab=1.35);  
plot(hhist$mids,log(hhist$density),xlab="Run length",ylab="Frequency"); 
fit=lm(log(hhist$density[1:8])~hhist$mids[1:8]); abline(fit,col="blue",lty=2); 
title(main="Heads run-length frequencies on semilog scale")

plot(thist$mids,log(thist$density),xlab="Run length",ylab="Frequency"); 
fit=lm(log(thist$density[1:8])~thist$mids[1:8]); abline(fit,col="blue",lty=2); 
title(main="Tails run-length frequencies on semilog scale")

### Now do Exercise 11.4: compare the simulation results with the 
### exact theoretical prediction based on the value of p
