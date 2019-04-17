rm(list=ls(all=TRUE)); # clear past actions from R's memory 

p=rep(0,5); 					        
cat("3 runs of geometric growth model","\n")
for (init in c(1,5,9)) {						
	p[1]=init; 				            
	for (j in 2:5) {					        	
	    p[j]=2*p[j-1]
	}
	cat("initial value = ",init,"\n"); 
	for (j in 1:5) {cat(j,p[j],"\n")};
}

					            
