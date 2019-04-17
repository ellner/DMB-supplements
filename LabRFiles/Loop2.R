rm(list=ls(all=TRUE)); # clear past actions from R's memory 

p=rep(0,5); 					        
for (init in c(1,5,9)) {						
	p[1]=init; 				            
	for (j in 2:5) {					        	
	    p[j]=2*p[j-1]
	    cat(init,j,p[j],"\n"); 	  
        }
}						            
