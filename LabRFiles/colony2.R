rm(list=ls(all=TRUE))
setwd("e:/classes/dynmod/labrfiles");
source("DMBpplane2017.R"); 

colony=function(f,h,parms){
	a=parms[1];
	o=parms[2];
	L=parms[3];
	w=parms[4];
	m=parms[5];
	i=parms[6];
	s=parms[7];
	dh=L*((h+f)/(w+h+f)) - h*(a-o*(f/(h+f))) - h*(s/(i+h+f));
	df=h*(a-o*(f/(h+f))) - m*f - f*(s/(i+h+f));
	return(c(df,dh));
}

y_label <<- "In Hive Workers"; x_label <<- "Foragers"; 

parms=c(a=0.25,o=0.75,L=2000,w=27000,m=0.24,i=0.402,s=0);

Rpplane(fun=colony, xlim=c(0,.0001),ylim=c(0,.0001),parms=parms); 
