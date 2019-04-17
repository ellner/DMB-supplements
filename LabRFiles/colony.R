rm(list=ls(all=TRUE))
setwd("C:/Users/jb2294/Downloads");
source("DMBpplane2017.R"); 

colony=function(F,H,parms){
	a=parms[1];
	o=parms[2];
	L=parms[3];
	w=parms[4];
	m=parms[5];
	i=parms[6];
	s=parms[7];
	dh=L*((H+F)/(w+H+F)) - H*(a-o*(F/(H+F))) - H*(s/(i+H+F));
	df=H*(a-o*(F/(H+F))) - m*F - F*(s/(i+H+F));
	return(c(dh,df));
}

y_label <- "In Hive Workers"; x_label <- "Foragers"; 

parms=c(a=0.25,o=0.75,L=2000,w=27000,m=0.24,i=0.402,s=0);
Rpplane(fun=colony, xlim=c(0,15000),ylim=c(0,6000),parms=parms); 
