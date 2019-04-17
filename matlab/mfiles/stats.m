%[mean_x,var_x,median_x,min_x,max_x]=stats(x)
%Dynamic Models in Biology, Stephen Ellner and John Guckenheimer
%Illustrates m-files defining matlab functions with several outputs
%This function "bundles" elementary statisitics of data stored as a vector
function [mean_x,var_x,median_x,min_x,max_x]=stats(x); 
    mean_x=mean(x); var_x=var(x); median_x=median(x);
    min_x=min(x); max_x=max(x);
return;
