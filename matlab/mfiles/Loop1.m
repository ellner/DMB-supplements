%Dynamic Models in Biology, Stephen Ellner and John Guckenheimer
%First example of Matlab loops

% initial population size
initsize=4; 

% create matrix to hold results sizes and store initial size 
popsize=zeros(10,1); popsize(1)=initsize;

% calculate population size at times 2 through 10, write to Command Window
for n=2:10;  
	popsize(n)=popsize(n-1)*2;
	x=log(popsize(n));
    q=[num2str(n), '  ', num2str(x)];
end;
