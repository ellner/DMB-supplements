%Dynamic Models in Biology, Stephen Ellner and John Guckenheimer
%Second example of Matlab loops: nested loops

%Initialize array to hold data
p=zeros(5,1);

%Set initial population sizes: three values
for init=linspace(1,10,3);      		
	p(1)=init; 	

    %Loop as in previous example    
	for n=2:5;					   	
	    p(n)=p(n-1)*2; x=log(p(n));		
		q=[num2str(n), '  ', num2str(x)];	
	    disp(q)				   
	end;	
    
    %Print blank line to separate data
    disp(' ')
end;						         
