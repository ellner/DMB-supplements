function yprime=diffus2(t,y,options,L,R);
yprime= [-R*y(1)+L*y(2) ; R*y(1)-L*y(2) ];
	 
