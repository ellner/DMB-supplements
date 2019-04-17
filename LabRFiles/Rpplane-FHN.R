
fhn=function(x,y,parms) {
	dx = parms["C"]*(x-(1/3)*x^3-y+parms["j"]);
	dy = (1/parms["C"])*(x+parms["a"]-parms["b"]*y)
	return(c(dx,dy))
}


#parms1 = c(a = -0.5, b = 0.5, C = 1, j = -0.5)
#Rpplane(fun=fhn,c(-3,3),c(-2,2),parms=parms1)

parms2 = c(a = -0.5, b = 1.25, C = 1.95, j = -0.45)
Rpplane(fun=fhn,c(-3,3),c(-2,2),parms=parms2)