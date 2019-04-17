function F=FCDF(x,nn,nd)
F=betainc(nn.*x./(nd+nn.*x),nn/2,nd/2);