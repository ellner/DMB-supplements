for i=1:nsteps
[uf,vf] = sfn(u,v);
u = u+h*uf;
v = v+h*vf;
if mod(i,5) == 1
pcolor(u); drawnow;
end 
end 
