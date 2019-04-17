global dx2  b e;

nsteps = 5000;
nxy = 100;
u = 5*ones(nxy,1)*([1:nxy]-0.4*nxy)/nxy;
v = 5*(([1:nxy]-0.4*nxy)')*ones(1,nxy)/nxy;

dx2 = 1/2;
h  = 0.04;
e = 0.1;
b = 0.67;

colormap hsv;
set(0,'DefaultSurfaceEdgecolor','none');

