global dx2  b e;

nsteps = 500;

nxy = 100;

v = 5*(([1:nxy]-0.4*nxy)')*ones(1,2*nxy)/nxy;
u = 5*ones(nxy,1)*([1:nxy]-0.4*nxy)/nxy;
u = [u,fliplr(u)] % + [1:nxy]'*[1:2*nxy]/nxy^2 -1;


dx2 = 1/9;
h  = 0.04;
e = 0.1;
b = 0.67;

colormap jet;
set(0,'DefaultSurfaceEdgecolor','none');

