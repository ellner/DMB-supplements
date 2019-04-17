function [uf,vf] = sfn(u,v)
global dx2 e b ;

ny = size(u,1);
nx = size(u,2);

uer = [u(:,1),u,u(:,nx)];
uec = [u(1,:);u;u(ny,:)];
ul = uec(3:ny+2,:)+uec(1:ny,:)+uer(:,1:nx)+uer(:,3:nx+2)-4*u;
     u3 = u.*u.*u;

uf = (u-u3/3-v)/e + dx2*ul;
vf = e*(u+b-0.5*v);
