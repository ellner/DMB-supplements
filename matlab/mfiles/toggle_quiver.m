[U,V] = meshgrid(0:.2:3);
Xq = -U + 3./(1+V.^2);
Yq = -V + 3./(1+U.^2);
quiver(U,V,Xq,Yq);
