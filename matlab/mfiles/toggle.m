function dy =  toggle(t,y,p)

    dy = zeros(2,1);
    dy(1) = - y(1) + p(1)./(1+y(2).^p(2));
    dy(2) = - y(2) + p(1)./(1+y(1).^p(3));
 
