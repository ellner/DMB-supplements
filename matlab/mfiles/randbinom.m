% Dynamic Models in Biology, Stephen Ellner and John Guckenheimer
% Vector of m binomial random deviate with parameters N and p

function b=randbinom(N,p,m);
    mu=N*p;
    if(N==0);
        b=0;
    elseif(N<1000);
        b=sum(rand(N,m)<p);
    elseif(mu<200);
        b=randpois(mu,m);
    else;
        b=round(mu+sqrt(mu*(1-p))*randn(m,1));
    end;
 