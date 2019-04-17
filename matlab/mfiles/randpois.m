% RANDPOIS(p,n)
%   This function creates a random variable with a Poisson distribution
%   with parameter 'a'.  If a second argument is used, a vector of
%   'n' poisson random variables is created.
%
% See also RAND, RANDN, RANDUNIFC, RANDEXPO, RANDGEO, RANDGAUSS

function out = randpois(a,n)

if nargin == 1
    i = 0;
    f = exp(-a);
    p = exp(-a);
    u = rand;
    while f<=u
        p = p * (a / (i + 1));
        f = f + p;
        i = i + 1;
    end
    out = i;
end

if nargin == 2
    for k=1:n
        i = 0;
        f = exp(-a);
        p = exp(-a);
        u = rand;
        while f<=u
            p = p * (a / (i + 1));
            f = f + p;
            i = i + 1;
        end
        out(k) = i;
    end
end