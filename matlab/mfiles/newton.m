function [x,df] = newton(f,x0,p)
  inc = 1e-7; %increment for computing derivatives
  conv = 1e-14 % Convergence criterion
  niter = 40 % Maximum number of iterations
  dim = length(x0); %dimension
  t = 0; % arbitrary time value when calling f
  x = x0; %initial x
  for i=1:40 %Newton iteration loop
    y = feval(f,t,x,p) %Evaluate f
    df = zeros(dim); % Compute df
    for j=1:dim
	v=zeros(dim,1); 
	v(j) = inc; %Increment vector for estimating derivative wrt jth coordinate
        df(:,j) = (feval(f,t,x+v,p) - y)/inc;
    end;
    if (norm(y) < 1e-14) %Test for convergence
        return; %Return if converged
    end;
    x = x - df\y %New value of x
  end;
  display('Convergence failed');
