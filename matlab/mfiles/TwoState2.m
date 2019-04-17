% Dynamic Models in Biology, Stephen Ellner and John Guckenheimer
% TwoState2.m simulates a group of particles moving 
% back and forth between two compartments, using Poisson and Normal
% approximations for the coin-tossing process (randbinom.m). 
% Also calls diffus2.m.

% parameters
R=0.3; L=0.2; dt=0.01; tmax=25; 

% initial condition
nL0=500; nR0=0; ntot=nL0+nR0; 

% movement probabilities per time step 
Rdt=R*dt; Ldt=L*dt;

% initialize variables for current population size
nL=nL0; nR=nR0; 

% initialize vectors storing the results over time
ntL=nL; ntR=nR; tvals=0;

t=0; 
while t<tmax;
   % use Poisson to do the coin-tossing 
   	  moveR=randbinom(nL,Rdt,1);
	  moveL=randbinom(nR,Ldt,1);
   % make the changes
     nL=nL + moveL - moveR;
     nR=nR + moveR - moveL;
   % store the results for plotting
   	t=t+dt; tvals=[tvals;t];
	ntL=[ntL;nL]; ntR=[ntR;nR];
end; 

% solve the ODEs for the continuous time chain (dt --> 0)
tspan=[0 tmax]; y0=[nL0;nR0]; 
[t2,y2]=ode45('diffus2',tspan,y0,[],L,R);

% plot stochastic simulation and deterministic solutions
plot(tvals,ntL/ntot,tvals,ntR/ntot); axis([0 tmax 0 1])
hold on; plot(t2,y2/ntot); hold off;	 

% label the plot
xlabel('Time','Fontsize',14); ylabel('Fraction of Particles','Fontsize',14); 
title('Two-compartment diffusion','Fontsize',16)
legend('Left','Right');
