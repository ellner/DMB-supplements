% Dynamic Models in Biology, Stephen Ellner and John Guckenheimer
% TwoState.m simulates a group of particles all moving 
% back and forth between two compartments, using a coin-tossing
% process. In time interval dt, each particle in the left box has
% probability R*dt of moving to the right box, and each one in the right
% box has probability L*dt of moving into the left box. For dt small
% and a large number of particles, the process should approximate
% a two-compartment differential equation
%
%   dQ1/dt= -R*Q1 + L*Q2      dQ2/dt= -L*Q2 + R*Q11
%
% Simulations of the Markov Chain are plotted along with solutions
% of the differential equation defined in the file diffus2.m

% parameters
R=0.3; L=0.2; dt=0.1; tmax=25; 

% initial condition
nL0=2500; nR0=0; ntot=nL0+nR0; 

% movement probabilities per time step 
Rdt=R*dt; Ldt=L*dt;

% initialize variables for current population size
nL=nL0; nR=nR0; 

% initialize vectors storing the results over time
ntL=nL; ntR=nR;  t=0; tvals=0;

while t<tmax;
   % how many move from L to R?
    if nL>0
	  moveR=sum(rand(nL,1)<Rdt);
	else
	  moveR=0;
	end
% how many move from R to L?
	if nR>0	
	  moveL=sum(rand(nR,1)<Ldt);
	else
	  moveL=0;
	end	
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
plot(tvals,ntL,tvals,ntR); axis([0 tmax 0 ntot])
hold on; plot(t2,y2); hold off;	 

% label the plot
xlabel('Time','Fontsize',14); ylabel('Particles','Fontsize',14); 
title('Two-compartment diffusion','Fontsize',16)
legend('Left','Right');

