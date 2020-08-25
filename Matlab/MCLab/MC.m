%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NS3dLab code shared under a license, see LICENSE.
% Monte-Carlo code shared under the MIT license. Simulates 
% random walkers infecting others in a public place
% The code accumulates dose for each individual
% This code was used in publication to demonstrate virus exposure process
% in Matlab as described in 
% https://www.sciencedirect.com/science/article/pii/S0925753520302630
% Various parameter sweeps can be easily carried out as mentioned below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear *
close all
counter = 0;

Nps = [1000 4000]; % Number of persons e.g. sweeps over 1000 4000
sicks = [0.01];    % fraction of infected
Pcoughs = [0 1/3600 2/3600];     % coughing or not e.g. none once/hour twice/hour
superswitch = [0]; % possibility to investigate superemitters 
Upaves = [0.01 0.1 0.5 1]; % walking speed 
taus = [100];      % ventilation timescale in seconds 
Ds = [0.05];       % diffusivity including aerosol mixing due to turbulence

!rm Results*.mat
HowOftenSample = 25;

% sweep over various parameters

for(aaa=1:length(Nps))
    for(bbb=1:length(sicks))
    for(ccc=1:length(Pcoughs))
        for(ddd=1:length(superswitch))
            for(eee=3:length(Upaves))
                for(fff=1:length(taus))
                    for(ggg=1:length(Ds))
%%%%%%%%%%%%%%%%%%%%%%%%%%
% System dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%
Ny=100; Nx=100; Np = Nps(aaa); 
Lx=100; Ly=100; % meters 100m*100m

% reference tables for finite difference diffusion simulation
iny = 1:Ny; inx = 1:Nx; 
north = iny+1; north(Ny)=Ny; 
south = iny-1; south(1) = 1;
east = inx+1; east(Nx)=Nx; 
west = inx-1; west(1)=1;
dx = Lx/Nx; dy = Ly/Ny;

C = zeros(Ny,Nx); % concentration, assume zerogradient bc's

S = zeros(Ny,Nx); % source term in units [# of particles / m^3 ]
                  % assume control volume = 1m^3  

simuTime = 3600;  % seconds
dt = 0.1;         % seconds
simuSteps = round(simuTime/dt); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters of sick individuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = 5; % lambda aerosols [1/s] e.g. 5 means 5 per second into control volume
sick = sicks(bbb); % fraction of sick persons
nofsick = round(sick*Np); 
superemitter = 0.1; % you can pick the fraction how many emit more than others 
nofsuper = superswitch(ddd)*max(1, round(superemitter*nofsick));

% x,y coordinates of Np persons, sp = 0 -> healthy sp=1 -> sick
xp = Lx*rand(Np,1); yp = Ly*rand(Np,1); sp = zeros(Np,1); sp(1:nofsick) = 1; se = zeros(Np,1); se(1:nofsuper)=1;
Upave = Upaves(eee); % velocity of each individual 
tau = taus(fff); % timescale [s] of removing aerosols upwards due to ventilation
D = Ds(ggg); % m^2/s  diffusivity in x,y plane direction due to turbulence

lambdacough = 40000; % how many particles come from a single cough, assume this is 
                    % immediately spread to the control volume during dt
                    
Vinh = 0.001*20/60;   % how many m^3 you inhale in a second, here 20 l per 60 sec
Vvol = dx*dy*1; % volume in which the droplets locate; now assume 1m^3 if dx=dy=1m
Pcough = Pcoughs(ccc); % Here. e.g. 1 cough per 3600 sec; then during dt probability dt*Pcough

xt = rand(Np,1)*Lx; yt = rand(Np,1)*Ly; % target points where walker is going 
dose = zeros(Np,1);                     % the healthy persons gain a dose 
                                        % i.e. number of aerosols
PDF = zeros(1,1+lambdacough*100);                     % PDF that accumulates the dose 

tic
% run single simulation with fixed parameters 
for(k=1:simuSteps) 
     solveSource; 
     solveScalar; 
     solvePersons;
   if(mod(k,100)==0)
    VisuResults;
    k
   end
 
  
end  
toc

counter = counter + 1
SaveResults;    
                    end
    end
end
    end
end
    end
end



