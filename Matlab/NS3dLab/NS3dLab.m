%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NS3dLab code shared under a license, see LICENSE.
% The code may or may not function on a GPU. Nevertheless, the code 
% may be executable also on newer Matlab versions on a regular CPU without 
% compatible GPU. 
% This code was used in publication to demonstrate droplet sedimentation 
% in turbulent flow and study dilution. 
% https://www.sciencedirect.com/science/article/pii/S0925753520302630
% The code is based on the pseudo-spectral solution of incompressible 
% Navier-Stokes equations with RK4 time integration utilizing the
% projection method. The model represents low speed isothermal air. 
% The fluid solver is a 3D extension of the NS2dLab
% code https://www.sciencedirect.com/science/article/pii/S0010465516300388
% A scalar transport equation is solved to mimic a passive scalar/aerosol
% concentration. 
% A lagrangian droplet tracking approach is used to simulate droplets in
% turbulence. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultAxesFontName','Times New Roman')
set(0,'DefaultTextFontName','Times New Roman')
set(0,'DefaultAxesFontSize',15)

CreateFields;
InitializeDrops;
InitializeScalar;
clear hhh

InitializeTurbulence;
CreateGpuArrays; 

% 
for(t=1:simutimeSteps)
SolveNavierStokes; 
SolveScalar;
AdvanceDrops;
%SanityCheck; % uncomment this line if you want to check if energy balance
              % works 
t
if(mod(t,1)==0)
FilterFields; % use explicit filtering as LES model with higher order filter
end

AnalyzeData;

VisuResults;
end


