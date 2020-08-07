Nx          = 128; 
Ny          = 128; 
Nz          = 128;
dt          = 0.06*3; % timestep
Lx          = pi;    % domain side length
Ly          = pi;    % domain side length
Lz          = pi;    % domain side length

nu = 1.6e-5; counter = 1; HowOftenSave = 5;


dx = Lx/Nx;     % grid spacing dx
dy = Ly/Ny;     % grid spacing dy
dz = Lz/Nz;

simutimeSeconds = 240;
simutimeSteps   = round(simutimeSeconds/dt); % how many timesteps altogether
visualize       = 1; % set 0 for speed-up   

%%%%%%%%%%%%%%%%%%%%
% create the 2D grid
%%%%%%%%%%%%%%%%%%%%
x = (Lx-dx)*(0:(Nx-1))/(Nx-1); % to be precise, dx is subtracted here 
y = (Ly-dy)*(0:(Ny-1))/(Ny-1); % otherwise the derivative of
                               % e.g. sin(x) wouldnt be
                               % differentiable on this periodic
                               % grid
z = (Lz-dy)*(0:(Nz-1))/(Nz-1); % otherwise the derivative of
                               % e.g. sin(x) wouldnt
                               %be differentiable on this periodic grid
[X,Y,Z] = meshgrid(x,y,z);  