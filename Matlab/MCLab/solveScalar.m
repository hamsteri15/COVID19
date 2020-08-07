% diffusion equation with source and sink terms for aerosol concentration 
% following the paper using expl.Euler and CD2 method (finite difference)
C = C + dt*D*(C(north,inx)+C(south,inx)+C(iny,east)+C(iny,west)-4*C(iny,inx))/dx^2 ;

C = C + S - dt*C/tau; % [C] number of aerosols in the control volume
