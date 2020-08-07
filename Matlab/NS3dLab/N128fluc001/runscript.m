addpath /m/home/home9/94/vavuorin/unix/NS3dLab/
clear *
!rm *.png
fluc = 0.01; % initial turbulence level (m/s)
SetParameters; 

% call the solver
NS3dLab;

% save different parameters to a struct to be post processed later on
% the dilution curve can be calculated from the HIST i.e. histogram of 
% scalar concentration
resu.HIST = gather(HIST); resu.dt = gather(dt); resu.HowOftenSave = gather(HowOftenSave);
resu.fluc = gather(fluc); resu.binvals = binvals; 
 
save('result.mat', '-struct', 'resu');