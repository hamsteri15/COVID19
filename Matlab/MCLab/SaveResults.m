% you need to modify the 
% save(strcat('ResultsMay14/Results_',num2str(aaa),'_',num2str(eee),'.mat'), '-struct', 'Resu');
% line so this line e.g. now only assumes that the parameter sweeps over 
% aaa and eee is carried out. If the for loops are over aaa bbb ccc .. 
% -> modify accordingly

Resu.C = C; Resu.maxinf = maxinf; Resu.aveinf = aveinf; Resu.aveC = aveC;
Resu.HowOftenSample = HowOftenSample;
Resu.dt = dt; Resu.Pcough = Pcough;
Resu.Np = Np; Resu.sick = sick; Resu.dt = dt; Resu.Lx=Lx; Resu.Ly=Ly;
Resu.Upave = Upave; Resu.tau = tau; Resu.D = D; 
Resu.xp = xp; Resu.yp = yp; Resu.sp = sp; Resu.dose = dose; 
save(strcat('ResultsMay14/Results_',num2str(aaa),'_',num2str(eee),'.mat'), '-struct', 'Resu');
