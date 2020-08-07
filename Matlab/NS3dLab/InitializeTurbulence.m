UseStart = 1; % the parameter value = 1 -> reads initial turbulence field Ustart.mat 

if(UseStart==1)
load Ustart.mat; 

Resize; % scales the initial data to another grid resolution

% initialize turbulence level  
TKE = mean(mean(mean( 0.5*(U-mean(mean(mean(U)))).^2 + 0.5*(V-mean(mean(mean(V)))).^2 + 0.5*(W-mean(mean(mean(W)))).^2)));
Upp = sqrt(TKE); % first for initial data
U = U*(fluc/Upp); V = V*(fluc/Upp); W = W*(fluc/Upp); % rescale to fluc 
U = U-mean(mean(mean(U))); V = V-mean(mean(mean(V))); 
W = W-mean(mean(mean(W)));
[U,V,W] = project(U,V,W,KX,KY,KZ,AA,OnePerK,KXXP,KYYP,KZZP,KXYP,KXZP,KYZP);

end