% distance dxp and dyp to target destination
dxp = (xt-xp);
dyp = (yt-yp);
dr = sqrt(dxp.^2 + dyp.^2);  

up = Upave*dxp./max(1e-9,dr); vp = Upave*dyp./max(1e-9,dr); % velocity components

% update position of each person
xp = xp + up*dt; 
yp = yp + vp*dt; 

reached=double( dr < 1 ); % reaches if we are within 1 m from destination 
                          % binary 1/0  

% target points
xt = reached.*rand(Np,1)*Lx+(1-reached).*xt;
yt = reached.*rand(Np,1)*Ly+(1-reached).*yt;
