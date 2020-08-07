Nd = 150;
ud = 0.1*rand(Nd,1); vd = 0.01*(rand(Nd,1)-0.5); wd = 0.01*(rand(Nd,1)-0.5); 
%dd(1:Nd/2,1) = (1+4*rand(Nd/2,1))*10^(-6);
%dd((Nd/2+1):Nd,1) = (5+10*rand(Nd/2,1))*10^(-6);
dd(1:130,1) = (1+9*rand(130,1))*10^(-6);
dd((131):Nd,1) = (10+10*rand(20,1))*10^(-6);
rhog = 1; 
rhod = 1000; 
gd = -9.81; 

xc = Lx/2; yc = Ly/2; zc = Lz/2; sig = ay; 

for(kkk=1:Nd)
    RRR=1000; 
    while(RRR > dnoz)
    xd(kkk) = xc + 2*(rand-0.5)*dnoz; 
    yd(kkk) = yc + 2*(rand-0.5)*dnoz;
    zd(kkk) = zc + 2*(rand-0.5)*dnoz;
    RRR = sqrt( (xd(kkk)-Lx/2).^2 + (yd(kkk)-Ly/2).^2 + (zd(kkk)-Lz/2).^2 );
    end
end
xd=xd'; yd=yd'; zd=zd';
%V = V+0.5*exp( (-(X-xc).^2 - (Y-yc).^2)/(2*sig^2)); 
%Pulse = double( abs(X-xc) < sig/2 & abs(Y-yc)<sig/2);
%U = U + Pulse*0.2; tt = 0;
%[U,V,P,dPdx,dPdy] = projectNeumann(U,V,P,OnePerK,nuP,east,west,north,south,inx,iny,dx,dy,tt,AA); 
