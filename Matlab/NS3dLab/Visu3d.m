[XX,YY,ZZ] = sphere(120); 

XX=XX;
YY=YY;
ZZ=ZZ;


%figure(10), clf, hold on
f=2000;
for(k=1:Nd)
    surf(XX.*min(dd(k),30*10^-6)*f+xd(k),YY.*min(dd(k),30*10^-6)*f+yd(k),ZZ.*min(dd(k),30*10^-6)*f+zd(k))
end


view(0,90)
shading interp
lightangle(-45,-45)
h.FaceLighting = 'phong';
h.AmbientStrength = 0.3;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';

colormap bone
