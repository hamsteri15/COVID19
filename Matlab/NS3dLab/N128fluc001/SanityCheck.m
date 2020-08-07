% This file can be used to check if the code works properly
% See the DNSLab paper by Vuorinen et al.
% Kinetic energy decay rate should approximately match
% dEkin/dt = -dissipation_phys - dissipation_num
% For well resolved LES dEkin/dt \approx -dissipation_phys   

if(t==1)
    sss=0;
end

if(mod(t,1)==0)
    sss=sss+1; 
figure(13), hold on 
dUdx = real(ifftn(iKX.*fftn(U))); 
dUdy = real(ifftn(iKY.*fftn(U)));
dUdz = real(ifftn(iKZ.*fftn(U)));
dVdx = real(ifftn(iKX.*fftn(V))); 
dVdy = real(ifftn(iKY.*fftn(V)));
dVdz = real(ifftn(iKZ.*fftn(V)));
dWdx = real(ifftn(iKX.*fftn(W))); 
dWdy = real(ifftn(iKY.*fftn(W)));
dWdz = real(ifftn(iKZ.*fftn(W)));
Ekin(sss) = 0.5*sum(sum(sum( U.^2 + V.^2 + W.^2 )));
ksi = sum(sum(sum(nu*(dUdx.^2 + dUdy.^2 + dUdz.^2 + dVdx.^2 + dVdy.^2 + dVdz.^2 +...
    dWdx.^2 + dWdy.^2 + dWdz.^2)))); 
if(sss>2)
%plot(t, sum(sum(sum( U.^2 + V.^2 + W.^2 ))),'ko')
hold on
plot(t, -(Ekin(sss)-Ekin(sss-1))/dt,'x')
plot(t, ksi,'ro')

end
end
