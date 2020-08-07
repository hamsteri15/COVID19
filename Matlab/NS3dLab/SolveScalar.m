
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% advection-diffusion stage with RK4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Told = T; 
Tc   = T; 
    
for(rk=1:4)

ConstructScalarIncrement; 

if(rk<4)
    T=Told+b(rk)*dT;
end
    Tc=Tc+a(rk)*dT;   
end

T = Tc; 