
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% advection-diffusion stage with RK4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Uold = U; Vold = V; Wold = W; 
Uc   = U; Vc   = V; Wc = W; 

    
for(rk=1:4)

ConstructVelocityIncrement; 
[dU,dV,dW]=project(dU,dV,dW,KX,KY,KZ,AA,OnePerK,KXXP,KYYP,KZZP,KXYP,KXZP,KYZP); 
if(rk<4)
    U=Uold+b(rk)*dU;
    V=Vold+b(rk)*dV;
    W=Wold+b(rk)*dW;
end
    Uc=Uc+a(rk)*dU;
    Vc=Vc+a(rk)*dV;
    Wc=Wc+a(rk)*dW;
    
end

U = Uc; V = Vc; W = Wc;
