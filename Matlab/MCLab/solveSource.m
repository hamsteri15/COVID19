% zero source term
S = 0*S;

% make person coordinate an integer 
xpint = max(min(round(Nx*xp/Lx),Nx),1); 
ypint = max(min(round(Ny*yp/Ly),Ny),1);
    
for(s=1:Np)

    % sick persons generate lambda*dt particles per volume during dt
    if(sp(s)==1)    
    S(xpint(s),ypint(s))=S(xpint(s),ypint(s))+lambda*dt + se(s)*lambda*9*dt;
    % coughing occurs at probablity Pcough*dt
    if(rand < Pcough*dt)
        S(xpint(s),ypint(s))=S(xpint(s),ypint(s))+lambdacough;
    end
    % healthy individuals accumulate a dose 
    else
    dose(s) = dose(s)+C(xpint(s), ypint(s))*(Vinh/Vvol)*dt;    
    end
end

