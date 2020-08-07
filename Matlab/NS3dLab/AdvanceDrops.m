subiter = 3; % do 3 subiterations 
dtd = dt/subiter; % droplet timestep 

% index labels of droplets
 indy = min(Ny,max(1,1+round((Ny-1)*yd/Ly))); indx = min(Nz,max(1, 1+round((Nx-1)*xd/Lx)));
 indz = min(Nz,max(1, 1+round((Nz-1)*zd/Lz)));

% gas velocity at droplet positions
 ug = U(sub2ind(size(U),indy,indx,indz)); 
 vg = V(sub2ind(size(V),indy,indx,indz)); 
 wg = W(sub2ind(size(W),indy,indx,indz)); 

for(kk=1:subiter)
SolveDrops; 
end

% create a treatment what happens to droplets at domain boundaries
if(1)
for(kkk=1:Nd)
    if(yd(kkk)>Ly)
        yd(kkk)=Ly; ud(kkk)=0; vd(kkk)=0; wd(kkk)=0;
    elseif(yd(kkk)<0)
        yd(kkk)=0; ud(kkk)=0; vd(kkk)=0; wd(kkk)=0;
    end
    
    if(xd(kkk)>Lx)
        xd(kkk)=0; ud(kkk)=0; vd(kkk)=0; wd(kkk)=0;
    elseif(xd(kkk)<0)
        xd(kkk)=Lx; ud(kkk)=0; vd(kkk)=0; wd(kkk)=0;
    end
    
    if(zd(kkk)>Lz)
        zd(kkk)=0; ud(kkk)=0; vd(kkk)=0; wd(kkk)=0;
    elseif(zd(kkk)<0)
        zd(kkk)=Lz; ud(kkk)=0; vd(kkk)=0; wd(kkk)=0;
    end
    
end
end