if(Nx ~= 128)
    for(k=1:128)
    Utmp=imresize(s1.a(:,:,k),2); 
    Vtmp=imresize(s1.b(:,:,k),2);
    Wtmp=imresize(s1.c(:,:,k),2);
    U(:,:, 2*(k-1)+1) = Utmp; 
    V(:,:, 2*(k-1)+1) = Vtmp; 
    W(:,:, 2*(k-1)+1) = Wtmp; 
    U(:,:, 2*(k-1)+2) = Utmp; 
    V(:,:, 2*(k-1)+2) = Vtmp; 
    W(:,:, 2*(k-1)+2) = Wtmp; 
    end
else  
    U=s1.a; V=s1.b; W=s1.c; 
end