%FILE: randomwalk2.m test
% Matrix size m x m, m eg 100
% Location matrix a   m x m
% a(jj,kk) is the number of visits in location (jj,kk), 1<=jj,kk<=m
% Number of walkers : p
% Number of time steps: nsteps
% from location (x,y) each walker moves to
%    (x+dx, y+dy) ; dx= rand integer in (-2,-1,0,1,2)
%    (x+dx, y+dy) ; dy= rand integer in (-2,-1,0,1,2)
% If x+dx is <1 or >m move in opposite direction
% If y+dy is <1 or >m move in opposite direction

close all
clear
m=50;
a=zeros(m,m);
nsteps=50000;
p=100;
myloc=randi(m,[p,2]);
for k1=1:nsteps
    for k2=1: p
        %location of walker k2
        dx=randi(5)-3; dy=randi(5)-3;
        x=myloc(k2,1); y=myloc(k2,2);
        xtmp=x+dx; ytmp=y+dy;
        if (xtmp<1)
            xtmp=x+abs(dx);
        end
         if ( xtmp>m)
            xtmp=x-abs(dx);
        end
        if (ytmp<1)
            ytmp=y+abs(dy);
        end

         if  (ytmp>m)
            ytmp=y-abs(dy);
        end
        myloc(k2,1)=xtmp; myloc(k2,2)=ytmp;
        %x=randi(m); y=randi(m);
        a(xtmp, ytmp)= a(xtmp, ytmp)+1;
    end
end
% Also bar3(a) gives picture
surf(1:m,1:m,a)
