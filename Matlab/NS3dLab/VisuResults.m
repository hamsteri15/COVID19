% update this file on demand
if(mod(t,400)==0)
    Ug = gather(U); 
figure(1), clf,    hold on
    imagesc([0,Lx],[0,Ly],Ug(:,:,Nz/2)), colormap copper
plot(xd(1:Ndactive),yd(1:Ndactive),'ko', 'MarkerFaceColor','k'), axis([0,gather(Lx),0,gather(Ly)]); axis equal, axis tight, drawnow  
hhh=text(0.1,2.8,strcat('Time = ', num2str(round(t*dt)), ' sec.')); 
set(hhh,'Fontsize',24)
hhh=xlabel('x [m]'); set(hhh,'Fontsize',24)
hhh=ylabel('y [m]'); set(hhh,'Fontsize',24)
drawnow, pause(0.1)
end



if(mod(t,10)==0 & t > 1)
    Tg = gather(T);
figure(3), clf,    hold on
axis([0,gather(Lx),0,gather(Ly)])
%imagesc([0,Lx],[0,Ly],1-Tg(:,:,Nz/2)), colormap bone
imagesc([0,Lx],[0,Ly],exp(-sum(T,3)/10)), colormap bone
%imagesc([0,Lx],[0,Ly],1-Tg(:,:,Nz/2)), colormap bone
%Visu3d; 

caxis([0, 1])
plot(xd(1:130),yd(1:130),'b.', 'MarkerFaceColor','b'); 
plot(xd(131:150),yd(131:150),'ro', 'MarkerFaceColor','r'); box
hhh=text(0.1,2.8,strcat('Time = ', num2str(round(t*dt)), ' sec.')); 
set(hhh,'Fontsize',24)
hhh=xlabel('x [m]'); set(hhh,'Fontsize',24)
hhh=ylabel('y [m]'); set(hhh,'Fontsize',24)
h=legend('d \leq 10 micron', '10 micron < d \leq 20 micron'); set(h,'Fontsize',24)
drawnow
eval( strcat('print -dpng ', ' Scalar', num2str(1000+t))); 
drawnow, pause(0.1)
figure(4), clf
hist(T(find(T>0.01)),100)

end

