set(0,'DefaultAxesFontName','Times New Roman')
set(0,'DefaultAxesFontSize',15)
set(0,'DefaultTextFontSize',15)
 
  if(mod(k,100)==0) % increase the number to speed up
    figure(1), clf, box
    plot(xp,yp,'.')
    h = title('Walkers in a public place');
    h = xlabel('x [m]');
    h = ylabel('y [m]');
    axis([0 Lx 0 Ly]), axis equal, axis tight, drawnow, pause(0.02)

    figure(2), clf
    hold on
    imagesc([0 Lx], [0 Ly], C/1000)%, caxis([0 0.02])
          h = title('Number of aerosols per liter');
    h = xlabel('x [m]');
    h = ylabel('y [m]');
    colorbar 
    axis([0 Lx 0 Ly]), axis equal, axis tight, drawnow, pause(0.1)

  end

%  PDF = PDF + hist(C(:), 0:0.01:lambdacough); % bin centers, max  =1000 /m^3
  
  if(mod(k,HowOftenSample)==0)
%      figure(3), clf
%      loglog([0:0.01:lambdacough]/1000, PDF)
%    h = title('PDF of number concentration');
%    h = xlabel('Number of aerosols per liter');
%    drawnow, pause(0.02)

    for(s=1:Np)
        dosemean = mean(dose( sp==0 ));
    end
    
    aveinf(1+round(k/HowOftenSample))= dosemean; 
    maxinf(1+round(k/HowOftenSample))= max(dose);
    aveC(1+round(k/HowOftenSample))  = mean(mean(C));
    
   if(mod(k,400000)==0)
    figure(4), box, hold on
        plot(k*dt, dosemean,'ko')
        plot(k*dt, max(dose),'kx')
%         h = title('PDF of number concentration');
        h = xlabel('Time [s]');
         h = ylabel('Total number of inhaled aerosols');
        h = legend('Average dose', 'Maximum dose','Location','NorthWest');
         
        drawnow, pause(0.1)
   end
   end