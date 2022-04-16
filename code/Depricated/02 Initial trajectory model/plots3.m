% function [] = plots3(C,C0,P,P0,B,tau,tfapp,amax,r0,mu,m,animations, n)

%% prelinimaries


% for gif image
imsize = 350;
im = cell(1,m);
xpos=10; ypos=150;

% time calculations
tF = flight_time(P, B, m, tfapp, r0, n);
time_days = tau.*(tF/(60*60*24));


% convert cylindrical coordinates to cartesian
[x, y, z] = cylindrical_to_cartesian(C);


%% Orbit representation


figure_orbits = figure(fig);
set(figure_orbits,'position',[xpos,ypos,1.2*imsize,imsize])

hold on
xlabel('x coordinate [AU]')
ylabel('y coordinate [AU]')
plot(0,0,'*k');
plot(x(1),y(1),'*k');
viscircles([0,0],initial.pos(1)/r0,'LineStyle','--','Color','r','LineWidth',0.3); % earth's orbit
viscircles([0,0],final.pos(1)/r0,'LineStyle','-.','Color','b','LineWidth',0.3); % mars' orbit
hold on

if animations==1
    for i=1:m
        marker = plot(x(i),y(i),'k*');
        hold on
%         legend(num2str(i)./m)
        title(strcat('time elapsed = ',{' '}, num2str(round(time_days(i))),' days'));
        frame = getframe(figure_orbits);
        im{i} = frame2im(frame);
        delete(marker)
        
        plot(x(i),y(i),'k.');
    end
    
    frame = getframe(figure_orbits);
    im{i} = frame2im(frame);
    writegif('orbit.gif',im,m,2/m);
end


legend('off')
title('orbit representation')
plot(x,y,'k','LineWidth',1);
plot(x(m),y(m),'*k');



%% plot propulsive acceleration

figure_propulsion = figure(fig+1);
set(figure_propulsion,'position',[xpos + 1.2*imsize,ypos,1.2*imsize,imsize])
title('propulsive acceleration')
hold on
[arho, atheta, az] = aplot(C,r0,tfapp,mu);
plot(time_days, acceleration(P, mu, B, r0, tfapp,n), 'k','LineWidth',1)
plot(time_days,arho, 'LineWidth',0.3)
plot(time_days,atheta, 'LineWidth',0.3)
plot(time_days,az, 'LineWidth',0.3)
yline(amax, '--k')
xlabel('flight time [days]')
ylabel('prop (final) [m/s^2]')
legend('a total','a \rho','a \theta','a z')


%% plot evolution of coordiates
figure_coordinates = figure(fig+2);
set(figure_coordinates,'position',[xpos + 2*1.2*imsize,ypos,1.2*imsize,imsize])
title('coordinate evolution')
hold on

% earth orbit
thetaE = 2*pi/365.25*time_days+initial.pos(2);
xE = initial.pos(1)/r0.*cos(thetaE);
yE = initial.pos(1)/r0.*sin(thetaE);

% mars orbit
% thetaM = 2*pi/687*(time_days-47);
thetaM = 2*pi/687*time_days+final.pos(2)-pi/2-0.075;
xM = final.pos(1)/r0.*cos(thetaM);
yM = final.pos(1)/r0.*sin(thetaM);

subplot(2,1,1)
hold on
plot(time_days, xM, 'LineStyle','-.','Color','b','LineWidth',0.3)
plot(time_days, xE, 'LineStyle','--','Color','r','LineWidth',0.3)
plot(time_days, x, 'k','LineWidth',1)
plot(time_days(1),x(1),'*k','DisplayName','')
plot(time_days(end),x(end),'*k','DisplayName','')
xlabel('flight time [days]')
ylabel('x [au]')
% legend('earth','mars','probe')

subplot(2,1,2)
hold on
plot(time_days, yM, '-.','LineWidth',0.3)
plot(time_days, yE, '--','LineWidth',0.3)
plot(time_days, y, 'k','LineWidth',1)
plot(time_days(1),y(1),'*k','DisplayName','')
plot(time_days(end),y(end),'*k','DisplayName','')
xlabel('flight time [days]')
ylabel('y [au]')



%%

function [arho, atheta, az] = aplot(C,r0,tfapp,mu)



[rho, theta, z] = extract_coordinates(C, r0, tfapp);


% find distance to sun
r = sqrt(rho.o.^2 + z.o.^2);

% equations of motion
arho = rho.DD - rho.o.*theta.D.^2 + mu.*rho.o./r.^3;
atheta = rho.o.*theta.DD + 2.*rho.D.*theta.D;
az = z.DD + mu.*z.o./r.^3;

end