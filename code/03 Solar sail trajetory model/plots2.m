function [] = plots2(C,C0,P,P0,B,tau,tfapp,amax,r0,mu,m,animations, n)


close all


% for gif image
imsize = 400;
im = cell(1,m);

% time calculations
tF = flight_time(P, B, m, tfapp, r0, n);
time = tau.*tF;
time_days = time./(60*60*24);


%% convert cylindrical coordinates to cartesian
[x, y, z] = cylindrical_to_cartesian(C);
[x0, y0, z0] = cylindrical_to_cartesian(C0);

%% plot propulsive acceleration

figure_propulsion = figure(1); xpos=10; ypos=200;
set(figure_propulsion,'position',[xpos,ypos,1.2*imsize,imsize])

subplot(2,1,1)
hold on
plot(time_days, acceleration(P0, mu, B, r0, tfapp,n))
yline(amax, '--k')
xlabel('flight time [days]')
ylabel('prop (initial) [m/s^2]')

subplot(2,1,2)
hold on
plot(time_days, acceleration(P, mu, B, r0, tfapp,n))
yline(amax, '--k')
xlabel('flight time [days]')
ylabel('prop (final) [m/s^2]')



%% plot evolution of coordiates
figure_coordinates = figure(2);
set(figure_coordinates,'position',[xpos + 1.2*imsize,ypos,1.2*imsize,imsize])
hold on

subplot(3,1,1)
hold on
plot(time_days, x0)
plot(time_days, x)
xlabel('flight time [days]')
ylabel('x [au]')
legend('initial','final')

subplot(3,1,2)
hold on
plot(time_days, y0)
plot(time_days, y)
xlabel('flight time [days]')
ylabel('y [au]')
legend('initial','final')

subplot(3,1,3)
hold on
plot(time_days, z0)
plot(time_days, z)
xlabel('flight time [days]')
ylabel('z [au]')
legend('initial','final')

%% Orbit representation

% earth orbit
thetaE = linspace(0,2*pi,10*m);
xE = cos(thetaE);
yE = sin(thetaE);

% mars orbit
xM = 1.524*cos(thetaE);
yM = 1.524*sin(thetaE);

figure_orbits = figure(3);
hold on
set(figure_orbits,'position',[xpos + 2*1.2*imsize,ypos,1.2*imsize,imsize])
xlabel('x coordinate [AU]')
ylabel('y coordinate [AU]')
plot(0,0,'*k');
plot(x(1),y(1),'*k');
scatter(xE,yE,1,'k','filled');
scatter(xM,yM,1,'b','filled');
for i=1:m
    marker = plot(x(i),y(i),'r*');
    hold on
    
    if animations==1
        legend(num2str(i))
        frame = getframe(figure_orbits);
        im{i} = frame2im(frame);
    end
    
    delete(marker)
    plot(x(i),y(i),'r.','DisplayName','');
end

plot(x(m),y(m),'*k');


if animations==1
    frame = getframe(figure_orbits);
    im{i} = frame2im(frame);
    writegif('orbit.gif',im,m,2/m);
end

% figure(3)
% hold on
% xlabel('x coordinate [m]')
% ylabel('y coordinate [m]')
% zlabel('z coordinate [m]')
% plot3(0,0,0,'*k');
% scatter(P(1,j,i+1),P(2,j,i+1),40,c(i,:),'filled');
% for j=1:steps
%     plot3(r1(1,j),r1(2,j),r1(3,j),style);
%     F(j) = getframe;
% end
% 

end

