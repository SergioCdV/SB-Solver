
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




%% plot propulsive acceleration

 ACE4 = acceleration(P, mu, B, r0, tfapp,n);

figure_propulsion = figure(fig);
set(figure_propulsion,'position',[xpos + 1.2*imsize,ypos,1.2*imsize,imsize])
% title('propulsive acceleration')
hold on
% [arho, atheta, az] = aplot(C,r0,tfapp,mu);
plot(time_days, ACE4,'LineWidth',1)
% plot(time_days,arho, 'LineWidth',0.3)
% plot(time_days,atheta, 'LineWidth',0.3)
% plot(time_days,az, 'LineWidth',0.3)
% yline(amax, '--k')
xlabel('flight time [days]')
% ylabel('prop (final) [m/s^2]')
% legend('a total','a \rho','a \theta','a z')

ylabel('propulsive acceleration [m/s^2]')

% legend('jul-16','jul-23','jul-31','aug-07')