
%% prelinimaries


% for gif image
imsize = 350;
im = cell(1,m);
xpos=10; ypos=150;

% time calculations
tF = flight_time(P, B, m, tfapp, r0, n);
time = tau.*tF;
time_days = time./(60*60*24);


%% convert cylindrical coordinates to cartesian
[x, y, z] = cylindrical_to_cartesian(C);


%% Orbit representation

fig = fig + 1;


figure_propulsion = figure(fig);
set(figure_propulsion,'position',[xpos,ypos,1.2*imsize,imsize])

title(strcat(num2str(iter1),num2str(iter2),num2str(iter3)))
hold on
[arho, atheta, az] = aplot(C,r0,tfapp,mu);
plot(time_days, acceleration(P, mu, B, r0, tfapp,n), 'k','LineWidth',2)
plot(time_days,arho, 'LineWidth',0.4)
plot(time_days,atheta, 'LineWidth',0.4)
plot(time_days,az, 'LineWidth',0.4)
yline(amax, '--k')
xlabel('flight time [days]')
ylabel('prop (final) [m/s^2]')
legend('a total','a \rho','a \theta','a z')



function [arho, atheta, az] = aplot(C,r0,tfapp,mu)



[rho, theta, z] = extract_coordinates(C, r0, tfapp);


% find distance to sun
r = sqrt(rho.o.^2 + z.o.^2);

% equations of motion
arho = rho.DD - rho.o.*theta.D.^2 + mu.*rho.o./r.^3;
atheta = rho.o.*theta.DD + 2.*rho.D.*theta.D;
az = z.DD + mu.*z.o./r.^3;

end