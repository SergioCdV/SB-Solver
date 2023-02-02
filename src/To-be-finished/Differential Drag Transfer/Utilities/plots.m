% function [] = plots3(C,C0,P,P0,B,tau,tfapp,amax,r0,mu,m,animations, n)

%% Set figures
% Setting up for gif image
imsize = 350;
im = cell(1,m);
xpos=10; ypos=150;

%% Final state vector
% Final time
time_days = tau.*tf*t0/86400;

% Final spacecraft trajectory in Cartesian coordinates
r = C(1,:).*(1-C(2,:).^2)./(1+C(2,:).*cos(C(3,:)));
x = r.*cos(C(3,:));
y = r.*sin(C(3,:));

% Initial and final orbits
theta = linspace(0,2*pi,length(time_days)); 

r = coe_initial(1).*(1-coe_initial(2)^2)./(1+coe_initial(2).*cos(theta));
x0 = r.*cos(theta);
y0 = r.*sin(theta);

r = coe_final(1).*(1-coe_final(2)^2)./(1+coe_final(2).*cos(theta));
xf = r.*cos(theta);
yf = r.*sin(theta);

%% Orbit representation
figure_orbits = figure;
set(figure_orbits,'position',[xpos,ypos,1.2*imsize,imsize])
hold on
xlabel('$X$ coordinate')
ylabel('$Y$ coordinate')
plot(0,0,'*k');
plot(x(1),y(1),'*k');
plot(x0,y0,'LineStyle','--','Color','r','LineWidth',0.3);   % Earth's orbit
plot(xf,yf,'LineStyle','-.','Color','b','LineWidth',0.3);   % Mars' orbit
grid on; 

% Animation
if (animations == 1)
    for i = 1:m
        marker = plot(x(i),y(i),'k*');
        hold on
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
title('Transfer orbit')
plot(x,y,'k','LineWidth',1);
plot(x(m),y(m),'*k');
grid on;

%% Propulsive acceleration plot
% Plot
figure_propulsion = figure;
subplot(3,1,1)
set(figure_propulsion,'position',[xpos + 1.2*imsize,ypos,1.2*imsize,imsize])
hold on;
plot(time_days, u(1,:), 'LineWidth', 0.3)
yline(T, '--k')
xlabel('Flight time [days]')
ylabel('$B$')
grid on;n;

% subplot(3,1,2)
% hold on
% plot(time_days, u(2,:), 'LineWidth', 0.3)
% yline(dT, '--k')
% xlabel('Flight time [days]')
% ylabel('$\dot{B}$')
% grid on;
% 
% subplot(3,1,3)
% hold on
% plot(time_days, u(3,:), 'LineWidth', 0.3)
% yline(ddT, '--k')
% xlabel('Flight time [days]')
% ylabel('$\ddot{B}$')
% grid on;
% sgtitle('Spacecraft $B^*$ state in time')

%% Position coordinates
figure_coordinates = figure;
set(figure_coordinates,'position',[xpos + 2*1.2*imsize,ypos,1.2*imsize,imsize])
title('Spacecraft position coordinates in time')
hold on

subplot(3,1,1)
hold on
yline(coe_initial(1), '--b')
yline(coe_final(1), '--r')
plot(time_days, C(1,:), 'k','LineWidth',1)
plot(time_days(1), C(1,1),'*k','DisplayName','')
plot(time_days(end),C(1,end),'*k','DisplayName','')
xlabel('Flight time [days]')
ylabel('$SMA$ [R]')
grid on;

subplot(3,1,2)
hold on
yline(coe_initial(2), '--b')
yline(coe_final(2), '--r')
plot(time_days, C(2,:), 'k','LineWidth',1)
plot(time_days(1), C(2,1),'*k','DisplayName','')
plot(time_days(end),C(2,end),'*k','DisplayName','')
xlabel('Flight time [days]')
ylabel('$ECC$ []')
grid on;

subplot(3,1,3)
hold on
yline(coe_initial(end), '--b')
yline(coe_final(end), '--r')
plot(time_days, rad2deg(C(4,:)), 'k','LineWidth',1)
plot(time_days(1), rad2deg(C(4,1)),'*k','DisplayName','')
plot(time_days(end),rad2deg(C(4,end)),'*k','DisplayName','')
xlabel('Flight time [days]')
ylabel('$$\nu$ [deg]')
grid on;