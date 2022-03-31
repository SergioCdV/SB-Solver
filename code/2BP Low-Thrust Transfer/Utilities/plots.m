% function [] = plots3(C,C0,P,P0,B,tau,tfapp,amax,r0,mu,m,animations, n)

%% Set figures
% Setting up for gif image
imsize = 350;
im = cell(1,m);
xpos=10; ypos=150;

%% Final state vector
% Final time
time_days = tau.*(tf)*(365);

% Final trajectory in Cartesian coordinates
[S] = cylindrical2cartesian(C(1:3,:),true);
x = S(1,:);
y = S(2,:); 
z = S(2,:);

%% Orbit representation
figure_orbits = figure;
set(figure_orbits,'position',[xpos,ypos,1.2*imsize,imsize])
hold on
xlabel('$X$ coordinate [AU]')
ylabel('$Y$ coordinate [AU]')
plot(0,0,'*k');
plot(x(1),y(1),'*k');
viscircles([0,0],initial(1),'LineStyle','--','Color','r','LineWidth',0.3); % Earth's orbit
viscircles([0,0],final(1),'LineStyle','-.','Color','b','LineWidth',0.3);   % Mars' orbit
hold on
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
set(figure_propulsion,'position',[xpos + 1.2*imsize,ypos,1.2*imsize,imsize])
title('Spacecraft acceleration in time')
hold on
plot(time_days, sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2), 'k','LineWidth',1)
plot(time_days, u, 'LineWidth',0.3)
yline(4, '--k')
xlabel('Flight time [days]')
ylabel('$\vec{a}$ [m/$s^2$]')
legend('$a$','$a_\rho$','$a_\theta$','$a_z$')
grid on;

figure 
hold on
plot(time, unwrap(atan2(u(2,:),u(1,:)))); 
hold off 
grid on;
xlabel('Time')
ylabel('$\theta$')
title('Thrust in-plane angle')

figure 
hold on
plot(time, unwrap(atan2(u(3,:),sqrt(u(1,:).^2+u(2,:).^2)))); 
hold off 
grid on;
xlabel('Time')
ylabel('$\phi$')
title('Thrust out-of-plane angle')

%% Position coordinates
figure_coordinates = figure;
set(figure_coordinates,'position',[xpos + 2*1.2*imsize,ypos,1.2*imsize,imsize])
title('Spacecraft position coordinates in time')
hold on

% Spacecraft orbit 
x = C(1,:).*cos(C(2,:));
y = C(1,:).*sin(C(2,:));

% Earth's orbit
thetaE = 2*pi/365.25*time_days+initial(2);
xE = initial(1).*cos(thetaE);
yE = initial(1).*sin(thetaE);

% Mars's orbit
thetaM = 2*pi/687*time_days+final(2)-pi/2-0.075;
xM = final(1).*cos(thetaM);
yM = final(1).*sin(thetaM);

subplot(2,1,1)
hold on
plot(time_days, xM, 'LineStyle','-.','Color','b','LineWidth',0.3)
plot(time_days, xE, 'LineStyle','--','Color','r','LineWidth',0.3)
plot(time_days, x, 'k','LineWidth',1)
plot(time_days(1), x(1),'*k','DisplayName','')
plot(time_days(end),x(end),'*k','DisplayName','')
xlabel('Flight time [days]')
ylabel('$x$ [AU]')
grid on;

subplot(2,1,2)
hold on
plot(time_days, yM, '-.','LineWidth',0.3)
plot(time_days, yE, '--','LineWidth',0.3)
plot(time_days, y, 'k','LineWidth',1)
plot(time_days(1), y(1),'*k','DisplayName','')
plot(time_days(end),y(end),'*k','DisplayName','')
xlabel('Flight time [days]')
ylabel('$y$ [AU]')
grid on;