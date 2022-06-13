%% Project: 
% Date: 19/05/22

%% Plots %%
% Function to display the plots of the optimization

% Inputs: - structure system, containing the physical information of the
%           2BP of interest
%         - scalar tf, the final time of flight 
%         - vector tau, the final temporal grid 
%         - array C, the final state evolution matrix
%         - array u, a 3xm matrix with the control input evolution
%         - scalar T, the maximum allowed acceleration
%         - vector initial, the initial landing Cartesian state vector
%         - vector final, the final landing Cartesian state vector
%         - structure setup, containing the setup of the figures

function plots(system, tf, tau, C, u, T, initial, final, setup)
    % Set up 
    animations = setup.animations;

    % Constants 
    t0 = system.time;       % Fundamental time unit of the system
    time = tf*tau; 

    a = system.ellipsoid/system.distance; 

    % Setting up for gif image
    imsize = 350;
    im = cell(1,size(C,2));
    xpos = 10; ypos = 150;

    % Final state vector
    time_days = tau.*tf*t0/86400;
    
    % Final spacecraft trajectory in Cartesian coordinates
    x = C(1,:);
    y = C(2,:); 
    z = C(3,:);
    
    % Orbit representation
    figure_orbits = figure;
    view(3)
    set(figure_orbits,'position',[xpos,ypos,1.2*imsize,imsize])
    hold on
    xlabel('$X$ coordinate')
    ylabel('$Y$ coordinate')
    zlabel('$Z$ coordinate')
    plot3(0,0,0,'*k');
    plot3(x(1),y(1),z(1),'*k');
    ellipsoid(0,0,0,a(1),a(2),a(3));   % Earth's orbit
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
    title('Landing trajectory')
    plot3(x,y,z,'k','LineWidth',1);
    plot3(x(end),y(end),z(end),'*k');
    grid on;

    % Propulsive acceleration plot
    figure_propulsion = figure;
    set(figure_propulsion,'position',[xpos + 1.2*imsize,ypos,1.2*imsize,imsize])
    title('Spacecraft acceleration in time')
    hold on
    plot(time_days, sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2), 'k','LineWidth',1)
    plot(time_days, u, 'LineWidth', 0.3)
    yline(T, '--k')
    xlabel('Flight time [days]')
    ylabel('$\mathbf{a}$')
    legend('$a$','$a_\rho$','$a_\theta$','$a_z$')
    grid on;

    figure 
    hold on
    plot(time, rad2deg(atan2(u(2,:),u(1,:)))); 
    hold off 
    grid on;
    xlabel('Time')
    ylabel('$\theta$')
    title('Thrust in-plane angle')

    figure 
    hold on
    plot(time, rad2deg(atan2(u(3,:),sqrt(u(1,:).^2+u(2,:).^2)))); 
    hold off 
    grid on;
    xlabel('Time')
    ylabel('$\phi$')
    title('Thrust out-of-plane angle')
end