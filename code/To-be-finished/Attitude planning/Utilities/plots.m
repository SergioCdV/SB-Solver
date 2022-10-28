%% Project: 
% Date: 19/05/22

%% Plots %%
% Function to display the plots of the optimization

% Inputs: - system, the structure of the problem containing relevant data
%         - scalar tf, the final time of flight 
%         - vector tau, the final temporal grid 
%         - array C, the final state evolution matrix
%         - array u, a 3xm matrix with the control input evolution
%         - scalar T, the maximum allowed acceleration
%         - structure setup, containing the setup of the figures

function plots(system, tf, tau, C, u, T, setup)
    % Set up 
    animations = setup.animations;
    
    % Constants
    time = tf*tau; 

    % Setting up for gif image
    imsize = 350;
    im = cell(1,size(C,2));
    xpos = 10; ypos = 150;

    % Final state vector
    time_days = tau.*tf/3600;
    
    % Final spacecraft trajectory in Cartesian coordinates
    n = system.V; 
    x = zeros(3,size(C,2)); 
    for i = 1:size(C,2)
        aux = quaternion_product([0; n], quaternion_inverse(C(1:4,i)));
        aux = quaternion_product(C(1:4,i), aux);
        x(:,i) = aux(2:4);
    end
    
    % Orbit representation
    figure_orbits = figure;
    view(3)
    set(figure_orbits,'position',[xpos,ypos,1.2*imsize,imsize])
    xlabel('$X$ coordinate')
    ylabel('$Y$ coordinate')
    zlabel('$Z$ coordinate')
    hold on
    quiver3(zeros(3,1),zeros(3,1),zeros(3,1),[1;0;0],[0;1;0],[0;0;1])
    [X, Y, Z] = sphere(128);
    h = surfl(X, Y, Z); 
    set(h, 'FaceAlpha', 0.1)
    shading interp
    plot3(x(1,1),x(2,1),x(3,1),'*k');         % Initial conditions
    plot3(x(1,:),x(2,:),x(3,:),'b');          % Initial conditions
    plot3(x(1,end),x(2,end),x(3,end),'*k');   % Final conditions
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

    % Torque acceleration plot
    figure_propulsion = figure;
    set(figure_propulsion,'position',[xpos + 1.2*imsize,ypos,1.2*imsize,imsize])
    hold on
    plot(time_days, sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2), 'k','LineWidth',1)
    plot(time_days, u, 'LineWidth', 0.3)
    yline(T, '--k')
    xlabel('Time [hours]')
    ylabel('$\mathbf{M}$')
    legend('$M$','$M_x$','$M_y$','$M_z$')
    grid on;

    figure 
    hold on
    plot(time, rad2deg(atan2(u(2,:),u(1,:)))); 
    hold off 
    grid on;
    xlabel('Time [hours]')
    ylabel('$\theta$')
    title('Torque in-plane angle')

    figure 
    hold on
    plot(time, rad2deg(atan2(u(3,:),sqrt(u(1,:).^2+u(2,:).^2)))); 
    hold off 
    grid on;
    xlabel('Time [hours]')
    ylabel('$\phi$')
    title('Torque out-of-plane angle')
end