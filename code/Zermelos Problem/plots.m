figure(1)
plot(C(1,:),C(2,:),'b'); 
grid on; 
xlabel('$x$ coordinate [m]')
ylabel('$y$ coordinate [m]'); 
title('Boat trajectory in time')

figure(3)
subplot(2,1,1)
hold on
plot(time, C(1,:),'b'); 
plot(time, C(2,:),'r'); 
hold off
grid on; 
legend('$x$', '$y$')
xlabel('Time [s]')
ylabel('$\mathbf{r}$ [m]')
subplot(2,1,2)
hold on
plot(time, C(3,:),'b'); 
plot(time, C(4,:),'r'); 
hold off
grid on; 
legend('$v_x$', '$v_y$')
xlabel('Time [s]')
ylabel('$\mathbf{v}$ [m/s]')
sgtitle('Phase space evolution in time')

figure(2)
plot(time,theta,'b'); 
grid on; 
xlabel('Time [s]')
ylabel('$\theta$ [rad]')
title('Steering angle evolution in time')