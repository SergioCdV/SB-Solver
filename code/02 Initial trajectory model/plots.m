function [] = plots(C,Capp,m,tau,P,amag)

%% Plots

Fig_coordinates = figure(1);
coordlabels = ["\rho [AU]" , "\theta [rad]" "z [AU]"];
for i=1:3
    subplot(1,3,i)
    hold on
    plot(tau,squeeze(Capp(1,:,i)),'k')
    plot(tau,squeeze(C(1,:,i)),'b')
    xlabel('\tau')
    title(coordlabels(i))
    legend('approx', 'final')
end
hold off

%%

Fig_rho_derivatives = figure(2);
derivlabels = ["0th derivative of \rho" , "1st derivative of \rho", "2nd derivative of \rho"];
for i=1:3
    subplot(1,3,i)
    hold on
    plot(tau,squeeze(Capp(i,:,1)))
    plot(tau,squeeze(C(i,:,1)))
    xlabel('\tau')
    title(derivlabels(i))
    legend('approx', 'final')
end

plot([0,1],[0,0],':k')
legend('approx', 'final', '')

%%

Fig_teta_derivatives = figure(3);
derivlabels = ["0th derivative of \theta" , "1st derivative of \theta", "2nd derivative of \theta"];
for i=1:3
    subplot(1,3,i)
    hold on
    plot(tau,squeeze(Capp(i,:,2)))
    plot(tau,squeeze(C(i,:,2)))
    xlabel('\tau')
    title(derivlabels(i))
    legend('approx', 'final')
end

plot(tau,zeros(1,m),':k')
legend('approx', 'final', '')

%%

Fig_velocities = figure(4);
vellabels = ["\rho speed" , "\theta speed", "z speed"];
for i=1:3
    subplot(1,3,i)
    hold on
    plot(tau,squeeze(Capp(2,:,i)),'k')
    plot(tau,squeeze(C(2,:,i)),'b')
    xlabel('\tau')
    title(vellabels(i))
    legend('approx', 'final')
end
hold off

%% 

Fig_acceleration = figure(5);
plot(tau,amag)
xlabel('\tau')
title("Propulsive acceleration magnitude [m/s^2]")


end