function DV = velocity_variation(P, mu, B, r0, tau, tfapp, n)

% Computation of the performance index to be minimized in optimization, the
% total DeltaV needed in the trajectory

% Compute magnitude of necessary propulsive acceleration in flight
amag = acceleration(P, mu, B, r0, tfapp, n);

% % summation for integral
% dt = tau(2) - tau(1);
% DV = dt*sum(amag);

% in case the time intervals are not constant
DV = 0;
for i = 2:length(tau)
    dt = (tau(i) - tau(i-1))*tfapp;
    DV = DV + abs(amag(i))*dt;
end





