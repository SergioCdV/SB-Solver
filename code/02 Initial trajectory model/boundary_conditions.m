function Papp = boundary_conditions(initial, final, tfapp, n_init, r0)

% BOUNDARY CONDITIONS function for definition of initial and final control
% points using the initial and final spacecraft velocities.
% Coordinates are adimensionalised using r0.

% Note that index n=0 (for the first control point) is substituted to n=1
% (corresponding to the first row of the array), and all subsequent rows
% are shifted by one (due to matlab arrays starting at 1 instead of 0).


% initialization of array (4 initial points (n+1)
% and 3 coordinates (rho,theta,z, adim))
Papp = zeros(n_init+1,3);

% adimensionalization vector:
adimv = [r0, 1, r0];


% zeroth order (first row)
Papp(1,:) = initial.pos./adimv;

% first order (second row)
Papp(2,:) = initial.vel.*tfapp./n_init./adimv + initial.pos./adimv;

% the n-1 points (thrid row)
Papp(3,:) = final.pos./adimv - final.vel.*tfapp./n_init./adimv;

% the n points (fourth row)
Papp(4,:) = final.pos./adimv;


end