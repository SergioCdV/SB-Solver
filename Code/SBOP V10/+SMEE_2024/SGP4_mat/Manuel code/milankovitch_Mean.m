function dx = milankovitch_Mean(t, x)
%
% This function computes the evolution of the Mean Milankvitch elements
% according to the paper: 
% Nonsingular vectorial formulation of the short-periodic corrections in
% Kozai's oblateness solution, Paolo Izzo, Lamberto Dell'Elce, Pini Gurfil,
% Aaron J. Rosengren, Celestial Mechanics and Dynamical Astronomy, 2022
% 134:12
%
% State vector
%   x(1) = ex component
%   x(2) = ey component
%   x(3) = ez component
%   x(4) = hx component
%   x(5) = hy component
%   x(6) = hz component
%   x(7) = mean longitude
%
% Parameters
%
muE = 398600;   % km3/s2
J2  = 1.0826158e-3;     %  
RE  = 6371;     % km
%
% Earth's spin axis
%
pvec = zeros([3,1]); 
pvec(3) = 1; 
%
% Eccentricity vector and modulus
%
evec = x(1:3);  
emod = norm(evec); 
%
% Angular momentum vector and modulus
%
Hvec = x(4:6); 
Hmod = norm(Hvec); 
Hdir = Hvec / Hmod; 
%
% Mean longitude
%
lmean   = x(7); 
%
% Auxiliary scalars
%
pp      = Hmod^2 / muE; 
sma     = pp / (1 - emod^2); 
nmean   = sqrt( muE / sma^3 ) ;  
zeta    = dot( Hdir, pvec ); 
%
RHS_e_coeff = - 3 * nmean * J2 * RE^2 / 4 / pp^2; 
RHS_H_coeff = 3 * Hmod * nmean * J2 * RE^2 / 2 / pp^2; 
RHS_l_coeff = 3 * nmean * J2 * RE^2 / 4 / pp^2; 
%
% Auxiliary vectors
%
RHS_e1 = cross( Hdir, evec ); 
RHS_e2 = cross( pvec, evec ); 
RHS_H  = cross( Hdir, pvec ); 
%
RHS_e  = ( 1 - 5 * zeta) * RHS_e1 + 2 * zeta * RHS_e2; 
%
%
% Equations of Motion 
%
dx(1:3)   = RHS_e_coeff * RHS_e; 
%
dx(4:6)   = RHS_H_coeff * zeta * RHS_H; 
%
dx(7)     = nmean + RHS_l_coeff * ( sqrt(1 - emod^2) * (3 * zeta^2 - 1 ) + ...
    5 * zeta^2 - 2 * zeta - 1); 
%
dx        = dx'; 
%
end