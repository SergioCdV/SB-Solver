

function [B] = dipole_model(mu, COE, v)
    % Constants of the model
    Mref = 7.83e15;                 % Reference dipole moment
    gamma = deg2rad( 11.44 );       % Tilt angle
    omega = 2*pi / (86400);         % Earth's angular motion

    % Orbital constants
    i = COE(4);             % Orbital inclination
    RAAN = COE(3);          % Right ascension of the ascending node

    % Orbital motion 
    h = sqrt( mu * COE(1) * (1-COE(2)^2));                  % Angular momentum
    r = COE(1) * (1-COE(2)^2) ./ (1 + COE(2)*cos(v));       % Position vector

    % Time law 
    t = cumtrapz(v, h./r.^2);

    % Frame transformation
    beta = omega * t - RAAN;

    eta = atan2(-sin(gamma), sin(i)*cos(gamma)-cos(i)*sin(gamma)*cos(beta));
    xi = atan2( -sin(gamma)./sin(eta), cos(i)*cos(gamma)+sin(i)*sin(gamma)*cos(beta) );
    theta = v - eta;

    % Magnetic model
    B =  Mref ./ r.^3 .* [sin( xi ) .* cos(theta); -cos( xi );  2 * sin( xi ) .* sin(theta)];
end