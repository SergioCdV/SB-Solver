%% Example_1 %%
% Milankovitch vs SGP4
% Author: Sergio Cuevas 
% Date: 12/12/2023

% This script aims to verify the performance of the Milankovitch elements
% propagator under J2 perturbation against the standard SGP4 model. 
% 
% To do so, the osculating, propagated TLE data are transformed to the
% corresponding angular momentum and eccentricity vectors, realised in the
% TEME frame. These are then compared to those obtained by means of the
% Milankovitch propagator. 
% 
% Error performance is then analysed statistically in terms of the relative
% error between the two models obtained for the angular momentum magnitude,
% the eccentricity vector norm and the attitude error between the perifocal
% reference frames. Note that if the two models are exactly equivalent, all
% erros should be zero.

clear; 
clc;

%% Input data
file_name = 'SGP4-VER.TLE';     % Test data file name
output_name = 'tcppver.out';    % Output verification data

fid = fopen(file_name);         % Open file
tlecnt = 1;                     % TLE counter
tles = cell(tlecnt);            % Preallocation of a TLE cell array structure
line1 = fgetl(fid);

while ( ischar(line1) )
    if (line1(1) == '1')
        line2 = fgetl(fid);
        tle = TLE(line1, line2);    % Create the TLE object
        tles{tlecnt} = tle;         % Save the TLE object for its later use
        tlecnt = tlecnt + 1;
    end
    line1 = fgetl(fid);
end

fclose(fid);                    % Close file

%% Verification and validation of the propagator
% Preallocation of variables
curobj = 0;     % Current object index
strind = 0;
prevrv = 0;
ind = 0;
cnt = 0;

% Loop over the output verification file to compute Cartesian elements for various TLEs and epochs
fid = fopen(output_name);
line1 = fgetl(fid);

while ( ischar(line1) )
    strind = strfind(line1, 'xx'); 

    if (strind > 0)
        % Get the object number
        curobj = str2double(strtrim(line1(1:strind-1)));

        % Find tle with the matching object number
        ind = 1;
        while (ind <= tlecnt && tles{ind}.objectNum ~= curobj)
            ind = ind+1;
        end
        tle = tles{ind};

    elseif (curobj > 0)
        %% SGP4 MODEL PROPAGATION
        verout = sscanf(line1, '%f');   % TLE output data for verification purposes

        % Propagate to the given epoch by means of the SGP4 model
        elapsed_epoch = verout(1);      % Elapsed time to propagate since the generation of the TLE (initial conditions)
        rv = tle.getRV(elapsed_epoch);  % TEME propagated Cartesian elements
        mu = tle.rec.mu;                % Gravitational constant used (check for correct units in the other propagator!!!)
        
        % Original code reused RV so if there is an error we need to carry old values
        if(tle.sgp4Error > 0)
            rv = prevrv;
        end

        % Compute the osculating, propagated SPG4 perifocal triad (angular momentum, eccentricity and Hamilton's vector) in the TEME frame
        cnt = cnt + 1;
        h(:,cnt) = cross(rv(:,1), rv(:,2));                                     % Specific angular momentum
        e(:,cnt) = cross(rv(:,2), h(:,cnt)) / mu - rv(:,1) / norm(rv(:,1));     % Eccentricity vector
        y(:,cnt) = cross(h(:,cnt), e(:,cnt));                                   % Hamilton's vector (semiminor axis vector)

        %% Compute the very same quantities by propagating the Milankovitcch elements
        % Get the initial Milankovitch elements conditions in the TEME frame from the TLE (which is approximately equivalent to Brouwer's elements)

        h_vec = sqrt(mu * tle.rec.a * tle.rec.radiusearthkm * (1-tle.rec.ecco^2)) * [0; 0; 1];         % Angular momentum vector in the perifocal frame          
        ecc = tle.rec.ecco * [cos(tle.rec.argpo); sin(tle.rec.argpo); 0];                              % Eccentricity vector in the perifocal frame
        l = tle.rec.nodeo + tle.rec.argpo + tle.rec.mo;                                                % True longitude (RAAN + AoP + M)

        % Euler rotation matrix
        Q1 = [cos(tle.rec.nodeo) sin(tle.rec.nodeo) 0; -sin(tle.rec.nodeo) cos(tle.rec.nodeo) 0; 0 0 1];       % RAAN rotation
        Q2 = [cos(tle.rec.inclo) 0 sin(tle.rec.inclo); 0 1 0; -sin(tle.rec.inclo) 0 cos(tle.rec.inclo)];       % Inclination rotation
        Q3 = [cos(tle.rec.argpo) sin(tle.rec.argpo) 0; -sin(tle.rec.argpo) cos(tle.rec.argpo) 0; 0 0 1];       % AoP rotation
        R = Q3 * Q2 * Q1;
        R = R.';                % Rotation from the perifocal to an inertial frame

        % Initial conditions in the "TEME" frame 
        h_vec = R * h_vec;
        ecc = R * ecc;

        % Propagate them 
        elapsed_epoch = elapsed_epoch * 60;         % Propagation time in seconds

        % Save the results for the further processing
        hm(:,cnt) = cross(rv(:,1), rv(:,2));                                     % This shall be replaced by the Milankovitch propagation and osculating transformation
        em(:,cnt) = cross(rv(:,2), h(:,cnt)) / mu - rv(:,1) / norm(rv(:,1));     % This shall be replaced by the Milankovitch propagation and osculating transformation
        ym(:,cnt) = cross(h(:,cnt), e(:,cnt));                                   % This shall be replaced by the Milankovitch propagation and osculating transformation
        %%
        prevrv = rv;
    end

    line1 = fgetl(fid);
end
fclose(fid);

%% Analysis
% Orbit geometry error quantities
h_error = 100 * (1 - sqrt(dot(h, hm, 1)) / sqrt(dot(h, h, 1)));       % Propagated angular momentum relative error to the SGP4 reference
mu_herr = mean(h_error);                                              % Mean of the angular momentum relative error to the SGP4 reference
sigma_herr = std(h_error);                                            % Standard deviation of the angular momentum relative error to the SGP4 reference

e_error = 100 * (1 - sqrt(dot(e, em, 1)) / sqrt(dot(e, e, 1)));       % Propagated angular momentum relative error to the SGP4 reference
mu_eerr = mean(e_error);                                              % Mean of the angular momentum relative error to the SGP4 reference
sigma_eerr = std(e_error);                                            % Standard deviation of the angular momentum relative error to the SGP4 reference

% Perifocal error quantities 
q_error = zeros(1,length(h_error));

for i = 1:cnt 
    P = [e(:,i) / norm(e(:,i)) y(:,i) / norm(y(:,i)) h(:,i) / norm(h(:,i))]; 
    Pm = [em(:,i) / norm(em(:,i)) ym(:,i) / norm(ym(:,i)) hm(:,i) / norm(hm(:,i))];
    dP = P * Pm.';
    q_error(i) = real( acos( max(0.5 * (trace(dP)-1), 1) ) );
end

%% Display
figure
subplot(2,1,1)
hold on
plot(1:cnt, h_error + 3 * sigma_herr)
yline(mu_herr, 'k--')
scatter(1:cnt, h_error, 'filled')
ylabel('Angular momemtum error [%]')
legend('3sigma', 'Mean', 'Error')
grid on;

subplot(2,1,2)
hold on
plot(1:cnt, e_error + 3 * sigma_eerr)
yline(mu_eerr, 'k--')
scatter(1:cnt, e_error, 'filled')
ylabel('Eccentricity error [%]')
legend('3sigma', 'Mean', 'Error')
grid on;
xlabel('TLE case')

figure
hold on
plot(1:cnt, q_error + 3 * std(q_error))
yline(mean(q_error), 'k--')
scatter(1:cnt, q_error, 'filled')
ylabel('Perifocal attitude error [deg]')
legend('3sigma', 'Mean', 'Error')
grid on;
xlabel('TLE case')