%
% Test propagation with SGP4 in Starlink data
% Starlink constellation, inc 70 deg., 09 Oct - 08 Nov 23
%
clear vars
close all 
%
% Parameters
%
muE     = 398600;   % km3/s2
TWOPI   = 2*pi;
MINUTES_PER_DAY = 1440;
MINUTES_PER_DAY_SQUARED = (MINUTES_PER_DAY * MINUTES_PER_DAY);
MINUTES_PER_DAY_CUBED = (MINUTES_PER_DAY * MINUTES_PER_DAY_SQUARED);
% 
% Load Data
T   = readtable('st_20231108_1726475408.csv');
%
% Choice: NORAD ID = 49142
%
auxInd             = ( T(:,3).Variables == 49142 ); 
Taux               = T(auxInd,:); 
%
%
% Load TLEs at times t0, t1, ..., tii, ..., tN
for ii = 1:height(Taux)
    tleBlock(1, :) = char(Taux(ii,24).Variables); 
    tleBlock(2, :) = char(Taux(ii,25).Variables); 
    %
    if ii == 2
        %
        % State at time t1
        %
        [t1,r1,v1]  = TLE2RV(tleBlock,muE); 
    end
    %
    TLE(ii)         = readTLE(tleBlock);
    %
end
%
% SGP4 propagation from initial epoch to t1
%
satdata.epoch           = TLE(1).epoch;
satdata.norad_number    = TLE(1).Cnum;
satdata.bulletin_number = TLE(1).ID;
satdata.classification  = TLE(1).SC; % almost always 'U'
satdata.revolution_number = TLE(1).rNo;
satdata.ephemeris_type  = TLE(1).Etype;
satdata.xmo             = TLE(1).M * (pi/180);
satdata.xnodeo          = TLE(1).raan * (pi/180);
satdata.omegao          = TLE(1).omega * (pi/180);
satdata.xincl           = TLE(1).inc * (pi/180);
satdata.eo              = TLE(1).ecc;
satdata.xno             = TLE(1).no * TWOPI / MINUTES_PER_DAY;
satdata.xndt2o          = TLE(1).TD1 * 1e-8 * TWOPI / MINUTES_PER_DAY_SQUARED;
satdata.xndd6o          = TLE(1).TD2 * TWOPI / MINUTES_PER_DAY_CUBED;
satdata.bstar           = TLE(1).BStar;

tsince = (TLE(2).epoch - TLE(1).epoch) * MINUTES_PER_DAY; % amount of time in which you are going to propagate satellite's state vector forward (+) or backward (-) [minutes] 
%
[rteme, vteme] = sgp4(tsince, satdata);
% 
% 
%