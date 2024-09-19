function TLE = readTLE(tline)
%
% This function loads the TLE elements into the structure TLE. The input is
% an array of char, 
%   tline   char [2,68]
%
% The information of a TLE is as follows:
% 19-32	04236.56031392	Element Set Epoch (UTC)
% 3-7	25544	Satellite Catalog Number
% 9-16	51.6335	Orbit Inclination (degrees)
% 18-25	344.7760	Right Ascension of Ascending Node (degrees)
% 27-33	0007976	Eccentricity (decimal point assumed)
% 35-42	126.2523	Argument of Perigee (degrees)
% 44-51	325.9359	Mean Anomaly (degrees)
% 53-63	15.70406856	Mean Motion (revolutions/day)
% 64-68	32890	Revolution Number at Epoch
%
% Parameters
muE = 398600;   % km3/s2
%
% Read first line
%
TLE.Cnum = tline(1, 3:7);      			        % Catalog Number (NORAD)
TLE.SC   = tline(1, 8);					        % Security Classification
TLE.ID   = tline(1, 10:17);			            % Identification Number
TLE.year = str2double(tline(1, 19:20));            % Year
TLE.doy  = str2double(tline(1, 21:32));            % Day of year
TLE.epoch = str2double(tline(1, 19:32));              % Epoch
TLE.TD1   = str2double(tline(1, 34:43));              % first time derivative
TLE.TD2   = str2double(tline(1, 45:50));              % 2nd Time Derivative
TLE.ExTD2 = tline(1, 51:52);                       % Exponent of 2nd Time Derivative
TLE.BStar = str2double(tline(1, 54:59));              % Bstar/drag Term
TLE.ExBStar = str2double(tline(1, 60:61));            % Exponent of Bstar/drag Term
TLE.BStar = TLE.BStar*1e-5*10^TLE.ExBStar;
TLE.Etype = tline(1, 63);                          % Ephemeris Type
TLE.Enum  = str2double(tline(1, 65:end));          % Element Number
%
% Read second line
%
TLE.inc     = str2double(tline(2 ,9:16));                % Orbit Inclination (degrees)
TLE.raan    = str2double(tline(2, 18:25));               % Right Ascension of Ascending Node (degrees)
TLE.ecc     = str2double(strcat('0.',tline(2, 27:33)));  % Eccentricity
TLE.omega   = str2double(tline(2, 35:42));               % Argument of Perigee (degrees)
TLE.M       = str2double(tline(2, 44:51));               % Mean Anomaly (degrees)
TLE.no      = str2double(tline(2, 53:63));               % Mean Motion
TLE.sma     = ( muE/(TLE.no*2*pi/86400)^2 )^(1/3);    % semi major axis (m)
TLE.rNo     = str2double(tline(2, 65:end));              % Revolution Number at Epoch


end