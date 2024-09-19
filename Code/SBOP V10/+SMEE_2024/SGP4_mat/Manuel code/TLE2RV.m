function [t,r,v] = TLE2RV(TLE_Block,GM)


longstr1 = TLE_Block(1,:);
longstr2 = TLE_Block(2,:);

%     // set the implied decimal points since doing a formated read
%     // fixes for bad input data values (missing, ...)for (j = 11:16)

year    = str2double(longstr1(19:20));
day     = str2double(longstr1(21:32));
t       = (year*365 + day)*24*3600;

inclo = str2double(longstr2(9:16)); inclo = inclo*pi/180;        % Inclination                               [deg]
node = str2double(longstr2(18:25)); node = node* pi/180;    % Right Ascension of the Ascending Node     [deg]
ecc = str2double(['.' longstr2(27:33)]);                          % Eccentricity (with decimal point)         []
argp = str2double(longstr2(35:42)); argp = argp *pi/180;    % Argument of the perigee                   [deg]
m = str2double(longstr2(44:51));    m = m*pi/180;           % Mean Anomaly                              [deg]
n = str2double(longstr2(53:63));    n = n*2*pi/(24*3600);   % Mean Motion                               [revs/day]
revnum = str2double(longstr2(64:68));                       % Revolution number at epoch                [Revs]

a = (GM/n^2)^(1/3);     % m
p = a*(1-ecc^2);        % Semi-Latus Rectum

E       = fzero(@(E) -m + E - ecc*sin(E),m/(1-ecc)); %Eccentric Anomaly
nu      = 2*atan(sqrt((1+ecc)/(1-ecc))*tan(E/2)); % True Anomaly
arglat  = nu+argp;
x0      = sqrt(2 - 2*(cos(node)*cos(argp+nu) + sin(node)*sin(argp+nu)*sin(inclo)));
truelon = fzero(@(x) -cos(x) + cos(node)*cos(argp+nu) + sin(node)*sin(argp+nu)*sin(inclo),x0);
lonper  = node + argp;

[r,v] = coe2rv (p,ecc,inclo,node,argp,nu,arglat,truelon,lonper, GM);

% [p,ecc,inclo,node,argp,nu,arglat,truelon,lonper] = rv2coe(r,v,GM);
end
