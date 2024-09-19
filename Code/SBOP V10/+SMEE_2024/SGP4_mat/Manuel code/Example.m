
clear all;
clc;
close all;
format long;

addpath('HPOP')
addpath('SGP4')

%path2mice = '..\08Guillermo_Tesis\Codenew\npi-manoeuvre-detection\Mice_Repository\MicePackage\';
path2mice = '/Users/msanrivo/Matlab/mice';
addpath(genpath(path2mice));
global ME

ME = load_spice_kernels ( path2mice );

HPOPconfig; % load configuration files

%% Example ECI to TLE transformation

etstr = '2020-08-29T10:12:11.633'; %2020242101211.663
et = cspice_str2et(etstr);
xECI = [5562.4879416571 3523.1684782116 2148.4429892465 -1.0301056926 5.0440695064 -5.5767521291]';

% Direct transformation
tle = etstate2tle(et,xECI,PARAMS);

% With rotation matrix as input
[T,dpsi,deps] = ECI2TEMErotparams(PARAMS,et);
tm = roteci2temeSGP(T,dpsi,deps);

tle2 = etstate2tle(et,xECI,tm);

disp(tle-tle2);

% transform back to ECI
xECIp = ettle2state(et,tle,PARAMS);

xECIp2 = ettle2state(et,tle,tm);

disp([xECI-xECIp xECI-xECIp2])


%% Example propagation using HPOP
% S/C parameters
PARAMS.AuxParam.mass    = 1;
PARAMS.AuxParam.area_drag = 1;
PARAMS.AuxParam.Cd      = 1;
PARAMS.AuxParam.area_solar  = 1;
PARAMS.AuxParam.Cr = 1;
 
% Define the environment for the simulation

PARAMS.AuxParam.n       = 80;
PARAMS.AuxParam.m       = 80;
PARAMS.AuxParam.sun     = 1;
PARAMS.AuxParam.moon    = 1;
PARAMS.AuxParam.planets = 1;
PARAMS.AuxParam.sRad    = 1;
PARAMS.AuxParam.drag    = 1;
PARAMS.AuxParam.SolidEarthTides = 1;
PARAMS.AuxParam.OceanTides = 1;
PARAMS.AuxParam.Relativity = 1;
PARAMS.AuxParam.Anelastic = 1;
PARAMS.AuxParam.JB = 1; % use JB08 instead of nrlmsise00
PARAMS.options = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',1000);
PARAMS.integrators{1} = 'ode113';

PARAMS.nx = 8; % number of state variables
S0 = xECI;
S0(7) = 1.3*0.2/260;
S0(8) = 2.0*0.2/260;

model_handle = @(et,state)HPOPModel(et,state,PARAMS);

deltat = 100;
etf = et + 3600;
dt = et:deltat:etf;

[times,states] = feval(PARAMS.integrators{1},model_handle,dt,S0',PARAMS.options);


PARAMS.AuxParam.JB = 0; % use nrlmsise00
model_handle = @(et,state)HPOPModel(et,state,PARAMS); % need to modify model_handle...
[times2,states2] = feval(PARAMS.integrators{1},model_handle,dt,S0',PARAMS.options);

disp(states2(end,:)'-states(end,:)');

