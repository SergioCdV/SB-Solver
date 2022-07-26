% clearvars

%% Version 4 of program ( 07/04/2021 )
% Following the paper "Initial design..." by Fan et. al



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% %%%%%%%%%%%%%%%%%%%%%%   PRE-OPTIMISATION   %%%%%%%%%%%%%%%%%%%%%%%%%%% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

%% Variables to be defined for each run

m = 10; % number of discretization points
time_distribution = 'linear'; % distribution of time intervals
sigma = 1; % if normal distribution is selected
amax = 1.5e-4; % maximum acceleration available [m/s^2]


animations = 0; % set to 1 for gif
fig = 1; % figure start number

% order of Bezier curve functions for each coordinate
% in the format n = [n_rho, n_theta, n_z]
n = [12,12,8];

n = [6,6,4];



%% Global constants
r0 = 149597870700; % 1 au [m] (for adimansionalisation)
mu = 1.32712440042e+20; % gravitational parameter sun [m^3 s^−2]


%% Initial definitions

% creation of time array for iterations
if time_distribution == 'linear'
    tau = linspace(0,1,m); % non-dimensional time
elseif time_distribution == 'normal'
    pd = makedist('Normal');
    pd.sigma = sigma;
    xpd = linspace(-3,3,m);
    tau = cdf(pd,xpd);
else
    error('An appropriate time array distribution must be specified')
end

% initial coordinates and orbital elements are defined
[initial, final] = initial_coordinates2(r0);



%% Initial guess (for 4 control points)
% tfapp is the initial estimate of flight time
% Papp is the initial estimate of 4 control points (from boundary conditions)
% Bapp is the initial calculation of bersntein polynomials from Papp
% Capp is the coordinate estimations from these. 

[tfapp, Papp, Bapp, Capp] = initial_approximation(initial, final, tau, m, mu, amax, r0);


%% Fitting for n+1 control points

% Once the initial curve has been calculated, a series of n+1 points are
% determiend which will fit this curve, where n is the approximation order
% of the Bézier optimization. These points will be optimised later.

 % the bernstein polynomials calculated here will be useful for the
 % optimisation process, and will not need to be calculated repeatedly.
[B, P0, C0] = initial_fitting(n, tau, Capp, m);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% %%%%%%%%%%%%%%%%%%%%%%%%   OPTIMISATION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% define upper and lower bounds (empty in this case)
P_lb = [];
P_ub = [];

%
% P_lb = P0(:,1).*ones(size(P0));
% P_lb(:,n(1)+1) = P0(:,n(1)+1);
% P_ub = P0(:,n(1)+1).*ones(size(P0));
% P_ub(:,1) = P0(:,1);

% Optimisation: define objective function (minimise total DeltaV
objective = @(P)velocity_variation(P, mu, B, r0, tau, tfapp, n);

% define non-linear constraints (on max acceleration and
% initial and final points)
nonlcon = @(P)constraints(P, P0, B, amax, mu, r0, tfapp, n);

% Modification of fmincon optimisation options and parameters
% (according to the details in the paper)
options = optimoptions('fmincon','TolCon',1e-9,'TolFun',1e-6,'Display','iter-detailed');
options.MaxFunctionEvaluations = 5e+03;

% optimise that mfckr
[P,tf_final,exitflag,output]  =   fmincon(objective,P0,[],[],[],[],P_lb,P_ub,nonlcon,options);
% (the empty cells are for the linear constraints, not used in this case)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% %%%%%%%%%%%%%%%%%%%%%%%%   RESUSLTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% (such as there are)

% Coordinate calculation
C = cell(1,3);

% extract coordinates
for k=1:3 % for coordinates
    for j=1:3 % for derivatives
        C{k}(:,j) = squeeze(B{k}(j,:,:))*P(k,1:(n(k)+1))';
    end
end

%%

display_results(P0, P, B, m, exitflag, output, tfapp, r0, n)
plots3



