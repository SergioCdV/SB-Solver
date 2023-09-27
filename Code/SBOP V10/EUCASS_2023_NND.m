%% Embedding Pontryagin's Principle in Neural Networks for Optimal Asteroid Landing %% 
% Date: 08/08/2022

%% EUCASS 2023 Neural Networks
% This script provides the particular function implemented to analyze the
% feasibility of replacing high-order direct transcription methods with
% neural networks

% Inputs:  
% Outputs: 

clear;
close all; 
set_graphics();

rng("default")                          % For reproducibility of the data partition

%% Input dataset and flags
study_case = 'Optimal landing';  % Select the optimal problem to be solved 
GNRR_flag = false;               % Use GNRR architectures
data_augmentation = true;       % Use data augmentation via SBOP
training_flag = true;

switch (study_case)
    case 'Optimal landing'
        load("landing_set"); 
    case 'LQR'
        load("lqr_set"); 
    otherwise 
        error('No valid control problem was selected');
end

% Dimensioning 
State = State.';
control = control.'; 
time = time.';

%% Generation of the training and data set 
% Matrix shuffling
shuffler = randperm(size(State,1));
State = State(shuffler);
control = control(shuffler);
time = time(shuffler); 

% Partition
c = cvpartition(length(State), "Holdout", 0.2);                % Test partition
training_index = training(c);                                  % Training set indices
trainset = State(training_index,:);                            % Training state dataset
train_control = control(training_index,:);                     % Training control dataset
train_time = time(training_index);                             % Training sampling grid dataset

test_index = test(c);                                          % Test set indices
testset = State(test_index);                                   % Test state dataset
test_control = control(test_index);                            % Test control dataset
test_time = time(test_index);                                  % Test sampling grid dataset

%% Data augmentation 
if (data_augmentation)
    % Regression 
    basis = 'Legendre';                   % Polynomial basis to be use. Alternatively: Legendre, Bernestein, Orthogonal Bernstein
    time_distribution = 'Legendre';       % Distribution of time intervals. Alternatively: Bernstein, Orthogonal Bernstein, Chebsyhev, Legendre, Linear, Newton-Cotes, Normal, Random, Trapezoidal
    n = 7;                                % Polynomial order in the state vector expansion

    for i = 1:size(trainset,1)

        C = trainset{i};
        m = size(C,2)-1;
        switch (study_case)
            case 'LQR'
                n = [7 7 7];
                L = 2;
            case 'Optimal landing'
                n = [7 7 7];
                L = 2;
        end

        if (L > 1)
            Caux = [C(1:3,:) 0.5 * (train_time{i}(end) -  train_time{i}(1)) * C(1+3:2*3,:)];
        else
            Caux = C(1:3,:);       
        end
    
        solver = Solver(basis, n, time_distribution, m);
        Grid = solver.gridding(m);
        B = solver.state_basis(L, n, basis, Grid.tau);
    
        P = zeros(length(n), max(n)+1); 
        for j = 1:length(n)
            if (L > 1)
                A = [B{j}(1:n(j)+1,:) B{j}(n(j)+2:2*(n(j)+1),:)];
            else
                A = B{j}(1:n(j)+1,:);
            end
            P(j,1:n(j)+1) = Caux(j,:)*pinv(A);
        end
    
        % Compute the new solution
        m = 1e3;
        Grid = solver.gridding(m);
        B = solver.state_basis(L, n, basis, Grid.tau);
        C = solver.evaluate_state(n, L, P, B);
        C(4:6,:) = C(4:6,:)/ (0.5 * (train_time{i}(end) -  train_time{i}(1)));
        C(7:9,:) = C(7:9,:)/ (0.5 * (train_time{i}(end) -  train_time{i}(1)))^2;
        switch (study_case)
            case 'LQR'
                A = [0.4170 0.3023 0.1863 0.5388 0.2045 0.6705 ...
                     0.7203 0.1468 0.3456 0.4192 0.8781 0.4173 ...
                     0.0001 0.0923 0.3968 0.6852 0.0274 0.5587];

                problem_params = [train_time{i}(end); reshape(A, [], 1)];
                OptProblem = Problems.LQR(C(1:6,1), C(1:6,end), L, 3, 3, problem_params);
            case 'Optimal landing'
                problem_params = [0.0033; 1; 2; 1; 0.0408];
                OptProblem = Problems.AsteroidLanding(C(1:6,1), C(1:6,end), L, 3, 3, problem_params);
        end
        u = OptProblem.ControlFunction(problem_params, [],  train_time{i}(1), train_time{i}(end), Grid.tau, C);
        t = (Grid.tau + 1) * 0.5 * ( train_time{i}(end) - train_time{i}(1));
    
        % Save the data 
        trainset{i} = C;
        train_control{i} = u; 
        train_time{i} = t; 
    end
else
    m = 1;
end

%% Training 
if (training_flag)
    switch (study_case)
        case 'LQR'
            fun = @(s,u,t)LQR_cost(s,u,t);
            num = size(trainset,1) * max(101,m+1);
        case 'Optimal landing'
            fun = @(s,u,t)landing_cost(s,u,t);
            num = size(trainset,1) * max(21,m+1);
    end

    % Traininig set 
    X = zeros(23, num);
    O = zeros(3, num);

    for i = 1:size(trainset,1)
        s = trainset{i}; 
        u = train_control{i}; 
        t = train_time{i}; 
    
        % Create the input vector
        new = [t; repmat([s(1:6,1); s(1:6,end)], 1, length(t)); s; feval(fun, s, u, t)];
        X(:, 1 + num / size(trainset,1) * (i-1): num / size(trainset,1) * i) = new;
        O(:, 1 + num / size(trainset,1) * (i-1): num / size(trainset,1) * i) = u;
    end
    
    % Create and train the network
    tic
    if (GNRR_flag)
        control_net = newgrnn(X, O, 0.4);         % Generalized regression NN
    else
        % Define the number of neurons in the net
        hiddenLayers = 4;                                              % Number of hidden layers
        neuronsLayer = 20;                                             % Number of neurons per layer
        neuronsLayer = repmat(neuronsLayer, 1, hiddenLayers);          % Layers definitions
        control_net = fitnet(neuronsLayer);                            % Initial nets

        control_net = train(control_net, X, O);
    end
    training_time = toc;
    
    % Save the net 
    % save study_case;
end

%% Testing 
% Test the net
perf = zeros(1,size(testset,1));
tic
for i = 1:size(testset,1)
    s = testset{i}; 
    u = test_control{i}; 
    t = test_time{i}; 

    % Create the input vector 
    X = [t; repmat([s(1:6,1); s(1:6,end)], 1, length(t)); s; feval(fun, s, u, t)];

    % Assess the performance
    y = control_net(X);
    perf(i) = mse(y(1:3,:), u);
end
testing_time = toc;

metric = [mean(perf); mean(perf) + 3 * std(perf)]; 

% Generate the error bar graph 
E = linspace(min(perf), max(perf), 100); 
f = zeros(length(E),1); 

for i = 1:length(perf)
    for j = 1:length(E)-1
        if (perf(i) < E(j+1))
            f(j) = f(j)+1;
            break; 
        end
    end
end

%% Results 
figure
scatter(1:size(testset,1), perf)
hold off
xlabel('Mission time [days]')
ylabel('$e$')
grid on;

% Create the input vector 
index = length(testset);
s = testset{index}; 
u = test_control{index}; 
t = test_time{index}; 
X = [t; repmat([s(1:6,1); s(1:6,end)], 1, length(t)); s; feval(fun, s, u, t)];
ur = control_net(X);
error = u-ur;

figure 
hold on 
plot(t, u);
scatter(t, ur, 'filled')
hold off
xlabel('$t$')
ylabel('$\mathbf{u}$')
grid on;
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));

figure 
plot(t, error)
xlabel('$t$')
ylabel('$\mathbf{\epsilon}$')
grid on;
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));

figure
bar(E, f/size(testset,1)*100);
xlabel('RMSE')
ylabel('Number of cases $[\%]$')
grid on; 

%% Auxiliary functions
% LQR cost function 
function [cost] = LQR_cost(s,u,t)
    % Trapezoidal integration 
    cost = cumtrapz(t, dot(s,s,1) + dot(u,u,1));
end

% LQR cost function 
function [cost] = landing_cost(s,u,t)
    % Trapezoidal integration 
    cost = cumtrapz(t, ones(1,size(u,2)));
end


% Set up cool graphics 
function set_graphics()
    %Set graphical properties
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
    set(groot, 'defaultAxesFontSize', 11); 
    set(groot, 'defaultAxesGridAlpha', 0.3); 
    set(groot, 'defaultAxesLineWidth', 0.75);
    set(groot, 'defaultAxesXMinorTick', 'on');
    set(groot, 'defaultAxesYMinorTick', 'on');
    set(groot, 'defaultFigureRenderer', 'painters');
    set(groot, 'defaultLegendBox', 'off');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultLegendLocation', 'best');
    set(groot, 'defaultLineLineWidth', 1); 
    set(groot, 'defaultLineMarkerSize', 3);
    set(groot, 'defaultTextInterpreter','latex');
end