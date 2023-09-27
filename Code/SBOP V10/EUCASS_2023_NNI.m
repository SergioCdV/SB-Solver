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
study_case = 'LQR';             % Select the optimal problem to be solved 
training_flag = true;
lambda_flag = true;

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

max_u = zeros(3,size(State,1));
for i = 1:size(State,1)
    u = control{i};
    max_u(:,i) = max(u,[],2);
end

index = max(max_u) < 1e1;
State = State(index);
control = control(index); 
time = time(index);

max_u = zeros(3,size(State,1));
for i = 1:size(State,1)
    u = control{i};
    max_u(:,i) = max(u,[],2);
end

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

%% Solve the CARE 
A(1:3,:) = [zeros(3) eye(3)];
A(4:6,:) = [0.4170 0.3023 0.1863 0.5388 0.2045 0.6705; ...
           0.7203 0.1468 0.3456 0.4192 0.8781 0.4173; ...
           0.0001 0.0923 0.3968 0.6852 0.0274 0.5587];
[S, ~, ~] = icare(A, [zeros(3); eye(3)], eye(6), eye(3), [], [], []);

train_lambda = cell(size(trainset));
test_lambda = cell(size(testset));

for i = 1:size(trainset,1)
    train_lambda{i} = [S * trainset{i}(1:6,:); S * trainset{i}(4:9,:)];
end

for i = 1:size(testset,1)
    test_lambda{i} = [S * testset{i}(1:6,:); S * testset{i}(4:9,:)];
end

%% Symplecticity test
% Omega = zeros(1, 101 * size(trainset,1));
% J = [zeros(6) eye(6); -eye(6) zeros(6)];
% 
% for i = 1:size(trainset,1)
%     for j = 1:size(trainset{i},2)
%         z = [trainset{i}(1:6,j); train_lambda{i}(1:6,j)];
%         dz = [trainset{i}(4:9,j); train_lambda{i}(4:9,j)];
%         Omega(1 + 101*(i-1) + j) = z.' * J * z;
%     end
% end

%% Training 
if (training_flag)
    fun = @(s,u,t)LQR_cost(s,u,t);
    num = size(trainset,1) * 101;

    % Traininig set 
    X = zeros(13, num);
    O = zeros(3, num);
    Oaux = zeros(12, num);

    for i = 1:size(trainset,1)
        s = trainset{i}; 
        u = train_control{i}; 
        t = train_time{i}; 
            
        max_u(:,i) = max(u,[],2);

        % Create the input vector
        new = [t; repmat(s(1:6,end), 1, length(t)); s(1:6,:)];
        X(:, 1 + num / size(trainset,1) * (i-1): num / size(trainset,1) * i) = new;
        O(:, 1 + num / size(trainset,1) * (i-1): num / size(trainset,1) * i) = u;

        if (lambda_flag)
            Oaux(:, 1 + num / size(trainset,1) * (i-1): num / size(trainset,1) * i) = [s(1:6,:); train_lambda{i}(1:6,:)];
        end
    end
    
    % Create and train the first network
    hiddenLayers = 2;                                             % Number of hidden layers
    neuronsLayer = 10;                                            % Number of neurons per layer
    neuronsLayer = repmat(neuronsLayer, 1, hiddenLayers);         % Layers definitions

    if (lambda_flag)
        lambda_net = fitnet(neuronsLayer);                        % Initial nets
        lambda_net.performFcn = 'myLoss';

        tic
        lambda_net = train(lambda_net, X, Oaux);
        training_time(1) = toc;

        Y = lambda_net(X);
        X = Y;
    end

    % Create and train the second network
    control_net = fitnet(neuronsLayer);
    
    tic
    control_net = train(control_net, X, O);
    training_time(2) = toc;
    
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
    X = [t; repmat(s(1:6,end), 1, length(t)); s(1:6,:)];
    if (lambda_flag)
        X = lambda_net(X);
    end

    % Assess the performance
    y = control_net(X);
    perf(i) = mse(y(1:3,:), u);
end
testing_time = toc;

metric = [mean(perf); mean(perf) + 3 * std(perf)]; 

% Generate the error bar graph 
E = linspace(min(perf), 2, 100); 
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
index = 15;
s = testset{index}; 
u = test_control{index}; 
t = test_time{index}; 
X = [t; repmat(s(1:6,end), 1, length(t)); s(1:6,:)];
if (lambda_flag)
    X = lambda_net(X);
end

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

if (lambda_flag)
    figure
    hold on 
    plot(t, test_lambda{index});
    scatter(t, X(10:end,:), 'filled')
    hold off
    xlabel('$t$')
    ylabel('$\mathbf{\lambda}$')
    grid on;
    xticklabels(strrep(xticklabels, '-', '$-$'));
    yticklabels(strrep(yticklabels, '-', '$-$'));
end


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

% Symplecticity 
function [cost] = Symplectic(x, y, dz, gamma)
    % Initial loss function 
    res = y-x;
    cost = sqrt(dot(res,res,1));

    % Symplecticity

    % Final cost function 
    cost = cost + gamma * res;
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