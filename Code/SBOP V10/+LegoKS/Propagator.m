

function [T, Y] = Propagator(func_handle, myProblem, int_options, prop_config)

    % Sanity checks 
    if ~exist('int_options', 'var')
        int_options = [];                        % No additional integration options
    end

    if ~exist('prop_config', 'var')
        prop_config = LegoKS.PropConfig();
    end
    
    % Propagation
    if isa(myProblem, 'PropProblems.NBodyProblemMH')
        [T, Y] = NBodyPropagator(func_handle, myProblem, int_options, prop_config);
    elseif isa(myProblem, 'PropProblems.TwoBodyProblem')
        [T, Y] = TwoBodyPropagator(func_handle, myProblem, int_options, prop_config);
    else
        error('The problem class is not supported. Aborting...')
    end
end

%% Auxiliary functions
% For N bodies 
function [T, Y] = NBodyPropagator(func_handle, myProblem, int_options, prop_config)
    % New initial conditions
    m = 8;                                  % Standard state space dimension (KS regularization)

    if mod(length( myProblem.IC ), m)
        warning('Initial conditions are not regularized... Aborting integration')
    else
        % Constants
        n = length( myProblem.IC ) / m;     % Number of states
        
        if ( prop_config.PhaseCnt )
            prop_config.PhaseInt = true;
        end

        % Integrate the fiber
        s0 = myProblem.IC;              % Initial conditions for the particle
        M = n * m;
    
        if (prop_config.PhaseInt)
            s0 = [s0; zeros(n,1)];                     % Add the phase to the initial conditions 
        end

        % Propagation
        if ( ~prop_config.OrthoProj && ~prop_config.EnrgyProj )
            [T, y] = func_handle( @(t,s)myProblem.Dynamics(prop_config, myProblem.params, t, s), myProblem.t_span, s0, int_options );
            Y = y.';

        else
            % Pre-allocation
            x = zeros(M, length(myProblem.t_span));
            x(:,1) = s0;

            for j = 1:length(myProblem.t_span)-1
                % Integration
                [~, saux] = func_handle( @(t,s)myProblem.Dynamics(prop_config, myProblem.params, t, s), [myProblem.t_span(j) myProblem.t_span(j+1)], s0, int_options );
                saux = saux(end,:);
                F = reshape(saux.', [], 1);

                % Projection onto the bilinear constraint 
                if ( prop_config.OrthoProj )
                    for i = 1:n
                        F(5+m*(i-1):m*i,:) = LegoKS.project_bilinear( F(1+m*(i-1):4+m*(i-1),:), F(5+m*(i-1):m*i,:) );
                    end
                end

                x(:,j+1) = F;
                s0 = saux.';
            end
            
            % Final assembly of the trajectory
            T = myProblem.t_span;
            Y = x;
        end
    end
end

% For two bodies
function [T, Y] = TwoBodyPropagator(func_handle, myProblem, int_options, prop_config)
    % New initial conditions
    m = 8;                                  % Standard state space dimension (KS regularization)

    if mod(length( myProblem.IC ), m)
        warning('Initial conditions are not regularized... Aborting integration')
    else
        % Constants
        n = length( myProblem.IC ) / m;     % Number of states
        
        % Pre-allocation
        Y = cell(1,n);                      % Set of trajectories
        T = cell(1,n);                      % Set of clocks


        if ( prop_config.EnrgyProj )
            prop_config.EnergyInt = true;
        end

        if ( prop_config.PhaseCnt )
            prop_config.PhaseInt = true;
        end

        % Integrate the fiber
        for i = 1:n
            s0 = myProblem.IC(1+m*(i-1):m*i);              % Initial conditions for the particle
            M = m;
    
            if (prop_config.EnergyInt)                                                         % Add energy to the initial conditions
                E0 = LegoKS.OscEnergy(myProblem.mu, s0, prop_config.SundmanTransformation );   % Energy initial conditions
                s0 = [s0; E0];                                                                 % Complete initial conditions

                M = M + 1;
            end
        
            if (prop_config.PhaseInt)
                % Add the phase to the initial conditions   
                s0 = [s0; 0];

                M = M + 1;
            end

            % Propagation
            if ( ~prop_config.OrthoProj && ~prop_config.EnrgyProj )
                [t, y] = func_handle( @(t,s)myProblem.Dynamics(prop_config, myProblem.params, t, s), myProblem.t_span, s0, int_options );
                Y{i} = y.';
                T{i} = t;
            else
                % Pre-allocation
                x = zeros(M, length(myProblem.t_span));
                x(:,1) = s0;

                for j = 1:length(myProblem.t_span)-1
                    % Integration
                    [~, saux] = func_handle( @(t,s)myProblem.Dynamics(prop_config, myProblem.params, t, s), [myProblem.t_span(i) myProblem.t_span(i+1)], s0, int_options );
                    saux = saux(end,:);
                    F = reshape(saux.', [], 1);

                    % Projection onto the bilinear constraint 
                    if ( prop_config.OrthoProj )
                        F(5:8,:) = LegoKS.project_bilinear( F(1:4,:), F(5:8,:) );
                    end

                    % Projection onto the energy manifold 
                    if ( prop_config.EnrgyProj )
                        F(1:8,:) = LegoKS.project_Hamiltonian( myProblem.mu, F(1:8,:), F(9,:) );
                    end

                    x(:,j+1) = F;
                    s0 = saux.';
                end
                
                % Final assembly of the trajectory
                T{i} = myProblem.t_span;
                Y{i} = x;
            end
        end

        if (n == 1)
            Y = Y{1}; 
            T = T{1};
        end
    end
end