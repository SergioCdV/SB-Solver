
classdef KeplerianProblem < PropProblems.TwoBodyProblem
    % This is the definition of the class for the Keplerian problem in KS
    % space

    properties
        myLc;                   % Characteristic longitude
        myFc;                   % Characteristic frequency
        SolType;                % String to select the type of propagated solution desired: either numerical or analytical
    end

    methods
        % Constructor
        function obj = KeplerianProblem(myMu, myIC, myTspan, mySolType)
            % Dimensional quantities
            [Lc, Fc, mu] = PropProblems.KeplerianProblem.DimensionalQts(myMu, myIC);     % Dimensional quantities
            myIC = myIC ./ [repmat(Lc, 3, 1); repmat(Lc * Fc, 3, 1)];                    % Non-dimensional Cartesian state vector
            myTspan = myTspan * Fc;                                                      % Non-dimensional time

            % Initial conditions
            ksIC = LegoKS.KS_mapping( myIC, true, "1", mu );                             % Initial KS fiber                                                         
    
            % Reshape initial conditions 
            ksIC = reshape(ksIC, [], 1);

            % Parent constructor
            super_arguments{1} = ksIC;
            super_arguments{2} = myTspan;
            super_arguments{3} = [];

            obj@PropProblems.TwoBodyProblem(super_arguments{:});

            % Dimensionalization
            obj.myLc = Lc; 
            obj.myFc = Fc;
            obj.params = [mu];
            obj.SolType = mySolType;
        end

        % Methods
        function [ap] = ForceModel(obj, t, s)
            ap = zeros(3, size(s,2));
        end

        function [t,y] = AnalyticalSol(obj, params, t, s, u)
        end
    end

    methods (Static)
        function [Lc, Fc, mu] = DimensionalQts(myMu, myIC)
            % Compute an initial longitude 
            Lc = sqrt( dot(myIC, myIC, 1) ); 
            Fc = sqrt( myMu ./ Lc.^3 );
            mu = 1;
        end
    end
end