

classdef Cost 
    % Fundamental definition of the cost function
    properties 
        path;           % Cost function file path 

        cost;           % Value of the cost function
        units;          % Units of the cost function

        CostFunction;   % Cost function
    end

    % Public methods
    methods 
        % Add cost function from txt file
        function [obj] = Cost(myCostFile)
            % Read the txt file
            eqs = fileread(strcat('Problem\',myCostFile));
            eqs = strsplit(eqs, newline);
            eqs = strtrim(eqs);

            % Remove empty lines 
            Eqs = {};
            k = 1;
            for i = 1:size(eqs,2)
                if (~isempty(eqs{i}))
                    Eqs{k} = eqs{i};
                    k = k+1;
                end
            end

            myFunc = cellfun(@(beta, t0, tf, x, u)(str2func(['@(beta,t0,tf,s,u)' beta, t0, tf, x, u])), Eqs);

            % Save the cost function 
            obj.CostFunction = myFunc; 
        end
        
        % Evaluate the cost function
        function [L, M] = evaluate_cost(obj, beta, t0, tf, s, u)
            % Evaluate the user defined cost function
            cost = feval(myFunc, beta, t0, tf, s, u);

            % Decompose the cost vector
            M = cost(1);    % Mayer penalty term
            L = cost(2);    % Lagrange integral term
        end
    end
end