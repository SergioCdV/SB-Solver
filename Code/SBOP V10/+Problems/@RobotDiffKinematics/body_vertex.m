function [v] = body_vertex(index)
    % Lengths 
    a = [0 -0.24365 -0.21325 0 0 0];
    d = [0.15185 0 0 0.1124 0.08535 0.0921];
    l = a + d;

    % Radius 
    r = [0.08 0.08 0.065 0.06 0.06 0.06];

    % Main computation
    l = l(index); 
    r = r(index);

    v = [-r/2, -r/2, +l/2; ...
         +r/2, -r/2, +l/2; ...
         -r/2, +r/2, +l/2; ...
         +r/2, +r/2, +l/2; ...
         -r/2, -r/2, -l/2; ...
         +r/2, -r/2, -l/2; ...
         -r/2, +r/2, -l/2; ...
         +r/2, +r/2, -l/2];

    % Safety factor 
    v = 1.05 * v.';
end