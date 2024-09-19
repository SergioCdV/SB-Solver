
function [s] = project_Hamiltonian(mu, s, E)
    for i = 1:5
        sigma = sqrt( mu ./ (-E .* dot(s(1:4,:), s(1:4,:), 1) + 2 * dot(s(5:8,:), s(5:8,:), 1)) );
        s = sigma .* s;
    end
end