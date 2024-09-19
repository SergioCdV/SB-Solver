
function [H] = Hamiltonian(mu, u, du, E)
    H = dot(du, du, 1) - E / 2 .* dot(u, u, 1) - mu/2;
end