

function [dE] = EnergyDer(t, s, ap)
    % Pre-allocation 
    dE = zeros(1, size(s,2)); 

    % Compute the work done by the non-Keplerian forces in KS space
    for i = 1:size(dE,2)
        dE(1,i) = 2 * s(5:8,i).' * LegoKS.KSmatrix( s(1:4,i) ).' * ap(:,i);
    end
end