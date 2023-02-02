function c = acceleration_constraints(r, r0, aO, ac)

% function to evaluate acceleration for each given instant
% to see if the constraints (defined by the solar sail) are met

axO = aO(1);
ayO = aO(2);
azO = aO(3);

radim = r0/r;

param = azO/ac/radim;

if (param >= 0) && (param <= 2/3)
    c = sqrt(axO^2+ayO^2) - azO/2/sqrt(2);
elseif (param <= 1 ) && (param > 2/3)
    c = sqrt(axO^2+ayO^2) - ac*radim*sqrt(1/16-(param-3/4)^2);
else
    error('There is an error in the definition of acceleration constraints.')
end


end