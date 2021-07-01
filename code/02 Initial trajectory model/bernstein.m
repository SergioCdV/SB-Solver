function B = bernstein(n, j, tau, derivative)

% calculate bernstein basis polynomial and derivatives

switch derivative
    case 0
        B = bernstein0(n, j, tau);
    case 1
        B = bernstein1(n, j, tau);
    case 2
        B = bernstein2(n, j, tau);
    otherwise
        error('A higher-order Bernstein polynomial derivative is required, but has not been computed.')
end

end