function B = bernstein0(n, j, tau)

% calculate bernstein basis polynomial

B = factorial(n).*tau.^j.*(1-tau).^(n-j)/...
    (factorial(j)*factorial(n-j));

end