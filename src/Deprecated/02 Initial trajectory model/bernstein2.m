function B = bernstein2(n, j, tau)

% calculate second derivative of bernstein basis polinomial

if j == 0
    B = n.*(n-1).*(1-tau).^(n-2);
elseif j == 1
    B = n.*(n-1).*(n-2).*tau.*(1-tau).^(n-3)-2.*n.*(n-1).*(1-tau).^(n-2);
elseif j == (n-1)
    B = n.*(n-1).*(n-2).*tau.^(n-3).*(1-tau)-2.*n.*(n-1).*tau.^(n-2);
elseif j == n
    B = n.*(n-1).*tau.^(n-2);
else
    B = factorial(n).*tau.^(j-2).*(1-tau).^(n-j)/...
    (factorial(j-2)*factorial(n-j)) - ...
    2*factorial(n).*tau.^(j-1).*(1-tau).^(n-j-1)/...
    (factorial(j-1)*factorial(n-j-1)) + ...
    factorial(n).*tau.^j.*(1-tau).^(n-j-2)/...
    (factorial(j)*factorial(n-j-2));
end

end