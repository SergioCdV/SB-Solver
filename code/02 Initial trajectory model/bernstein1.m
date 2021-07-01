function B = bernstein1(n, j, tau)

% calculate first derivative of bernstein basis polinomial

% if j == 0
%     B = -n.*(1-tau).^(n-1);
% elseif j == n
%     B = n.*tau.^(n-1);
% else
%     B = factorial(n).*tau.^(j-1).*(1-tau).^(n-j)/...
%     (factorial(j-1)*factorial(n-j)) - ...
%     factorial(n).*tau.^j.*(1-tau).^(n-j-1)/...
%     (factorial(j)*factorial(n-j-1));
% end
% 
% end


if j == 0
    B = -n.*(1-tau).^(n-1);
elseif j < n
    B = factorial(n).*tau.^(j-1).*(1-tau).^(n-j)/...
    (factorial(j-1)*factorial(n-j)) - ...
    factorial(n).*tau.^j.*(1-tau).^(n-j-1)/...
    (factorial(j)*factorial(n-j-1));
elseif j == n
    B = n.*tau.^(n-1);
else
    error('wrong n and j input into bernstein polinomial (1st derivative)')
end

end