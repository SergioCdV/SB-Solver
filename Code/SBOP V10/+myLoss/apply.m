function perfs = apply(t,y,e,param)
%MSE.APPLY Calculate performances

% Copyright 2012-2015 The MathWorks, Inc.
  lambda = 10; 
  J = [zeros(6) eye(6); -eye(6) zeros(6)]; 
  omega = (t([1:6 10:end],1) * J * y([1:6 10:end],1))^2;
  perfs = e .* e + lambda * omega;
end
