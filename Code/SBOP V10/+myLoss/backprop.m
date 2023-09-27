function dy = backprop(t,y,e,param)
%MSE.BACKPROP Backpropagate derivatives of performance

% Copyright 2012-2015 The MathWorks, Inc.

  lambda = 10; 
  J = [zeros(6) eye(6); -eye(6) zeros(6)]; 
  de = 2 * (t([1:6 10:end],1) * J * y([1:6 10:end],1)) * J * t([1:6 10:end],1);
  dy = -2 .* e - lambda * de;
end
