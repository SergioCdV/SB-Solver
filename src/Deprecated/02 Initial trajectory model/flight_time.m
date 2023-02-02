function tf = flight_time(P, B, m, tfapp, r0, n)

% FLIGHT_TIME function for to calculate total flight time from the set of
% points.


% The coordinates are calculated using the points and the Bernstein
% matrices

% extract coordinates
for k=1:3 % for coordinates
    for j=1:3 % for derivatives
        C{k}(:,j) = squeeze(B{k}(j,:,:))*P(k,1:(n(k)+1))';
    end
end


[rho, theta, z] = extract_coordinates(C, r0, tfapp);
v = sqrt(rho.D.^2 + (rho.o.*theta.D).^2 + z.D.^2);

dt = zeros(1,m-1);

for i=2:m % for intervals (discretization points)
   dr = sqrt((rho.o(i)-rho.o(i-1))^2 + (z.o(i)-z.o(i-1))^2 + (rho.o(i)*(theta.o(i)-theta.o(i-1)))^2);
   dt(i) = dr / abs(v(i));
end

% multiplied by the initial estimate to recover dimension
tf = sum(dt);


end