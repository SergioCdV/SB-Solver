function amag = acceleration(P, mu, B, r0, tfapp, n)

% calculates the acceleration magnitude for each instant

% initialize cell array
C = cell(1,3);

% extract coordinates
for k=1:3 % for coordinates
    for j=1:3 % for derivatives
        C{k}(:,j) = squeeze(B{k}(j,:,:))*P(k,1:(n(k)+1))';
    end
end

% [rho,   theta,   z]   = extract_coordinates(C(1,:,:),r0,1);
% [rho.D,  theta.D,  z.D]  = extract_coordinates(C(2,:,:),r0/tfapp,1/tfapp);
% [rho.DD, theta.DD, z.DD] = extract_coordinates(C(3,:,:),r0/tfapp^2,1/tfapp^2);


[rho, theta, z] = extract_coordinates(C, r0, tfapp);

% find distance to sun
r = sqrt(rho.o.^2 + z.o.^2);

% equations of motion
arho = rho.DD - rho.o.*theta.D.^2 + mu.*rho.o./r.^3;
atheta = rho.o.*theta.DD + 2.*rho.D.*theta.D;
az = z.DD + mu.*z.o./r.^3;

% magnitude of acceleration
amag = sqrt(arho.^2 + atheta.^2 + az.^2);

end