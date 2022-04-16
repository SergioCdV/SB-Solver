function [B, P0, C0] = initial_fitting(n, tau, Capp, m)

% initialziation of arrays
B = cell(1,3);
C0 = cell(1,3);

for k=1:3
    B{k} = zeros(3,m,n(k)+1);
    C0{k} = zeros(m,3);
end

P0 = zeros(3,max(n)+1); % 3 coordinates, order n


%% The Bernstein-basis polinomials for the increasd order are calculated

for k = 1:3 % for coordinates
    for i = 1:3 % for derivatives
        for j=0:n(k)
           B{k}(i,:,j+1) = bernstein(n(k),j,tau,(i-1));
        end
    end
end

%% The points are fitted

for k = 1:3 % for coordinates
    P0(k,1:(n(k)+1)) = (squeeze(B{k}(1,:,:))\squeeze(Capp(1,:,k))')';
end


for k=1:3 % for coordinates
    for j=1:3 % for derivatives
        C0{k}(:,j) = squeeze(B{k}(j,:,:))*P0(k,1:(n(k)+1))';
    end
end

end