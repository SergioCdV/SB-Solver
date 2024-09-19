%
% Script to check motion Starlink constellation, inc 70 deg., 09 Oct - 08 Nov 23
% 
%
clearvars
close all
% 
% Parameters
muE = 398600;   % km3/s2
% Load Data
T   = readtable('st_20231108_1726475408.csv');
%
% Process satellites one by one
%
arrayIDs = zeros([height(T),1]); 
%
% satS(:,1)     = NORAD CATALOG ID
% satS(:,2)     = MJD
% satS(:,3:5)   = hvec
% satS(:,6:8)   = evec
% satS(:,9)     = Mean longitude
% satS(:,10)    = Inclination
% satS(:,11)    = RAAN
% satS(:,12)    = Argument of perigee
% satS(:,13)    = Period
% satS(:,14)    = Mean motion
%
[IDordered,iT,iID] = unique(T(:,"NORAD_CAT_ID"));
nTLEperO           = accumarray(iID,1); 
maxRep             = max(nTLEperO); 
satS               = zeros([height(IDordered), 14, maxRep]); 
%
for ii = 1:height(IDordered)
    %
    % Loop in the individual satellites
    %
    %
    % Norad Cat ID
    %
    nCID = T(iT(ii),3).Variables; 
    %
    auxInd         = ( iID(iT(ii)) == iID ); 
    %
    for jj = 1:nTLEperO(ii)
        %
        % Compute the Milankovitch elements from TLEs
        %        
        Taux           = T(auxInd,:); 
        %
        tleBlock(1, :) = char(Taux(jj,24).Variables); 
        tleBlock(2, :) = char(Taux(jj,25).Variables); 
        %
        [MJD, rvec, vvec]      = TLE2RV(tleBlock,muE); 
        %
        % Angular Momentum
        %
        hvec   = cross(rvec,vvec); 
        %
        % Eccentricity Vector
        %
        evec   = cross(vvec, hvec)/muE - rvec/norm(rvec); 
        %
        % Mean longitude
        % 
        lmean   = Taux(jj,13).Variables + Taux(jj,14).Variables + Taux(jj,15).Variables; 
        %
        % Load data in array
        %
        satS(ii,1,jj)     = nCID; 
        satS(ii,2,jj)     = MJD;
        satS(ii,3:5, jj)  = hvec;
        satS(ii,6:8, jj)  = evec; 
        satS(ii,9,jj)     = lmean; 
        satS(ii,10, jj)   = Taux(jj,12).Variables; 
        satS(ii,11, jj)   = Taux(jj,13).Variables; 
        satS(ii,12, jj)   = Taux(jj,14).Variables; 
        satS(ii,13, jj)   = Taux(jj,29).Variables; 
        satS(ii,14, jj)   = Taux(jj,10).Variables; 
        %
    end
end
%
% Plots
%
%
for ii = 1:height(IDordered)
    for jj = 1:nTLEperO(ii)
        %
        hmod    = sqrt(satS(ii,3, jj).^2 + satS(ii,4, jj).^2 + satS(ii,5, jj).^2); 
        emod    = sqrt(satS(ii,6, jj).^2 + satS(ii,7, jj).^2 + satS(ii,8, jj).^2); 
        %
        figure(1)
        %
        hold on
        %
        plot3(satS(ii,3, jj)/hmod, satS(ii,4, jj)/hmod, satS(ii,5, jj)/hmod, 'ko'); 
        %
        figure(2)
        %
        hold on
        %
        plot3(satS(ii,6, jj)/emod, satS(ii,7, jj)/emod, satS(ii,8, jj)/emod, 'ro'); 
        %
        figure(3)
        %
        hold on
        %
        plot(emod, hmod, 'bo')
        %
    end
end