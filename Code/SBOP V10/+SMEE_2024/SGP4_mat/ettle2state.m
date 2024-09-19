
function x = ettle2state(et,tle,PARAMS)

tsince = 1e-16;

ae = 1.0;
tothrd = (2.0./3.0);
XJ3 = -2.53881e-6;
e6a = 1.0E-6;
xkmper = 6378.135;
ge = 398600.8; % Earth gravitational constant
CK2 = (1.0826158e-3 ./ 2.0);
CK4 = (-3.0 .* -1.65597e-6 ./ 8.0);

% Variables
xno = tle(1);
eo = tle(2);
xincl = tle(3);
xnodeo = tle(4);
omegao = tle(5);
xmo = tle(6);
try
bstar = tle(7);
catch
bstar = 0;
end

% Constants
s = ae + 78 ./ xkmper;
qo = ae + 120 ./ xkmper;
xke = sqrt((3600.0 .* ge) ./ (xkmper.^3));
qoms2t = ((qo - s).^2).^2;
temp2 = xke ./ (xno);
a1 = temp2.^tothrd;
cosio = cos (xincl);
theta2 = (cosio.^2);
x3thm1 = 3.0 .* theta2 - 1.0;
eosq = (eo.^2);
betao2 = 1.0 - eosq;
betao = sqrt(betao2);
del1 = 1.5 .* CK2 .* x3thm1 ./ ((a1.^2) .* betao .* betao2);
ao = a1 .* ( 1.0 - del1.*((1.0./3.0) + del1 .* (1.0 + (134.0./81.0) .* del1)));
delo = 1.5 .* CK2 .* x3thm1 ./ ((ao.^2) .* betao .* betao2);
xnodp = (xno)./(1.0 + delo);
aodp = ao./(1.0 - delo);
% Initialization
% For perigee less than 220 kilometers, the isimp flag is set and
% the equations are truncated to linear variation in sqrt a and
% quadratic variation in mean anomaly.  Also, the c3 term, the
% delta omega term, and the delta m term are dropped.
% isimp = 0;
% if ((aodp * (1.0 - satdata.eo)/ ae) < (220.0/xkmper + ae))
%     isimp = 1;
% end
isimp = 0.*xno;
isimp((aodp .* (1.0 - eo)./ ae) < (220.0./xkmper + ae)) = 1;
% For perigee below 156 km, the values of s and qoms2t are altered.
s4 = 0.*xno + s;
qoms24 = qoms2t;
perige = (aodp .* (1.0 - eo) - ae) .* xkmper;
% if (perige < 156)
%     s4 = perige - 78.0;
%     if (perige <= 98)
%         s4 = 20.0;
%     end
%     qoms24 = (((120.0 - s4) * ae / xkmper)^4.0);
%     s4 = s4 / xkmper + ae;
% end
s4(perige < 156) = perige(perige < 156) - 78.0;
s4(perige <= 98) = 20.0;
qoms24 = qoms24 + zeros(size(perige));
ae = ae + zeros(size(perige));
xkmper = xkmper + zeros(size(perige));
qoms24(perige < 156) = (((120.0 - s4(perige < 156)) .* ae(perige < 156) ./ xkmper(perige < 156)).^4.0);
s4(perige < 156) = s4(perige < 156) ./ xkmper(perige < 156) + ae(perige < 156);

pinvsq = 1.0 ./ ( (aodp.^2) .* (betao2.^2) );
tsi = 1.0 ./ (aodp - s4);
eta = aodp .* (eo) .* tsi;
etasq = (eta.^2);
eeta = (eo) .* eta;
psisq = abs( 1.0 - etasq);
coef = qoms24 .* (tsi.^4.0);
coef1 = coef ./ (psisq.^3.5);
c2 = coef1 .* xnodp .* (aodp .* (1.0 + 1.5 .* etasq + eeta .* (4.0 + etasq)) + 0.75 .* CK2 .* tsi ./ psisq .* x3thm1 .* (8.0 + 3.0 .* etasq .* (8.0 + etasq)));
c1 = (bstar) .* c2;
sinio = sin(xincl);
a3ovk2 = -XJ3 ./ CK2 .* (ae.^3.0);
c3 = coef .* tsi .* a3ovk2 .* xnodp .* ae .* sinio ./ (eo);
x1mth2 = 1.0 - theta2;
c4 = 2.0 .* xnodp .* coef1 .* aodp .* betao2 .* ( eta .* (2.0 + 0.5 .* etasq) + (eo) .* (0.5 + 2.0 .* etasq) - 2.0 .* CK2 .* tsi ./ (aodp .* psisq) .* ( -3.0 .* x3thm1 .* ( 1.0 - 2.0 .* eeta + etasq .* (1.5 - 0.5.*eeta)) + 0.75 .* x1mth2 .* (2.0 .* etasq - eeta .* (1.0 + etasq)) .* cos(2.0 .* (omegao))));
c5 = 2.0 .* coef1 .* aodp .* betao2 .* (1.0 + 2.75 .* (etasq + eeta) + eeta .* etasq);
theta4 = (theta2.^2);
temp1 = 3.0 .* CK2 .* pinvsq .* xnodp;
temp2 = temp1 .* CK2 .* pinvsq;
temp3 = 1.25 .* CK4 .* pinvsq .* pinvsq .* xnodp;
xmdot = xnodp + 0.5 .* temp1 .* betao .* x3thm1 + 0.0625 .* temp2 .* betao .* (13.0 - 78.0 .* theta2 + 137.0 .* theta4);
x1m5th = 1.0 - 5.0 .* theta2;
omgdot = -0.5 .* temp1 .* x1m5th + 0.0625 .* temp2 .* (7.0 - 114.0 .* theta2 + 395.0 .* theta4) + temp3 .* (3.0 - 36.0 .* theta2 + 49.0 .* theta4);
xhdot1 = -temp1 .* cosio;
xnodot = xhdot1 + (0.5 .* temp2 .* (4.0 - 19.0 .* theta2) + 2.0 .* temp3 .* (3.0 - 7.0 .* theta2)) .* cosio;
omgcof = (bstar) .* c3 .* cos(omegao);
xmcof = -(2.0./3.0) .* coef .* (bstar) .* ae ./ eeta;
xnodcf = 3.5 .* betao2 .* xhdot1 .* c1;
t2cof = 1.5 .* c1;
xlcof = 0.125 .* a3ovk2 .* sinio .* (3.0 + 5.0 .* cosio) ./ (1.0 + cosio);
aycof = 0.25 .* a3ovk2 .* sinio;
delmo = ((1.0 + eta .* cos(xmo)).^3);
sinmo = sin(xmo);
x7thm1 = 7.0 .* theta2 - 1.0;
% if (isimp==0)	
%     c1sq = (c1^2);
%     d2 = 4.0 * aodp * tsi * c1sq;
%     temp = d2 * tsi * c1 / 3.0;
%     d3 = (17.0 * aodp + s4)*temp;
%     d4 = 0.5 * temp * aodp * tsi * (221.0 * aodp + 31.0 * s4) * c1;
%     t3cof = d2 + 2.0*c1sq;
%     t4cof = 0.25 * (3.0 * d3 + c1 * (12.0 * d2 + 10.0 * c1sq));
%     t5cof = 0.2 * (3.0 * d4 + 12.0 * c1 * d3 + 6.0 * d2 * d2 + 15.0 * c1sq * (2.0 * d2 + c1sq));
% end
ind = find( isimp == 0 );
    c1sq(ind) = (c1(ind).^2);
    d2(ind) = 4.0 .* aodp(ind) .* tsi(ind) .* c1sq(ind);
    temp(ind) = d2(ind) .* tsi(ind) .* c1(ind) ./ 3.0;
    d3(ind) = (17.0 .* aodp(ind) + s4(ind)).*temp(ind);
    d4(ind) = 0.5 .* temp(ind) .* aodp(ind) .* tsi(ind) .* (221.0 .* aodp(ind) + 31.0 .* s4(ind)) .* c1(ind);
    t3cof(ind) = d2(ind) + 2.0.*c1sq(ind);
    t4cof(ind) = 0.25 .* (3.0 .* d3(ind) + c1(ind) .* (12.0 .* d2(ind) + 10.0 .* c1sq(ind)));
    t5cof(ind) = 0.2 .* (3.0 .* d4(ind) + 12.0 .* c1(ind) .* d3(ind) + 6.0 .* d2(ind) .* d2(ind) + 15.0 .* c1sq(ind) .* (2.0 .* d2(ind) + c1sq(ind)));

% Update for secular gravity and atmospheric drag.
xmdf = xmo + xmdot .* tsince;
omgadf = omegao + omgdot .* tsince;
xnoddf = xnodeo + xnodot .* tsince;
omega = omgadf;
xmp = xmdf;
tsq = (tsince.^2);
xnode = xnoddf + xnodcf .* tsq;
tempa = 1.0 - c1 .* tsince;
tempe = (bstar) .* c4 .* tsince;
templ = t2cof .* tsq;
% if (isimp == 0)
%     delomg = omgcof * tsince;
%     delm = xmcof*(((1.0 + eta * cos(xmdf))^ 3.0) - delmo);
%     temp = delomg + delm;
%     xmp = xmdf + temp;
%     omega = omgadf - temp;
%     tcube = tsq * tsince;
%     tfour = tsince * tcube;
%     tempa = tempa - d2 * tsq - d3 * tcube - d4 * tfour;
%     tempe = tempe + (satdata.bstar) * c5 * (sin(xmp) - sinmo);
%     templ = templ + t3cof * tcube + tfour * (t4cof + tsince * t5cof);
% end
ind = find(isimp == 0);
    delomg(ind) = omgcof(ind) .* tsince;
    delm(ind) = xmcof(ind).*(((1.0 + eta(ind) .* cos(xmdf(ind))).^ 3.0) - delmo(ind));
    temp(ind) = delomg(ind) + delm(ind);
    xmp(ind) = xmdf(ind) + temp(ind);
    omega(ind) = omgadf(ind) - temp(ind);
    tcube = tsq .* tsince;
    tfour = tsince .* tcube;
    tempa(ind) = tempa(ind) - d2(ind) .* tsq - d3(ind) .* tcube - d4(ind) .* tfour;
    tempe(ind) = tempe(ind) + (bstar(ind)) .* c5(ind) .* (sin(xmp(ind)) - sinmo(ind));
    templ(ind) = templ(ind) + t3cof(ind) .* tcube + tfour .* (t4cof(ind) + tsince .* t5cof(ind));

a = aodp .* (tempa.^2);
e = (eo) - tempe;
xl = xmp + omega + xnode + xnodp.*templ;
beta = sqrt(1.0 - (e.^2));
xn = xke ./ (a.^1.5);
% Long period periodics
axn = e .* cos(omega);
temp = 1.0 ./ (a .* (beta.^2));
xll = temp .* xlcof .* axn;
aynl = temp .* aycof;
xlt = xl + xll;
ayn = e .* sin(omega) + aynl;
% Solve Kepler's Equation
capu = fmod2p(xlt - xnode);
temp2 = capu;
i=1;
% while(1)
%     sinepw = sin(temp2);
%     cosepw = cos(temp2);
%     temp3 = axn * sinepw;
%     temp4 = ayn * cosepw;
%     temp5 = axn * cosepw;
%     temp6 = ayn * sinepw;
%     epw = (capu - temp4 + temp3 - temp2) / (1.0 - temp5 - temp6) + temp2;
%     temp7 = temp2;
%     temp2 = epw;
%     i = i+1;
% 	if ((i>10) || (abs(epw - temp7) <= e6a))
% 		break
% 	end
% end
ind = 1:length(xno);
while(1)
    sinepw(ind) = sin(temp2(ind));
    cosepw(ind) = cos(temp2(ind));
    temp3(ind) = axn(ind) .* sinepw(ind);
    temp4(ind) = ayn(ind) .* cosepw(ind);
    temp5(ind) = axn(ind) .* cosepw(ind);
    temp6(ind) = ayn(ind) .* sinepw(ind);
    epw(ind) = (capu(ind) - temp4(ind) + temp3(ind) - temp2(ind)) ./ (1.0 - temp5(ind) - temp6(ind)) + temp2(ind);
    temp7(ind) = temp2(ind);
    temp2(ind) = epw(ind);
    i = i+1;
    ind = find((abs(epw - temp7) > e6a));
	if ((i>10) || isempty(ind))%(abs(epw - temp7) <= e6a))
		break
	end
end
% Short period preliminary quantities
ecose = temp5 + temp6;
esine = temp3 - temp4;
elsq = (axn.^2) + (ayn.^2);
temp = 1.0 - elsq;
pl = a .* temp;
r = a .* (1.0 - ecose);
temp1 = 1.0 ./ r;
rdot = xke .* sqrt(a) .* esine .* temp1;
rfdot = xke .* sqrt(pl) .* temp1;
temp2 = a .* temp1;
betal = sqrt(temp);
temp3 = 1.0 ./ (1.0 + betal);
cosu = temp2 .* (cosepw - axn + ayn .* esine .* temp3);
sinu = temp2 .* (sinepw - ayn - axn .* esine .* temp3);
u = actan(sinu, cosu);
sin2u = 2.0 .* sinu .* cosu;
cos2u = 2.0 .* (cosu.^2) - 1.0;
temp = 1.0 ./ pl;
temp1 = CK2 .* temp;
temp2 = temp1 .* temp;
% Update for short periodics
rk = r .* (1.0 - 1.5 .* temp2 .* betal .* x3thm1) + 0.5 .* temp1 .* x1mth2 .* cos2u;
uk = u - 0.25 .* temp2 .* x7thm1 .* sin2u;
xnodek = xnode + 1.5 .* temp2 .* cosio .* sin2u;
xinck = (xincl) + 1.5 .* temp2 .* cosio .* sinio .* cos2u;
rdotk = rdot - xn .* temp1 .* x1mth2 .* sin2u;
rfdotk = rfdot + xn .* temp1 .* (x1mth2 .* cos2u + 1.5 .* x3thm1);
% Orientation vectors
MV.v(1,:) = -sin(xnodek) .* cos(xinck);
MV.v(2,:) = cos(xnodek) .* cos(xinck);
MV.v(3,:) = sin(xinck);

NV.v(1,:) = cos(xnodek);
NV.v(2,:) = sin(xnodek);
NV.v(3,:) = 0;

for i=1:3
	UV.v(i,:) = MV.v(i,:) .* sin(uk) + NV.v(i,:) .* cos(uk);
	VV.v(i,:) = MV.v(i,:) .* cos(uk) - NV.v(i,:) .* sin(uk);
end

% position + velocity
for i=1:3
	posv(i,:) = rk .* UV.v(i,:);
	velv(i,:) = rdotk .* UV.v(i,:) + rfdotk .* VV.v(i,:);
end

[pos, vel] = Convert_Sat_State(posv, velv);

if isstruct(PARAMS)
    
t = cspice_et2utc(et,'ISOC',8);

tline = t;%fgetl(fid);
year = str2double(tline(1:4));
month = str2double(tline(6:7));
day = str2double(tline(9:10));
hour = str2double(tline(12:13));
minute = str2double(tline(15:16));
sec = str2double(tline(18:end));
MJD_UTC = Mjday(year, month, day, hour, minute, sec);

[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(PARAMS,MJD_UTC,'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
MJD_UT1 = MJD_UTC + UT1_UTC/86400;
MJD_TT  = MJD_UTC + TT_UTC/86400;
T = (MJD_TT-51544.5)/36525;


[reci, veci] = teme2eciSGP(pos,vel,T,dpsi,deps);

else

reci = PARAMS * pos;
veci = PARAMS * vel;

end

x = [reci;veci];


end
