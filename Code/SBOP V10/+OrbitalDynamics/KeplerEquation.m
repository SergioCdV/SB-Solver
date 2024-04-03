
function [dt] = KeplerEquation(n, e, nu_0, nu_f)
    % Initial mean anomaly 
    cos_E = (e + cos(nu_0)) / (1 + e * cos(nu_0));
    sin_E = (sqrt(1-e^2)*sin(nu_0)) / (1 + e * cos(nu_0));
    E = atan2(sin_E, cos_E);
    M0 = E-e*sin(E);
    M0 = mod(M0,2*pi);

    % Final mean anomaly 
    cos_E = (e + cos(nu_f)) / (1 + e * cos(nu_f));
    sin_E = (sqrt(1-e^2)*sin(nu_f)) / (1 + e * cos(nu_f));
    E = atan2(sin_E, cos_E);
    Mf = E-e*sin(E);
    Mf = mod(Mf,2*pi);

    % Time step 
    dM = Mf-M0;
    if (dM < 0)
        Mf = Mf + 2 * pi;
        dM = Mf-M0;
    end
    dt = dM/n;
end