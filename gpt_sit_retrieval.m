function [hi_corr_m, hi_raw_m, meta] = gpt_sit_retrieval(Ts, Ta, u10, sza, Ls, opts)
% Implements the "Putting it all together" step from Rudjord et al. (2022)
% Inputs: Ts, Ta [K], u10 [m/s], sza [deg], Ls [0..1]; same-size arrays
% opts: struct with fields (optional):
%   .RH (default 0.9), .Sw (seawater salinity, default 34.5), .applyCalib (true)
% Outputs:
%   hi_corr_m : calibrated ice thickness [m]
%   hi_raw_m  : uncalibrated thickness [m]
%   meta      : struct of intermediate fields (optional use)

% ---- defaults
if nargin < 6, opts = struct; end
if ~isfield(opts,'RH'), opts.RH = 0.90; end               % relative humidity (f)  (Key/Yu; used across terms)
if ~isfield(opts,'Sw'), opts.Sw = 34.5; end               % seawater salinity for Tf   (Tf = -0.055*Sw)  [psu]
if ~isfield(opts,'applyCalib'), opts.applyCalib = true; end

% ---- constants (from paper)
sigma = 5.6704e-8;                                        % Stefan每Boltzmann [W m^-2 K^-4]  (Eq. 15每16)
eps_s = 0.988;                                            % surface emissivity (snow/ice)     (Eq. 15)
Pa    = 1013.25;                                          % air pressure [hPa]               (Eq. 20, 23)
Rgas  = 287.1;                                            % [J kg^-1 K^-1]                   (Eq. 20)
L_v   = 2.5e6;                                            % latent heat [J kg^-1]            (Eq. 23)
S0    = 1361;                                             % solar constant [W m^-2]          (Eq. 3)
qs_sat= 0.003767;                                         % saturated specific humidity      (text under Eq. 21每22)
e_star= 0.7855*(1 + 0.2232*0^2.75);                       % clear-sky effective emissivity   (Eq. 17; C=0)
ks    = 0.31;                                             % snow thermal conductivity [W m^-1 K^-1] (text near Eq. 7)
k0    = 2.034;                                            % pure ice conductivity            (Eq. 8)
betaK = 0.13;                                             % coefficient in ki(S,Ti)          (Eq. 8)
T0    = 273.15;                                           % pure water freezing [K]          (Eq. 8)

% land proximity factor (Eq. 13)
L = 2.5.*Ls + 1.0;

% freezing point of seawater (used on ocean side of ice)  Tf = -0.055*Sw  (～C) => [K]
Tf = 273.15 + (-0.055*opts.Sw);

% RH-related quantities
f  = opts.RH;
q  = f*qs_sat;                                            % specific humidity (approx)       (Eq. 21每22)
Tv = (1 + 0.608*q).*Ta;                                   % virtual temperature              (Eq. 21)
rho_a = (100.*Pa) ./ (Rgas.*Tv);                          % air density                      (Eq. 20)
cp = 1004.5*(1.0 + 0.9433*q);                             % moist air cp                     (Eq. 22)

% 10m->2m wind, bulk coefficients (Eq. 19, 24; Bentamy 2003)
u2 = u10 ./ 1.27;                                         % (Eq. 19)
aB=-0.146785; bB=-0.292400; cB=-2.206648; dB=1.6112292;   % (Eq. 24 coefficients)
Ce = (aB.*exp(bB.*(u10 + cB)) + dB./(u10 + 1)).*1e-3;     % guard small u10 with +1 (paper has +1 in denom) (Eq. 24, bracketed)
Cs = 0.98.*Ce;                                            % (text after Eq. 19)

% Saturation vapor pressure poly (Maykut 1982; Eq. 25)
poly_esat = @(T) ( 2.7798202e-6*T.^4 ...
                 - 2.6913393e-3*T.^3 ...
                 + 0.979208*T.^2 ...
                 - 158.6377*T ...
                 + 9653.1925 );
es0 = poly_esat(Ts);                                      % saturation at surface Ts (Eq. 26)
esa = poly_esat(Ta);                                      % saturation at 2 m Ta   (Eq. 27)
ea  = f.*esa;                                             % actual air vapor pressure (text below Eq. 27)
e0  = f.*es0;                                             % surface vapor pressure  (used in Shine shortwave term)

% Shortwave incoming flux Finsw (Shine 1984; Eq. 3)
mu = cosd(sza);
Finsw = S0.*( mu.^2 ./ (1.2*mu + (1.0 + mu).*1e-3.*e0 + 0.0455) );

% Helper: Grenfell albedo & transmittance parameterizations
    function alb = grenfell_albedo(hi_cm, hs_m)
        % if no snow -> "bare ice" coeffs depend on ice thickness (in meters)
        % if snow -> "snow-covered" coeffs depend on snow depth
        if hs_m <= 0
            A1=0.130; B1=15.46; C1=0.820; D1=0.1216;   % bare ice (Eq. 4 text)
            h = max(hi_cm/100, 0);                      % convert cm->m for exponent scale
        else
            A1=0.2213; B1=77.48; C1=0.198; D1=0.0;     % snow-covered (Eq. 4 text)
            h = max(hs_m, 0);
        end
        alb = 1 - A1*exp(-B1*h) - C1*exp(-D1*h);        % (Eq. 4)
    end

    function tr = grenfell_i0(hi_m, hs_m)
        if hs_m <= 0
            % bare ice constants (Eq. 5 text)
            A2=0.1925; B2=12.96; C2=0.515; D2=1.227;
        else
            % snow-covered: A2..D2 are functions of hs (Eq. 6)
            hs = max(hs_m,0);
            A2 = 0.2257*exp(-16.73*hs) + 0.4174*exp(-43.89*hs);
            B2 = 0.7280*exp(-0.1862*hs) + 0.3532*exp(-13.04*hs);
            C2 = 0.1561*exp(-92.79*hs);
            D2 = 1./(0.06 + 0.0995*exp(-94.20*hs));
        end
        tr = A2.*exp(-B2*hi_m) - C2.*exp(-D2*hi_m);     % (Eq. 5)
    end

% Snow depth piecewise-linear with land factor (continuous form, Eq. 14; L from Eq. 13)
    function [a_s, b_s] = snow_linear_coeffs(h_init_cm, Lfac)
        hi = h_init_cm;
        if hi <= 5
            a_s = 0.0;             b_s = 0.0;                           % hs = 0
        elseif hi < 20
            a_s = 0.05*Lfac;       b_s = -0.0025*Lfac;                  % hs = (0.05*hi - 0.0025)*L
        else
            a_s = 0.1*Lfac;        b_s = -0.0125*Lfac;                  % hs = (0.1*hi - 0.0125)*L
        end
        % units: hi in cm; hs in m  (a_s has units m/cm, b_s in m)
        a_s = a_s/100; % convert to m per cm of ice
    end

% Sea-ice salinity piecewise-linear (Cox & Weeks 1974; Eqs. 9每10)  S = s0 + s1*hi_cm
    function [s0,s1] = salinity_linear(hi_init_cm)
        if hi_init_cm <= 40
            s0 = 12.24; s1 = -19.39/100;     % hi in cm -> per cm slope
        else
            s0 = 7.88;  s1 = -1.59/100;
        end
    end

% thermal conductivity of ice ki (Untersteiner 1964; Eq. 8) using top-layer Ti＞Ts
    function ki_val = ice_k(S_psu, TsK)
        % Guard vs division by small (Ts-T0)
        dT = max(abs(TsK - T0), 0.5);          % avoid blow-up; paper notes weak dependence below -4～C
        ki_val = k0 + betaK * S_psu ./ dT;
        % Optional clamp to reasonable range
        ki_val = min(max(ki_val, 1.6), 3.5);
    end

% Saturation vapor pressure helper for shortwave term already defined as poly_esat

% Core: compute flux pieces that do NOT depend on hi (for given albedo/i0 they do, so those go in loop)
Fupl = eps_s .* sigma .* Ts.^4;                 % (Eq. 15)
Fdnl = e_star .* sigma .* Ta.^4;                % (Eq. 16)
% Sensible & latent need Cs/Ce,u2, rho_a, cp, Ta, Ts, ea, es0
Fs = rho_a .* cp .* Cs .* u2 .* (Ta - Ts);      % (Eq. 18)
Fe = rho_a .* L_v .* Ce .* u2 .* (ea - es0) .* (0.622 ./ Pa);   % (Eq. 23)

% Seven initial guesses in cm (paper)
h_init_list = [2.5, 7.5, 15, 30, 50, 70, 90];

% Allocate arrays
sz = size(Ts);
h_cand_all = nan([sz, numel(h_init_list)]);

% Loop over initial guesses to determine segments & solve analytic hi
for k = 1:numel(h_init_list)
    hi0 = h_init_list(k);                            % cm
    % piecewise coeffs for hs = a*hi + b (in meters), where hi in cm
    [a_s, b_s] = snow_linear_coeffs(hi0, L);
    % build a provisional hs from hi0 just to choose albedo/transmittance regime
    hs0 = a_s.*hi0 + b_s;                            % meters
    % choose albedo & i0 based on regime (bare vs snow)
    % NOTE: albedo uses hi in cm (converted inside), transmittance uses hi in m
    % We'll use a provisional hi (convert to m) to compute albedo & i0 feedback terms
    alb = arrayfun(@(hs) grenfell_albedo(hi0, hs), hs0);  % scalar per pixel regime
    i0  = arrayfun(@(hs) grenfell_i0(hi0/100, hs), hs0);
    % Shortwave absorbed (Eq. 2): Fsw = (1 - a)(1 - i0) Finsw
    Fsw = (1 - alb).*(1 - i0).*Finsw;

    % Target conductive flux from energy balance (Eq. 1): Fc = -(Fsw - Fupl + Fdnl + Fs + Fe)
    Fc_target = -(Fsw - Fupl + Fdnl + Fs + Fe);

    % Salinity linear segment (S = s0 + s1*hi_cm); we evaluate ki at S(hi0) as constant within this branch
    [s0,s1] = salinity_linear(hi0);
    S_eval  = s0 + s1*hi0;                            % psu at initial guess
    ki_eval = ice_k(S_eval, Ts);                      % W/m/K  (use top layer Ts)

    % Analytic solution for hi (cm) with hs = a*hi + b:
    % Fc = ki*ks*(Tf - Ts) / (ks*hi_m + ki*hs)  with hi_m in meters, hs=a*hi_cm + b
    % Denominator = ks*(hi_cm/100) + ki*(a*hi_cm + b) = (ks/100 + ki*a)*hi_cm + ki*b
    Knum = ki_eval .* ks .* (Tf - Ts);               % numerator
    Acoef = (ks/100) + ki_eval.*a_s;
    Bcoef = ki_eval.*b_s;
    hi_cm = ( (Knum ./ Fc_target) - Bcoef ) ./ Acoef; % cm  〞〞  derived from energy balance rearrangement

    % keep only physically meaningful (positive) and finite
    hi_cm = max(hi_cm, 0);
    h_cand_all(:,:,k) = hi_cm;
end

% Select candidate thickness closest to its own initial guess (paper's selection rule)
D = abs(h_cand_all - reshape(h_init_list,1,1,[]));
[~, idx] = min(D, [], 3);
linIdx = sub2ind([sz, numel(h_init_list)], repmat((1:sz(1))',1,sz(2)), repmat(1:sz(2),sz(1),1), idx);
hi_sel_cm = h_cand_all(linIdx);

% raw and calibrated outputs
hi_raw_m  = hi_sel_cm / 100;                % meters
if opts.applyCalib
    kcal = 0.171;                           % paper's scale factor (hi_corr = hi / k) 
    hi_corr_m = hi_raw_m / kcal;
else
    hi_corr_m = hi_raw_m;
end

% pack some helpful meta if needed
if nargout > 2
    meta.Fupline = Fupl; meta.Fdnl = Fdnl; meta.Fs = Fs; meta.Fe = Fe; 
    meta.Finsw = Finsw; meta.L = L; meta.u2 = u2; meta.Ce = Ce; meta.Cs = Cs;
end
end

