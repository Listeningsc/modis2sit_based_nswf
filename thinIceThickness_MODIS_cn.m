function [hi_corr_m, hi_raw_m, meta] = thinIceThickness_MODIS(Ts, Ta, u10, sza, Ls, opts)
% 基于 Rudjord 等 (2022) 论文 “Putting it all together” 步骤的薄冰厚度反演
% 输入参数：
%   Ts  : 冰/雪表面温度 [K]（来自 MODIS 反演）
%   Ta  : 2 米高度空气温度 [K]（来自 NCEP/NCAR 再分析）
%   u10 : 10 米高度风速 [m/s]（来自 NCEP/NCAR 再分析）
%   sza : 太阳天顶角 [度]
%   Ls  : 近岸因子（0~1，高斯平滑后的陆地掩膜）
%   opts: 可选结构体，包括：
%         .RH           相对湿度（默认 0.9）
%         .Sw           海水盐度 [psu]（默认 34.5）
%         .applyCalib   是否应用整体定标因子 k=0.171（默认 true）
%
% 输出参数：
%   hi_corr_m : 定标后的冰厚 [m]
%   hi_raw_m  : 未定标的冰厚 [m]
%   meta      : 中间变量（可选输出）

% --------------------- 参数初始化 ---------------------
if nargin < 6, opts = struct; end
if ~isfield(opts,'RH'), opts.RH = 0.90; end         % 相对湿度（论文假设值）
if ~isfield(opts,'Sw'), opts.Sw = 34.5; end         % 海水盐度 [psu]
if ~isfield(opts,'applyCalib'), opts.applyCalib = true; end

% --------------------- 常数定义（与论文一致） ---------------------
sigma = 5.6704e-8;     % 斯特藩-玻尔兹曼常数 [W m^-2 K^-4] (式15-16)
eps_s = 0.988;         % 表面发射率（冰/雪）(式15)
Pa    = 1013.25;       % 空气压强 [hPa] (式20, 23)
Rgas  = 287.1;         % 干空气气体常数 [J kg^-1 K^-1] (式20)
L_v   = 2.5e6;         % 水汽潜热 [J kg^-1] (式23)
S0    = 1361;          % 太阳常数 [W m^-2] (式3)
qs_sat= 0.003767;      % 饱和比湿 (式21-22)
e_star= 0.7855*(1 + 0.2232*0^2.75); % 有效大气发射率 (式17, C=0)
ks    = 0.31;          % 雪的导热率 [W m^-1 K^-1] (式7)
k0    = 2.034;         % 纯冰导热率 (式8)
betaK = 0.13;          % 导热率与盐度关系系数 (式8)
T0    = 273.15;        % 纯水冰点 [K] (式8)

% --------------------- 近岸系数 (式13) ---------------------
L = 2.5.*Ls + 1.0;

% --------------------- 海水冰点 (Tf = -0.055*Sw) ---------------------
Tf = 273.15 + (-0.055*opts.Sw);

% --------------------- 湿空气参数 ---------------------
f  = opts.RH;                         % 相对湿度
q  = f*qs_sat;                        % 比湿
Tv = (1 + 0.608*q).*Ta;               % 虚温 (式21)
rho_a = (100.*Pa) ./ (Rgas.*Tv);       % 空气密度 (式20)
cp = 1004.5*(1.0 + 0.9433*q);          % 湿空气定压比热 (式22)

% --------------------- 风速与交换系数 ---------------------
u2 = u10 ./ 1.27;  % 10米风速换算到2米高度 (式19)
% Bentamy (2003) 潜热交换系数 (式24)
aB=-0.146785; bB=-0.292400; cB=-2.206648; dB=1.6112292;
Ce = (aB.*exp(bB.*(u10 + cB)) + dB./(u10 + 1)).*1e-3;
Cs = 0.98.*Ce; % 感热交换系数 (式18)

% --------------------- 饱和水汽压计算 (式25-27) ---------------------
poly_esat = @(T) ( 2.7798202e-6*T.^4 ...
                 - 2.6913393e-3*T.^3 ...
                 + 0.979208*T.^2 ...
                 - 158.6377*T ...
                 + 9653.1925 );
es0 = poly_esat(Ts);  % 表面
esa = poly_esat(Ta);  % 2米高度
ea  = f.*esa;         % 空气实际水汽压
e0  = f.*es0;         % 表面水汽压

% --------------------- 短波入射辐射 (Shine 1984, 式3) ---------------------
mu = cosd(sza);
Finsw = S0.*( mu.^2 ./ (1.2*mu + (1.0 + mu).*1e-3.*e0 + 0.0455) );

% --------------------- Grenfell (1979) 反照率与透过率 ---------------------
    function alb = grenfell_albedo(hi_cm, hs_m)
        if hs_m <= 0
            A1=0.130; B1=15.46; C1=0.820; D1=0.1216; % 裸冰 (式4)
        else
            A1=0.2213; B1=77.48; C1=0.198; D1=0.0;   % 有雪冰 (式4)
        end
        h = hs_m>0 && hs_m || hi_cm/100; % 根据雪深或冰厚选输入
        alb = 1 - A1*exp(-B1*h) - C1*exp(-D1*h);
    end

    function tr = grenfell_i0(hi_m, hs_m)
        if hs_m <= 0
            A2=0.1925; B2=12.96; C2=0.515; D2=1.227; % 裸冰 (式5)
        else
            hs = max(hs_m,0);
            A2 = 0.2257*exp(-16.73*hs) + 0.4174*exp(-43.89*hs);
            B2 = 0.7280*exp(-0.1862*hs) + 0.3532*exp(-13.04*hs);
            C2 = 0.1561*exp(-92.79*hs);
            D2 = 1./(0.06 + 0.0995*exp(-94.20*hs));
        end
        tr = A2.*exp(-B2*hi_m) - C2.*exp(-D2*hi_m);
    end

% --------------------- 雪深分段线性模型 (式14) ---------------------
    function [a_s, b_s] = snow_linear_coeffs(h_init_cm, Lfac)
        if h_init_cm <= 5
            a_s = 0.0; b_s = 0.0;
        elseif h_init_cm < 20
            a_s = 0.05*Lfac; b_s = -0.0025*Lfac;
        else
            a_s = 0.1*Lfac; b_s = -0.0125*Lfac;
        end
        a_s = a_s/100; % cm -> m
    end

% --------------------- 盐度分段线性模型 (式9-10) ---------------------
    function [s0,s1] = salinity_linear(hi_init_cm)
        if hi_init_cm <= 40
            s0 = 12.24; s1 = -19.39/100;
        else
            s0 = 7.88;  s1 = -1.59/100;
        end
    end

% --------------------- 冰导热率 (式8) ---------------------
    function ki_val = ice_k(S_psu, TsK)
        dT = max(abs(TsK - T0), 0.5);
        ki_val = k0 + betaK * S_psu ./ dT;
        ki_val = min(max(ki_val, 1.6), 3.5);
    end

% --------------------- 感热和潜热的基础通量 (式15-18, 23) ---------------------
Fupl = eps_s .* sigma .* Ts.^4;
Fdnl = e_star .* sigma .* Ta.^4;
Fs = rho_a .* cp .* Cs .* u2 .* (Ta - Ts);
Fe = rho_a .* L_v .* Ce .* u2 .* (ea - es0) .* (0.622 ./ Pa);

% --------------------- 七个初始猜测冰厚（cm） ---------------------
h_init_list = [2.5, 7.5, 15, 30, 50, 70, 90];
sz = size(Ts);
h_cand_all = nan([sz, numel(h_init_list)]);

% --------------------- 遍历每个初始值，计算候选冰厚 ---------------------
for k = 1:numel(h_init_list)
    hi0 = h_init_list(k);
    [a_s, b_s] = snow_linear_coeffs(hi0, L);
    hs0 = a_s.*hi0 + b_s; % m
    alb = arrayfun(@(hs) grenfell_albedo(hi0, hs), hs0);
    i0  = arrayfun(@(hs) grenfell_i0(hi0/100, hs), hs0);
    Fsw = (1 - alb).*(1 - i0).*Finsw;
    Fc_target = -(Fsw - Fupl + Fdnl + Fs + Fe);
    [s0,s1] = salinity_linear(hi0);
    S_eval  = s0 + s1*hi0;
    ki_eval = ice_k(S_eval, Ts);
    Knum = ki_eval .* ks .* (Tf - Ts);
    Acoef = (ks/100) + ki_eval.*a_s;
    Bcoef = ki_eval.*b_s;
    hi_cm = ( (Knum ./ Fc_target) - Bcoef ) ./ Acoef;
    hi_cm = max(hi_cm, 0);
    h_cand_all(:,:,k) = hi_cm;
end

% --------------------- 选择最接近自身初值的解 ---------------------
D = abs(h_cand_all - reshape(h_init_list,1,1,[]));
[~, idx] = min(D, [], 3);
linIdx = sub2ind([sz, numel(h_init_list)], repmat((1:sz(1))',1,sz(2)), repmat(1:sz(2),sz(1),1), idx);
hi_sel_cm = h_cand_all(linIdx);

% --------------------- 输出冰厚（m） ---------------------
hi_raw_m  = hi_sel_cm / 100;
if opts.applyCalib
    kcal = 0.171; % 定标因子
    hi_corr_m = hi_raw_m / kcal;
else
    hi_corr_m = hi_raw_m;
end

% --------------------- 输出中间变量（可选） ---------------------
if nargout > 2
    meta.Fupline = Fupl; meta.Fdnl = Fdnl; meta.Fs = Fs; meta.Fe = Fe; 
    meta.Finsw = Finsw; meta.L = L; meta.u2 = u2; meta.Ce = Ce; meta.Cs = Cs;
end
end
