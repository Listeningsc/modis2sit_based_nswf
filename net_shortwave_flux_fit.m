
function net_s = net_shortwave_flux_fit(d, month)
% net_short  Compute net shortwave radiation per Maykut (1986) Table 5.6
% 针对文献参数进行拟合得出短波辐射值，与短波参数化结果进行对比
%   d.  Months 1,2,7,8,11,12 always 0.

    % months with zero net shortwave
    no_net_shortwave_flux_month = [1, 2, 7, 8, 11, 12];

    % optical depths
    x = [0; 0.05; 0.1; 0.2; 0.4; 0.8; 3.0];

    % check zero‐radiation months
    if ismember(month, no_net_shortwave_flux_month)
        net_s = 0;
        return
    end

    % select the measured y-values for this month
    switch month
        case 3
            y = [7;   5;  4;  4;  4;  3;  1];
        case 4
            y = [83; 56; 52; 49; 46; 42; 17];
        case 5
            y = [209;141;131;124;114;104;42];
        case 6
            y = [281;189;175;166;153;140;59];
        case 9
            y = [89; 60; 56; 53; 48; 45; 16];
        case 10
            y = [24; 16; 15; 14; 13; 12;  4];
        otherwise
            error('net_short:UnsupportedMonth', ...
                  'No data defined for month %d.', month);
    end

    % initial guesses: amplitude = max(y)-min(y), offset = min(y)
    a0 = max(y) - min(y);
    b0 = min(y);

    % set up fit type and options
    ft   = fittype(@tar_func, 'independent', 'x', 'coefficients', {'a','b'});
    opts = fitoptions('Method', 'NonlinearLeastSquares', ...
                      'StartPoint', [a0, b0]);

    % perform the fit
    cf   = fit(x, y, ft, opts);

    % evaluate at depth d
    net_s = cf.a * exp(-d) + cf.b;
end
% tar_func.m
function y = tar_func(x, a, b)
% tar_func   Exponential fit function: a·exp(−x) + b
%   y = tar_func(x, a, b) returns a*exp(-x) + b, vectorized in x.
    y = a .* exp(-x) + b;
end

