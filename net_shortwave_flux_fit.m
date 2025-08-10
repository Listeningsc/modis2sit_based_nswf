function net_s = net_shortwave_flux_fit(d, month)
% 计算净短波辐射（Maykut 1986 表5.6），d 可为标量或向量（光学厚度）

    % 1) 零辐射月份
    if ismember(month, [1 2 7 8 11 12])
        net_s = zeros(size(d));
        return
    end

    % 2) 光学厚度与对应月份观测 y
    x = [0; 0.05; 0.10; 0.20; 0.40; 0.80; 3.00];
    switch month
        case 3,  y = [  7;   5;   4;   4;   4;   3;  1];
        case 4,  y = [ 83;  56;  52;  49;  46;  42; 17];
        case 5,  y = [209; 141; 131; 124; 114; 104; 42];
        case 6,  y = [281; 189; 175; 166; 153; 140; 59];
        case 9,  y = [ 89;  60;  56;  53;  48;  45; 16];
        case 10, y = [ 24;  16;  15;  14;  13;  12;  4];
        otherwise
            error('net_shortwave_flux_fit:UnsupportedMonth', ...
                  'No data for month %d.', month);
    end

    % 3) 线性最小二乘拟合 y ≈ [exp(-x)  1] * [a; b]
    U = exp(-x);
    theta = [U, ones(size(U))] \ y;   % \ 是最小二乘解
    a = theta(1); 
    b = theta(2);

    % 4) 计算输出（支持 d 向量）
    net_s = a .* exp(-d) + b;
end
