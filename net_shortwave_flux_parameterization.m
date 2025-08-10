function nswhf = net_shortwave_flux_parameterization(solar_zenith_angle,sit)
    % 此为日间 “净短波辐射通量(nswhf)”的参数化计算方法
    % 该计算方法取自 Shine(1984)，以及Grenfell(1979)
    % 具体公式为Fsw = (1-α)(1-i)F[in]sw  其中α为反照率  i为下垫面透过率 F[in]sw为入射短波辐射通量
    if sit<0.05
        % 积雪不存在，所以snow=0
        snow = 0; 
        snow_d = 0;
    elseif sit>=0.05 && sit<=0.2
        % 积雪存在，所以snow=1
        snow = 1;
        snow_d = 0.05*sit;
    else
        snow = 1;
        snow_d = 0.1*sit;
    end
     
    %% 该块为入射短波辐射通量计算
    S0 = 1361;
    % S0为太阳入射常数 单位Wm-2
    u = cosd(solar_zenith_angle);
    % 太阳天顶角的余弦
    f = 0.9; % 相对湿度，这是经验常数
    es0 = 3;
    e0 = f*es0;
    F_in_sw = S0*u^2/(1.2*u+10e-3*(1+u)*e0+0.0455);
    %% 该块为反照率计算，和穿透率计算
    if snow
        % X1的表示都是计算反照率的参数
        % X2的表示都是计算穿透率的参数
        h1 = snow_d;
        h2 = snow_d;
        A1 = 0.2213;
        A2 = 0.2257*exp(-16.73*h2)+0.4174*exp(-43.89*h2);
        B1 = 77.48;
        B2 = 0.728*exp(-0.1862*h2)+0.3532*exp(-13.04*h2);
        C1 = 0.198;
        C2 = 0.1561*exp(-92.79*h2);
        D1 = 0;
        D2 = 1/(0.06+0.0995*exp(-94.0*h2));
    else
        A1 = 0.13;
        A2 = 0.1925;
        B1 = 15.46;
        B2 = 12.96;
        C1 = 0.82;
        C2 = 0.515;
        D1 = 0.1216;
        D2 = 1.227;
        h1 = sit;
        h2 = sit;
    end
    a = 1-A1*exp(-B1*h1)-C1*exp(-D1*h1);
    i = A2*exp(-B2*h2)-C2*exp(-D2*h2);
    
    nswhf = (1-a)*(1-i)*F_in_sw;
end

