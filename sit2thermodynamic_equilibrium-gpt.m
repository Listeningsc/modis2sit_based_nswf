function sit2thermodynamic_equilibrium(result_doy,ist_parent_path,ist_coord_parent_path,cloud_parent_path,sensor,night,target_year,era_path,era_mat_path,str_ib_date)
    era_output = era_data_exact(era_path,target_year,era_mat_path);
    [era_cloud_lat, era_cloud_lon, era_cloud_hcc, era_cloud_mcc, ...
     era_lon, era_lat, era_wind_u, era_wind_v, era_t2m] = era_output{:};
 
    %% modis 数据路径
    % 表面温度路径
    % day_path 表示是白天还是晚上
    if night
        day_path = 'night';
    else
        day_path = 'day';
    end
    target_ist_path = fullfile(ist_parent_path,sensor,day_path,num2str(target_year));
    % 坐标路径
    ist_coord_path = fullfile(ist_coord_parent_path,sensor,day_path,num2str(target_year));
    % 云掩膜路径
    cloud_path = fullfile(cloud_parent_path,sensor,day_path,num2str(target_year)); 
    % 温度数据
    ist_hdf_list = dir(fullfile(target_ist_path,'*.hdf'));
    % 温度1KM坐标数据
    ist_coord_list = dir(fullfile(ist_coord_path,'*.hdf'));
    %modis云掩膜数据
    cloud_hdf_list = dir(fullfile(cloud_path,'*.hdf'));
    
    % 根据得出的数据信息，进而对对应的IB数据日期进行提取
    % 每个数据天下，会有多个IB数据
    matches_ist = contains({ist_hdf_list.name}, result_doy);
    filtered_name = {ist_hdf_list(matches_ist).name}';
    target_month = str2double(str_ib_date(5:6));
    for i = 1:numel(filtered_name)
        target_ist_name = filtered_name{i};
        timestamp = target_ist_name(8:19);
        
        %era5辅助数据 现将温度的doy天数判断是否需要将中间月份剔除，并且将选取的时间点都提取出来
        target_doy = str2double(target_ist_name(12:14));
        target_hour = str2double(target_ist_name(16:17));
        if target_doy>89
            target_doy = target_doy - (30+31+30+31+31);
        end
        target_air_t2m = era_t2m(:,:,target_doy);
        target_spe_u = era_wind_u(:,:,target_doy);
        target_spe_v = era_wind_v(:,:,target_doy);
        
        wind_10m_speed = sqrt(target_spe_u.^2 + target_spe_v.^2);
        wind_2m_speed  = wind_10m_speed/1.27;
        % ===== ERA 2m 气温/风速插值器（适配 1440x84：lon在第1维，lat在第2维） =====
        Ta2D = target_air_t2m;      % 期望 size: [1440 x 84]  (维1=lon, 维2=lat)
        U22D = wind_2m_speed;       % 同上
       % （若加载为 [84 x 1440]，可加一行安全转置）
        if size(Ta2D,1)==84 && size(Ta2D,2)==1440
            Ta2D = Ta2D.';  U22D = U22D.';
        end
        
        % 正确抽轴：lon 取第1列（随行变化），lat 取第1行（随列变化）
        lon_vec = double(era_lon(:,1));   % 1440x1
        lat_vec = double(era_lat(1,:));   % 1x84
        
       % --- 经度处理：映射到 [-180,180)，并按同维度重排数据（第1维） ---
        lon_vec(lon_vec >= 180) = lon_vec(lon_vec >= 180) - 360;
        [lon_vec, lon_ord] = sort(lon_vec, 'ascend');
        Ta2D = Ta2D(lon_ord, :);
        U22D = U22D(lon_ord, :);
        % --- 纬度处理：若为降序则翻转，并按同维度重排数据（第2维） ---
        if numel(lat_vec) > 1 && lat_vec(2) < lat_vec(1)
            lat_vec = fliplr(lat_vec);
            Ta2D = Ta2D(:, end:-1:1);
            U22D = U22D(:, end:-1:1);
        end
           % --- 保险：去重（极少见，但可避免0/360重复之类情况） ---
        [lon_vec, ia_lon] = unique(lon_vec, 'stable');  Ta2D = Ta2D(ia_lon, :);  U22D = U22D(ia_lon, :);
        [lat_vec, ia_lat] = unique(lat_vec, 'stable');  Ta2D = Ta2D(:, ia_lat);  U22D = U22D(:, ia_lat);

        % 生成 griddedInterpolant（V 的尺寸必须是 [numel(lon_vec) x numel(lat_vec)]）
        F_Ta = griddedInterpolant({lon_vec, lat_vec}, Ta2D, 'linear', 'nearest');
        F_U2 = griddedInterpolant({lon_vec, lat_vec}, U22D, 'linear', 'nearest');
        
        matches_coord = contains({ist_coord_list.name}, timestamp);
        matches_cloud = contains({cloud_hdf_list.name}, timestamp);
        if sum(matches_coord)==0 || sum(matches_cloud)==0
            continue
        else
            fprintf('success in %s',timestamp);
            target_coord_name = ist_coord_list(matches_coord).name;
            target_cloud_name = cloud_hdf_list(matches_cloud).name; 

            % 注意，这里的输出参数，都是剔除过 多余列数的数据
            % 注意ist_coarse_lon 和 ist_lon 都是-180-180度的范围
            ist_output = ist_data_exact(target_ist_path,ist_coord_path,cloud_path,target_ist_name,target_coord_name,target_cloud_name);
            [mod_ist,~,~,ist_lat,ist_lon,ist4cloud,myd_x,myd_y,target_cloud_file] = ist_output{:};
            [ist_r,ist_c]= size(mod_ist);
             % 读取无云条件下的海冰表面温度
             % 得出的输出参数，第一项为根据Modis云掩膜得出的无云下的冰表面温度；第二项为根据block，云识别后得出的无云表面温度
             % surf_non_cloud_ist是n*3的数据矩阵
            [~,surf_non_cloud_ist] = ist_no_cloud_exact(mod_ist,ist4cloud,ist_lon,ist_lat,ist_r,ist_c,...
                                                        era_cloud_lon,era_cloud_lat, era_cloud_hcc,era_cloud_mcc,...
                                                            target_doy,target_hour,myd_x,myd_y);           
                                                        
           %% FIX: 空集保护
            if isempty(surf_non_cloud_ist)
                warning('No clear-sky IST for %s, skip.', timestamp);
                continue
            end
            
            % 先对夜间数据进行处理，得出夜间海冰厚度，作为计算日间海冰厚度的初始值,并再设定范围，进行迭代计算
            expanded_solar_zenith = solar_zenith_exact(target_cloud_file,mod_ist);
           %% 白天计算 得出应用的太阳天顶角的数据矩阵2030*1350
            inx_zenith = ist_lon < -120 & ist_lon > -170 & ist_lat>66 & ist_lat<80;
                        % 量化并排序
            nc = sortrows([round(surf_non_cloud_ist(:,1),4), ...
                           round(surf_non_cloud_ist(:,2),4), ...
                           surf_non_cloud_ist(:,3)], [1 2]);           % N×3

            sz = sortrows([round(ist_lon(inx_zenith),4), ...
                           round(ist_lat(inx_zenith),4), ...
                           expanded_solar_zenith(inx_zenith)], [1 2]); % M×3
            % 用 ismember 的定位索引保证一一对应与相同顺序
            [tf, loc] = ismember(nc(:,1:2), sz(:,1:2), 'rows');
            target_non_cloud_ist = nc(tf,:);         % N'×3
            target_solar_zenith  = sz(loc(tf),:);    % N'×3（与上同序）

            %% FIX: 再次空集保护
            if isempty(target_non_cloud_ist)
                warning('No matched SZA-IST pixels for %s, skip.', timestamp);
                continue
            end     
            % ---- AWI 前一天 SIT 作为初值 ----
            sit_awi = awi_exact_sit(str_ib_date);

            %% FIX: 启动并行池（只开一次，不在每景结束时关闭）
            if isempty(gcp('nocreate'))
                parpool(4);   % 让 MATLAB 自选并行度；你也可填固定数，如 parpool(4)
            end
            npts = size(target_non_cloud_ist,1);
            sit_fit   = nan(npts,1);   % 第3列
            sit_param = nan(npts,1);   % 第4列

            % 广播常量（避免每个 worker 拷贝/构造）
            const_F_Ta = parallel.pool.Constant(F_Ta);
            const_F_U2 = parallel.pool.Constant(F_U2);
            const_AWI  = parallel.pool.Constant(sit_awi);               
            constEll   = parallel.pool.Constant(wgs84Ellipsoid('km')); 
            % ---- 并行每个像元 ----
            parfor j = 1:npts
                lon_j = target_non_cloud_ist(j,1);
                lat_j = target_non_cloud_ist(j,2);
                Ts_j  = target_non_cloud_ist(j,3);   % K
                sza_j = target_solar_zenith(j,3);    % deg

                % 最近 AWI 初值
                awi  = const_AWI.Value;    % [lon, lat, sit]
                ang  = distance(lat_j, lon_j, awi(:,2), awi(:,1), constEll.Value);
                [~, j_inx]     = min(deg2km(ang));
                reference_sit  = awi(j_inx,3);

                % 基线：使用你的拟合（若改用参数化，请替换下一行）
                nswf_fit = net_shortwave_flux_fit(reference_sit, target_month);
                tmp3     = ist2sit([lon_j, lat_j, Ts_j], const_F_Ta.Value, const_F_U2.Value, nswf_fit);
                if isnan(tmp3), tmp3 = -1; end

                % 扫描参数化
                sit_max   = max(0.9, 1.5*reference_sit);
                sit_range = [0.05:0.025:0.5, 0.55:0.05:sit_max].';
                cand      = nan(numel(sit_range),1);
                for k = 1:numel(sit_range)
                    nswf_k = net_shortwave_flux_parameterization(sza_j, sit_range(k));
                    cand(k) = ist2sit([lon_j, lat_j, Ts_j], const_F_Ta.Value, const_F_U2.Value, nswf_k);
                end
                good = ~isnan(cand);
                if any(good)
                    svec = sit_range(good);  cvec = cand(good);
                    [~,best_idx] = min(abs(svec - cvec));   %% FIX: 单一 argmin，避免多元素赋值
                    tmp4 = cvec(best_idx);
                else
                    tmp4 = NaN;
                end

                sit_fit(j)   = tmp3;   % 切片写入
                sit_param(j) = tmp4;
            end

            target_non_cloud_sit = [target_non_cloud_ist, sit_fit, sit_param];
            % ---- 输出-----
            output_path = 'W:\ljx_aux\MODIS\output_sit';
             if ~exist(output_path,'dir'); mkdir(output_path); end
            output_name = [timestamp,'.mat'];
            save(fullfile(output_path,output_name),'target_non_cloud_sit');
        end
    end
end

function era_output = era_data_exact(era_path,target_year,era_mat_path)
                                                                                                        
 %%  era5 辅助数据
    % era5 不同高度云信息
    target_era_cloud = fullfile(era_path,[num2str(target_year),'_era_cloud.nc']);
    era_cloud_lat = ncread(target_era_cloud,'latitude');
    era_cloud_lon = ncread(target_era_cloud,'longitude');
    era_cloud_hcc = ncread(target_era_cloud,'hcc');
    era_cloud_mcc = ncread(target_era_cloud,'mcc');

    era_lon = load([era_mat_path,[num2str(target_year),'_lon.mat']]).mesh_lon;
    era_lat = load([era_mat_path,[num2str(target_year),'_lat.mat']]).mesh_lat;
    era_wind_u = load([era_mat_path,[num2str(target_year),'_wind_u.mat']]).wind_u;
    era_wind_v = load([era_mat_path,[num2str(target_year),'_wind_v.mat']]).wind_v;
    era_t2m = load([era_mat_path,[num2str(target_year),'_t2m.mat']]).t2m;
                        
    era_output = {era_cloud_lat,era_cloud_lon,era_cloud_hcc,era_cloud_mcc,era_lon,era_lat,era_wind_u,era_wind_v,era_t2m};
end

function ist_output = ist_data_exact(target_ist_path,ist_coord_path,cloud_path,target_ist_name,target_coord_name,target_cloud_name)
    proj = projcrs(3413);  % WGS84北极极射投影
    target_ist_file = fullfile(target_ist_path,target_ist_name);
    target_coord_file = fullfile(ist_coord_path,target_coord_name);
    target_cloud_file = fullfile(cloud_path,target_cloud_name);

    % 读取ist表面温度 1km分辨率, 
    % 注意mod_ist的格网个数是2030*1354/2030*1350，这是由modis的不同位置决定的，再后边需要对1354的列进行删除最后四列
    mod_ist = double(hdfread(target_ist_file,'MOD_Swath_Sea_Ice/Data Fields/Ice_Surface_Temperature'))*0.01;
    % 这里是粗糙位置 5km
    % 粗糙经纬度的格网个数是 406*271/406*270
    ist_coarse_lon = double(hdfread(target_ist_file,'MOD_Swath_Sea_Ice/Geolocation Fields/Longitude'));
    ist_coarse_lat = double(hdfread(target_ist_file,'MOD_Swath_Sea_Ice/Geolocation Fields/Latitude'));

    % 读取ist对应的经纬度 是1km分辨率
    % 这里格网个数是对应mod_ist，2030*1354/2030*1350
    ist_lat = double(hdfread(target_coord_file,'MODIS_Swath_Type_GEO/Geolocation Fields/Latitude'));
    ist_lon = double(hdfread(target_coord_file,'MODIS_Swath_Type_GEO/Geolocation Fields/Longitude'));

    % 读取云掩膜云层覆盖数据
    % 这里格网个数同样对应mod_ist，2030*1354/2030*1350,但是注意ist4cloud是6*2030*1354/6*2030*1350的三维矩阵，
    ist4cloud = double(hdfread(target_cloud_file,'mod35/Data Fields/Cloud_Mask'));
    [~,ist_c] = size(mod_ist);
    if ist_c == 1354
        % 剔除不可识别的 多余数据
        mod_ist(:,end-3:end)=[];
        ist_lat(:,end-3:end)=[];
        ist_lon(:,end-3:end)=[];
        ist4cloud(:,:,end-3:end)=[];
        ist_coarse_lon(:,end)=[];
        ist_coarse_lat(:,end)=[];
    end

    [myd_x, myd_y] = projfwd(proj, ist_lat, ist_lon);  % 与myd03_lon/myd03_lat同尺寸

    ist_output = {mod_ist,ist_coarse_lon,ist_coarse_lat,ist_lat,ist_lon,ist4cloud,myd_x,myd_y,target_cloud_file};
end

function expanded_solar_zenith = solar_zenith_exact(target_cloud_file,mod_ist)
    % 注意，提取得出的太阳天顶角的分辨率是粗糙的，是5KM*5KM的
    % 并且注意，coarse_coord的形状不一定与ist4solar_zenith不一致，这就涉及到一副MODIS影像的大小，会出现2030*1350或者2030*1354的大小
    % 但是天顶角的数据都是406*270的size，而输入参数中都是经过剔除多余列数的数据
    % 注意 天顶角的定义是自己天顶方向和太阳的夹角，因此，若天顶角大于90度，就说明是晚上。
    ist4solar_zenith = double(hdfread(target_cloud_file,'mod35/Data Fields/Solar_Zenith'))*0.01;
    [out_x, out_y] = size(mod_ist);
    [rz, cz] = size(ist4solar_zenith);
    bx = floor(out_x / rz);
    by = floor(out_y / cz);

    expanded_solar_zenith = kron(ist4solar_zenith, ones(bx, by));
    expanded_solar_zenith = expanded_solar_zenith(1:out_x, 1:out_y);
end

function [surf_clear_ist,surf_non_cloud_ist] = ist_no_cloud_exact(mod_ist,ist4cloud,ist_lon,ist_lat,ist_r,ist_c,era_cloud_lon,era_cloud_lat, era_cloud_hcc,era_cloud_mcc,target_doy,target_hour,myd_x,myd_y) 
    %% 该函数是利用MODIS的云掩膜数据以及era5的云层数据 进行提取无云处的表面温度 
    % 第一个8bit 
    % 注意在bits中开始是从右向左数的，编号从0开始,但是在该提取的数数据中是第8位
    cloud1 = ist4cloud(1,:,:);
    cloud1_2d = reshape(cloud1,ist_r,ist_c);
    cloud1_8bit = dec2bin(cloud1_2d,8);
    %% 以下是只利用myd35云掩膜数据进行剔除的表面温度数据
    % 首先判断bit0(cloud1_8bit的第8位)，是否为1 即从右往左数第一个， 1的意义是确定下来的，如果不是1 那么删除该pixel
    % 接下判断bit1-2(cloud1_8bit的第6-7位),'00'代表云，'11'代表晴0，'01','10'代表不确定。晴的概率分别是0,0.33,0.67,1
    % bit3(cloud1_8bit的第5位) 代表是夜间还是白天 '0'代表night, '1'代表白天
    % bit6-7(cloud1_8bit的第1-2位) 代表是否是陆地还是水体  '00'water,'01'coastal,'10'desert,'11'Land
    % 即mask_clear_inx最后确定下来的就是 水体上方夜间晴空的确定数据
    mask_clear_inx = ismember(cloud1_8bit(:,8),'1') & ismember(cloud1_8bit(:,6:7),'11','rows')... 
                                                   & ismember(cloud1_8bit(:,5),'0') ...
                                                   & ismember(cloud1_8bit(:,1:2),'00','rows');
    mask_clear_2d = reshape(mask_clear_inx,ist_r,ist_c);
    surf_clear_ist = [ist_lon(mask_clear_2d),ist_lat(mask_clear_2d),mod_ist(mask_clear_2d)];
    
    lat_inx = surf_clear_ist(:,2)<66 | surf_clear_ist(:,2)>85;
    surf_clear_ist(lat_inx,:)=[];
    ist_inx= surf_clear_ist(:,3)<243;
    surf_clear_ist(ist_inx,:)=[];
    
    %% 以下是根据云掩膜以及ERA5再分析数据中的高中云层数据进行剔除
    % 即mask_cloud_inx最后确定下来的就是 只要有云的地方
    % 按照10*10的block处理云层，大于10的情况下，需要剔除
    % 该操作是依照 Makynen 2017文章中剔除云层以及反演海冰厚度的操作
    
    % 首先对10*10的block进行判断，若是>10% 则判定为云
    block_size = 10;
    mask_cloud_inx = ismember(cloud1_8bit(:,8),'1') & ismember(cloud1_8bit(:,6:7),'00','rows');
    mask_cloud_2d = double(reshape(mask_cloud_inx,ist_r,ist_c));
    mask_cloud_nan = mask_cloud_2d;
    
    nr = floor(ist_r/block_size)*block_size;
    nc = floor(ist_c/block_size)*block_size;

    M = mask_cloud_nan(1:nr, 1:nc);    % 裁成可整除大小
    % 计算每个 10x10 块的像元计数（和原逻辑一致）
    S = reshape(M, block_size, nr/block_size, block_size, nc/block_size);
    S = squeeze(sum(sum(S,1),3));      % (nr/block) x (nc/block)
    bad = (S > 10);                     % 原标准：>10 判为云

    % 把坏块“膨胀回”像元分辨率
    bad_rep = kron(bad, true(block_size));
    bad_rep = logical(bad_rep);
    bad_rep = bad_rep(1:nr, 1:nc);               % 保险：裁到与 tmp 同尺寸
    % 写回到 mask_cloud_nan
    tmp = mask_cloud_nan(1:nr,1:nc);
    tmp(bad_rep) = NaN;
    mask_cloud_nan(1:nr,1:nc) = tmp;

    % 余下边界（整除外的最后几列/行）按你的原逻辑保持为 NaN
    if nc < ist_c, mask_cloud_nan(:,nc+1:end) = NaN; end
    if nr < ist_r, mask_cloud_nan(nr+1:end,:) = NaN; end

    d2h = (target_doy-1)*24 + target_hour;

    % 再判断陆地的情况下，也设置为nan  该表达式为不是海水
    mask_land_inx = ~ismember(cloud1_8bit(:,1:2),'00','rows');
    mask_land_2d = reshape(mask_land_inx,ist_r,ist_c);
    % 将陆地的值设置为nan 
    mask_cloud_nan(mask_land_2d) = nan;
    % 再判断小于10的情况下的高云和中云的覆盖情况
    % 网格化ERA5
    target_cloud_hcc = era_cloud_hcc(:,:,d2h);   % [nLon x nLat]
    target_cloud_mcc = era_cloud_mcc(:,:,d2h);

    lon_vec = era_cloud_lon(:)';   % 1 x nLon
    lat_vec = era_cloud_lat(:)';   % 1 x nLat

    if any(lon_vec > 180)
        lon_vec = lon_vec - 360;
        [lon_vec, ix] = sort(lon_vec);
        target_cloud_hcc = target_cloud_hcc(ix, :);
        target_cloud_mcc = target_cloud_mcc(ix, :);
    end
    if numel(lat_vec) > 1 && lat_vec(2) < lat_vec(1)
        lat_vec = fliplr(lat_vec);
        target_cloud_hcc = target_cloud_hcc(:, end:-1:1);
        target_cloud_mcc = target_cloud_mcc(:, end:-1:1);
    end

    F_hcc = griddedInterpolant({lon_vec, lat_vec}, target_cloud_hcc, 'linear','none');
    F_mcc = griddedInterpolant({lon_vec, lat_vec}, target_cloud_mcc, 'linear','none');

    hcc_interp = F_hcc(ist_lon, ist_lat);
    mcc_interp = F_mcc(ist_lon, ist_lat);

    inx_cc = hcc_interp > 0.33 | mcc_interp > 0.51;
    mask_cloud_nan(inx_cc) = nan;

    % --- 剔除小晴空区块 ---
    min_clear_size = 3;
    min_pixels = min_clear_size^2;
    clear_mask = ~isnan(mask_cloud_nan);
    mask_cloud_nan(~bwareaopen(clear_mask, min_pixels, 8)) = nan;
    non_cloud_inx = ~isnan(mask_cloud_nan);
    surf_non_cloud_ist = [ist_lon(non_cloud_inx),ist_lat(non_cloud_inx),mod_ist(non_cloud_inx)];
    % 将符合纬度的挑选出来以及modis表面温度的有效值是243-273值
    indx = surf_non_cloud_ist(:,2)<66 | surf_non_cloud_ist(:,2)>85 | surf_non_cloud_ist(:,1)<-170 | surf_non_cloud_ist(:,1)>-120 |surf_non_cloud_ist(:,3)<243 | surf_non_cloud_ist(:,3)>273;
    surf_non_cloud_ist(indx,:)=[];
%     if isempty(surf_non_cloud_ist)
%         ...
%     end   
end

function mod2sit = ist2sit(surf_aux,F_Ta,F_U2,nswf)
    % 注意这里的surf_aux只有一组数据
    a = 5.67e-8;        % 玻尔兹曼常数
    p = 1.3;                % 空气密度
    C = 1.44e-3;            % 感热与潜热中性体传导系数
    Li = 2.86e6;
    epsF = 5;     % 阈值
    
    lon = surf_aux(1); lat = surf_aux(2); Ts = surf_aux(3);
    Ta  = F_Ta(lon,lat);         % K
    U2  = F_U2(lon,lat);         % m/s
    if isnan(Ta) || isnan(U2), mod2sit = NaN; return; end    
    
    lwd_u = 0.99*a*Ts.^4;
    lwd_d = 0.7855*a*Ta.^4;
    Fs = p*1004*C*U2*(Ta - Ts);
    Fe = 0.622*p*Li*C*U2*(0.8*SVP_kaleschke(Ta-273.15) - SVP_kaleschke(Ts-273.15))/1000;
    Fc = nswf + lwd_d - lwd_u + Fs + Fe;
    if abs(Fc) < epsF, mod2sit = NaN; return; end

    Ki = 2.034 + 0.13*8/((Ts+271.35)/2-273);
    mod2sit = Ki * (271.35 - Ts) / Fc;
    if mod2sit < 0 || mod2sit > 1, mod2sit = NaN; end
end

function sit_awi = awi_exact_sit(str_ib_date)
    awi_parent_path = 'W:\ljx_aux\SMOS_SIT_prod\AWI-SMOS\SMOS';
    
    dt  = datetime(str_ib_date,'InputFormat','yyyyMMdd');
    dtm = dt - days(1);  % 前一天
    target_year = datestr(dtm,'yyyy');
    target_mon  = datestr(dtm,'mm');
    target_awi_date = datestr(dtm,'yyyymmdd');
    
    target_awi_name = ['SMOS_Icethickness_v3.3_north_',target_awi_date,'.nc'];
    target_awi_file = fullfile(awi_parent_path,target_year,target_mon,target_awi_name);
    awi_lat = ncread(target_awi_file,'latitude');
    awi_lon = ncread(target_awi_file,'longitude');
    awi_sit = ncread(target_awi_file,'sea_ice_thickness');
    
    awi_inx = awi_lat<80 & awi_lat>66 & awi_lon > -170 & awi_lon < -120 & awi_sit>=0 & awi_sit < 1;
    sit_awi = [awi_lon(awi_inx),awi_lat(awi_inx),awi_sit(awi_inx)];
    
end

function P = SVP_kaleschke(td)
    P = 6.11 * 10.^(9.5*td./(265.5+td));
end
