function sit2thermodynamic_equilibrium(result_doy,ist_parent_path,ist_coord_parent_path,cloud_parent_path,sensor,night,target_year,era_path,era_mat_path,str_ib_date)
    era_output = era_data_exact(era_path,target_year,era_mat_path);
    [era_cloud_lat, era_cloud_lon, era_cloud_hcc, era_cloud_mcc, ...
     era_lon, era_lat, era_wind_u, era_wind_v, era_t2m] = era_output{:};
    [era_r,era_c] = size(era_lon);
    
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
    ist_hdf_list = dir([target_ist_path,'\*.','hdf']);
    % 温度1KM坐标数据
    ist_coord_list = dir([ist_coord_path,'\*.','hdf']);
    %modis云掩膜数据
    cloud_hdf_list = dir([cloud_path,'\*.','hdf']);
    
    % 根据得出的数据信息，进而对对应的IB数据日期进行提取
    % 每个数据天下，会有多个IB数据
    matches_ist = contains({ist_hdf_list.name}, result_doy);
    filtered_name = {ist_hdf_list(matches_ist).name}';
    ib_data_range = length(filtered_name);
    target_month = str2double(str_ib_date(5:6));
    for i = 1:length(ib_data_range)
        target_ist_name = filtered_name{i};
        timestamp = target_ist_name(8:19);
        matches_coord = contains({ist_coord_list.name}, timestamp);
        matches_cloud = contains({cloud_hdf_list.name}, timestamp);
        if sum(matches_coord)==0 || sum(matches_coord)==0
            continue
        else
            fprintf('success in %s',timestamp);
            target_coord_name = ist_coord_list(matches_coord).name;
            target_cloud_name = cloud_hdf_list(matches_cloud).name;
           
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
            % 转换为2m的风速
            wind_2m_speed = wind_10m_speed/1.27;

            mod_aux = sortrows(double([reshape(era_lon,[era_r*era_c 1]),reshape(era_lat,[era_r*era_c 1]),... 
                                    reshape(target_air_t2m,[era_r*era_c 1]), reshape(wind_2m_speed, [era_r*era_c 1])]));
            
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
                                                        
            sort_non_cloud_ist = roundn(sortrows(surf_non_cloud_ist),-4);                                                                                       
            % 先对夜间数据进行处理，得出夜间海冰厚度，作为计算日间海冰厚度的初始值,并再设定范围，进行迭代计算
           %% 先是夜间计算
%             nswf = 0;
%             sit_night  = ist2sit(surf_non_cloud_ist,mod_aux,nswf);
%             sit_inx = sit_night < 0 | sit_night > 1;
%             sit_night(sit_inx,:)=[];

            expanded_solar_zenith = solar_zenith_exact(target_cloud_file,mod_ist);
           %% 白天计算 得出应用的太阳天顶角的数据矩阵2030*1350
            inx_zenith = ist_lon < -120 & ist_lon > -170 & ist_lat>66 & ist_lat<80;
            sort_solar_zenith = roundn(sortrows([ist_lon(inx_zenith),ist_lat(inx_zenith),expanded_solar_zenith(inx_zenith)]),-4);
            
            ist2zenith_inx = ismember(sort_solar_zenith(:,1:2),sort_non_cloud_ist(:,1:2),'rows');
            % 得出有效冰面温度下的太阳天顶角
            target_solar_zenith = sort_solar_zenith(ist2zenith_inx,:);
            % 得出和有效太阳天顶角下的冰面温度
            zenith2ist_inx = ismember(sort_non_cloud_ist(:,1:2),target_solar_zenith(:,1:2),'rows');
            target_non_cloud_ist = sort_non_cloud_ist(zenith2ist_inx,:);
            % target_solar_zenith；target_non_cloud_ist 二者经纬度和序列都一致
            % 输出文件为target_non_cloud_sit,第三列拟合短波辐射海冰厚度；第四列为参数化短波辐射海冰厚度
            target_non_cloud_sit = target_non_cloud_ist;
            % 提取白天的前一天的 awi 海冰厚度,作为初始海冰厚度值,
            % 输出AWI海冰厚度n*3数据矩阵
            sit_awi = awi_exact_sit(str_ib_date);
            % 将sit_awi海冰厚度对应到每一个有效冰面温度，由于只是初始值，所以利用最邻近法进行赋值
            wgs84 = referenceEllipsoid('WGS84');
            for j = size(target_non_cloud_ist,1)
                j_coord = target_non_cloud_ist(j,1:2);
                j_angdist = distance(j_coord(2),j_coord(1),sit_awi(:,2),sit_awi(:,1),wgs84);
                j_dist = deg2km(j_angdist);
                j_inx = j_dist == min(j_dist);
                reference_sit = sit_awi(j_inx,3);
                
                nswf_fit = net_shortwave_flux_fit(reference_sit,target_month);
                target_non_cloud_sit(j,3) = ist2sit(target_non_cloud_ist(j,:),mod_aux,nswf_fit);
                 
                sit_max = max(0.9,1.5*reference_sit);
                sit4nswf_range = [0.05:0.025:0.5,0.55:0.05:sit_max];
                sit_candidate4nswf_parameterization = sit4nswf_range;
                for k = 1:length(sit4nswf_range)
                    nswf_parameterization = net_shortwave_flux_parameterization(target_solar_zenith(k),sit4nswf_range(k));
                    sit_candidate4nswf_parameterization(k) = ist2sit(target_non_cloud_ist(j,:),mod_aux,nswf_parameterization);
                end
                sit_diff = abs(sit4nswf_range-sit_candidate4nswf_parameterization);
                select_inx = sit_diff == min(sit_diff);
                target_non_cloud_sit(j,4) = sit_candidate4nswf_parameterization(select_inx);
            end
        end
    end
    output_path = 'W:\ljx_aux\MODIS\output_sit';
    output_name = [timestamp,'.mat'];
    save(fullfile(output_path,output_name),'target_non_cloud_sit')
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
    [output_size_x, output_size_y] = size(mod_ist);
    expanded_solar_zenith = NaN(output_size_x, output_size_y);
    for i = 1:406
        for j = 1:270
            % 计算该天顶角值在输出矩阵中的对应位置
            row_start = (i - 1) * 5 + 1;  % 行的起始索引
            row_end = i * 5;  % 行的结束索引
            col_start = (j - 1) * 5 + 1;  % 列的起始索引
            col_end = j * 5;  % 列的结束索引

            % 将天顶角数据赋值到对应的5x5KM区域 
            expanded_solar_zenith(row_start:row_end, col_start:col_end) = ist4solar_zenith(i, j);
        end
    end
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
    num_rows = floor(ist_r / block_size);  % 划分的行数
    num_cols = floor(ist_c / block_size);  % 划分的列数

    d2h = (target_doy-1)*24 + target_hour;

    grid_count = zeros(num_rows, num_cols); % 用于记录每个格网块内的云层覆盖像元数
    for i = 1:num_rows
        for j = 1:num_cols
            % 当前格网块的范围
            row_min = (i-1) * block_size + 1;
            row_max = i * block_size;
            col_min = (j-1) * block_size + 1;
            col_max = j * block_size;

            % 查找当前格网块内的云层覆盖像元
            in_grid = mask_cloud_nan(row_min:row_max, col_min:col_max);
            cloud_block_coverage = sum(in_grid(:));
            grid_count(i, j) = cloud_block_coverage;
            % cloud_block_coverage >10 就说明是云的
            if cloud_block_coverage > 10
                mask_cloud_nan(row_min:row_max, col_min:col_max) = nan;
            end
        end
    end   
    % 再判断陆地的情况下，也设置为nan  该表达式为不是海水
    mask_land_inx = ~ismember(cloud1_8bit(:,1:2),'00','rows');
    mask_land_2d = reshape(mask_land_inx,ist_r,ist_c);
    % 将陆地的值设置为nan 
    mask_cloud_nan(mask_land_2d) = nan;
    % 再判断小于10的情况下的高云和中云的覆盖情况
    % 网格化ERA5
    [era_lon_grid, era_lat_grid] = meshgrid(era_cloud_lon, era_cloud_lat); 
    % 调整云数据维度（1440×101×time → 101×1440×time）
    target_cloud_hcc = permute(era_cloud_hcc(:,:,d2h), [2,1]);
    target_cloud_mcc = permute(era_cloud_mcc(:,:,d2h), [2,1]);

    proj = projcrs(3413);  % WGS84北极极射投影
    [era_x, era_y] = projfwd(proj, era_lat_grid, era_lon_grid);  % 输入101×1440网格

    % 插值到MODIS像元
    F_hcc = scatteredInterpolant(era_x(:), era_y(:), target_cloud_hcc(:), 'linear', 'none');
    hcc_interp = F_hcc(myd_x, myd_y);  % myd_x和myd_y已提前投影
    M_hcc = scatteredInterpolant(era_x(:), era_y(:), target_cloud_mcc(:), 'linear', 'none');
    mcc_interp = M_hcc(myd_x, myd_y);  % myd_x和myd_y已提前投影
    inx_cc = hcc_interp > 0.33 | mcc_interp > 0.51;
    mask_cloud_nan(inx_cc) = nan;
    % 将后几列没有处理的数据赋值nan
    mask_cloud_nan(:,num_cols*10+1:end) = nan;
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

function mod2sit = ist2sit(surf_aux,mod_aux,nswf)
    % 注意这里的surf_aux只有一组数据
    a = 5.67e-8;        % 玻尔兹曼常数
    p = 1.3;                % 空气密度
    C = 1.44e-3;            % 感热与潜热中性体传导系数
    Li = 2.86e6;
    
    % 将剔除后的表面温度转化为海冰厚度
    R = 6371; % 地球半径（km）
    % 转换为三维坐标
    x = R .* cosd(mod_aux(:,2)) .* cosd(mod_aux(:,1));
    y = R .* cosd(mod_aux(:,2)) .* sind(mod_aux(:,1));
    z = R .* sind(mod_aux(:,2));
    era_coords = [x, y, z];
    tree = KDTreeSearcher(era_coords); % 构建KDTree

    surf_lon = surf_aux(1);
    surf_lat = surf_aux(2);
    x_surf = R .* cosd(surf_lat) .* cosd(surf_lon);
    y_surf = R .* cosd(surf_lat) .* sind(surf_lon);
    z_surf = R .* sind(surf_lat);
    surf_coords = [x_surf, y_surf, z_surf];
    [idx, ~] = knnsearch(tree, surf_coords, 'K', 1); % 批量查询最近邻
    surf_aux(4:5) = mod_aux(idx, 3:4); % 提取温度和风速
    inx = surf_aux(4)>273.15-5;
    surf_aux(inx,:)=[];

    lwd_u = 0.99*a*surf_aux(3).^4;   % 长波 上行 
%         e = 0.7855*(1+0.2232*0.4^0.75);
    lwd_d = 0.7855*a*surf_aux(4).^4; % 长波下行     
    % 感热通量
    Fs = p*1004*C*surf_aux(5)*(surf_aux(4) - surf_aux(3));
    % 潜热通量 Fe
%         Fe_MXY = 0.622*p*Li*C*sort_surf_T(:,5).*SVP_MXY(1000,sort_surf_T(:,6)-273.15,sort_surf_T(:,3)-273.15)./1000; 
    % 另一种 潜热通量计算方式
%         Fe = 0.622*p*2.49e3*C*S_Loca(:,5).*(0.9*SVP(S_Loca(:,4))-SVP(S_Loca(:,3)))./1000;
    % 第三种计算潜热通量计算方式
    Fe = 0.622*p*Li*C*surf_aux(5)*(0.8*SVP_kaleschke(surf_aux(4)-273.15)-SVP_kaleschke(surf_aux(3)-273.15))/1000;

    % lwd_d - lwd_u + Fe + Fs + Fc = 0;
    Fc = nswf + lwd_d - lwd_u + Fs + Fe;
    Ki = 2.034+0.13*8/((surf_aux(3)+271.35)/2-273);
    mod2sit = Ki * (271.35 - surf_aux(3))/Fc;
%     surf_aux(:,6) = mod2sit;
    if mod2sit <0 || mod2sit>1
        mod2sit = nan;
    end
%     mod2sit_inx = mod2sit<0 | mod2sit>1;
%     surf_aux(sit_inx,:)=[];
%     surf_aux(3:5)=[];
end

function sit_awi = awi_exact_sit(str_ib_date)
    awi_parent_path = 'W:\ljx_aux\SMOS_SIT_prod\AWI-SMOS\SMOS';
    target_year = str_ib_date(1:4);
    target_mon = str_ib_date(5:6);
    target_awi_date = num2str(str2double(str_ib_date)-1);
    target_awi_name = ['SMOS_Icethickness_v3.3_north_',target_awi_date,'.nc'];
    target_awi_file = fullfile(awi_parent_path,target_year,target_mon,target_awi_name);
    awi_lat = ncread(target_awi_file,'latitude');
    awi_lon = ncread(target_awi_file,'longitude');
    awi_sit = ncread(target_awi_file,'sea_ice_thickness');
    
    awi_inx = awi_lat<80 & awi_lat>66 & awi_lon <-170 & awi_lon>-120 & awi_sit>=0 & awi_sit < 1;
    sit_awi = [awi_lon(awi_inx),awi_lat(awi_inx),awi_sit(awi_inx)];
    
end

function P = SVP_kaleschke(td)
    P = 6.11 * 10.^(9.5*td./(265.5+td));
end
