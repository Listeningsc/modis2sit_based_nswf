function sit2thermodynamic_equilibrium(result_doy,ist_parent_path,ist_coord_parent_path,cloud_parent_path,sensor,night,target_year,era_path,era_mat_path,str_ib_date)
    era_output = era_data_exact(era_path,target_year,era_mat_path);
    [era_cloud_lat, era_cloud_lon, era_cloud_hcc, era_cloud_mcc, ...
     era_lon, era_lat, era_wind_u, era_wind_v, era_t2m] = era_output{:};
 
    %% modis ����·��
    % �����¶�·��
    % day_path ��ʾ�ǰ��컹������
    if night
        day_path = 'night';
    else
        day_path = 'day';
    end
    target_ist_path = fullfile(ist_parent_path,sensor,day_path,num2str(target_year));
    % ����·��
    ist_coord_path = fullfile(ist_coord_parent_path,sensor,day_path,num2str(target_year));
    % ����Ĥ·��
    cloud_path = fullfile(cloud_parent_path,sensor,day_path,num2str(target_year)); 
    % �¶�����
    ist_hdf_list = dir(fullfile(target_ist_path,'*.hdf'));
    % �¶�1KM��������
    ist_coord_list = dir(fullfile(ist_coord_path,'*.hdf'));
    %modis����Ĥ����
    cloud_hdf_list = dir(fullfile(cloud_path,'*.hdf'));
    
    % ���ݵó���������Ϣ�������Զ�Ӧ��IB�������ڽ�����ȡ
    % ÿ���������£����ж��IB����
    matches_ist = contains({ist_hdf_list.name}, result_doy);
    filtered_name = {ist_hdf_list(matches_ist).name}';
    target_month = str2double(str_ib_date(5:6));
    for i = 1:numel(filtered_name)
        target_ist_name = filtered_name{i};
        timestamp = target_ist_name(8:19);
        
        %era5�������� �ֽ��¶ȵ�doy�����ж��Ƿ���Ҫ���м��·��޳������ҽ�ѡȡ��ʱ��㶼��ȡ����
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
        % ===== ERA 2m ����/���ٲ�ֵ�������� 1440x84��lon�ڵ�1ά��lat�ڵ�2ά�� =====
        Ta2D = target_air_t2m;      % ���� size: [1440 x 84]  (ά1=lon, ά2=lat)
        U22D = wind_2m_speed;       % ͬ��
       % ��������Ϊ [84 x 1440]���ɼ�һ�а�ȫת�ã�
        if size(Ta2D,1)==84 && size(Ta2D,2)==1440
            Ta2D = Ta2D.';  U22D = U22D.';
        end
        
        % ��ȷ���᣺lon ȡ��1�У����б仯����lat ȡ��1�У����б仯��
        lon_vec = double(era_lon(:,1));   % 1440x1
        lat_vec = double(era_lat(1,:));   % 1x84
        
       % --- ���ȴ���ӳ�䵽 [-180,180)������ͬά���������ݣ���1ά�� ---
        lon_vec(lon_vec >= 180) = lon_vec(lon_vec >= 180) - 360;
        [lon_vec, lon_ord] = sort(lon_vec, 'ascend');
        Ta2D = Ta2D(lon_ord, :);
        U22D = U22D(lon_ord, :);
        % --- γ�ȴ�����Ϊ������ת������ͬά���������ݣ���2ά�� ---
        if numel(lat_vec) > 1 && lat_vec(2) < lat_vec(1)
            lat_vec = fliplr(lat_vec);
            Ta2D = Ta2D(:, end:-1:1);
            U22D = U22D(:, end:-1:1);
        end
           % --- ���գ�ȥ�أ����ټ������ɱ���0/360�ظ�֮������� ---
        [lon_vec, ia_lon] = unique(lon_vec, 'stable');  Ta2D = Ta2D(ia_lon, :);  U22D = U22D(ia_lon, :);
        [lat_vec, ia_lat] = unique(lat_vec, 'stable');  Ta2D = Ta2D(:, ia_lat);  U22D = U22D(:, ia_lat);

        % ���� griddedInterpolant��V �ĳߴ������ [numel(lon_vec) x numel(lat_vec)]��
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

            % ע�⣬�������������������޳��� ��������������
            % ע��ist_coarse_lon �� ist_lon ����-180-180�ȵķ�Χ
            ist_output = ist_data_exact(target_ist_path,ist_coord_path,cloud_path,target_ist_name,target_coord_name,target_cloud_name);
            [mod_ist,~,~,ist_lat,ist_lon,ist4cloud,myd_x,myd_y,target_cloud_file] = ist_output{:};
            [ist_r,ist_c]= size(mod_ist);
             % ��ȡ���������µĺ��������¶�
             % �ó��������������һ��Ϊ����Modis����Ĥ�ó��������µı������¶ȣ��ڶ���Ϊ����block����ʶ���ó������Ʊ����¶�
             % surf_non_cloud_ist��n*3�����ݾ���
            [~,surf_non_cloud_ist] = ist_no_cloud_exact(mod_ist,ist4cloud,ist_lon,ist_lat,ist_r,ist_c,...
                                                        era_cloud_lon,era_cloud_lat, era_cloud_hcc,era_cloud_mcc,...
                                                            target_doy,target_hour,myd_x,myd_y);           
                                                        
           %% FIX: �ռ�����
            if isempty(surf_non_cloud_ist)
                warning('No clear-sky IST for %s, skip.', timestamp);
                continue
            end
            
            % �ȶ�ҹ�����ݽ��д����ó�ҹ�亣����ȣ���Ϊ�����ռ亣����ȵĳ�ʼֵ,�����趨��Χ�����е�������
            expanded_solar_zenith = solar_zenith_exact(target_cloud_file,mod_ist);
           %% ������� �ó�Ӧ�õ�̫���춥�ǵ����ݾ���2030*1350
            inx_zenith = ist_lon < -120 & ist_lon > -170 & ist_lat>66 & ist_lat<80;
                        % ����������
            nc = sortrows([round(surf_non_cloud_ist(:,1),4), ...
                           round(surf_non_cloud_ist(:,2),4), ...
                           surf_non_cloud_ist(:,3)], [1 2]);           % N��3

            sz = sortrows([round(ist_lon(inx_zenith),4), ...
                           round(ist_lat(inx_zenith),4), ...
                           expanded_solar_zenith(inx_zenith)], [1 2]); % M��3
            % �� ismember �Ķ�λ������֤һһ��Ӧ����ͬ˳��
            [tf, loc] = ismember(nc(:,1:2), sz(:,1:2), 'rows');
            target_non_cloud_ist = nc(tf,:);         % N'��3
            target_solar_zenith  = sz(loc(tf),:);    % N'��3������ͬ��

            %% FIX: �ٴοռ�����
            if isempty(target_non_cloud_ist)
                warning('No matched SZA-IST pixels for %s, skip.', timestamp);
                continue
            end     
            % ---- AWI ǰһ�� SIT ��Ϊ��ֵ ----
            sit_awi = awi_exact_sit(str_ib_date);

            %% FIX: �������гأ�ֻ��һ�Σ�����ÿ������ʱ�رգ�
            if isempty(gcp('nocreate'))
                parpool(4);   % �� MATLAB ��ѡ���жȣ���Ҳ����̶������� parpool(4)
            end
            npts = size(target_non_cloud_ist,1);
            sit_fit   = nan(npts,1);   % ��3��
            sit_param = nan(npts,1);   % ��4��

            % �㲥����������ÿ�� worker ����/���죩
            const_F_Ta = parallel.pool.Constant(F_Ta);
            const_F_U2 = parallel.pool.Constant(F_U2);
            const_AWI  = parallel.pool.Constant(sit_awi);               
            constEll   = parallel.pool.Constant(wgs84Ellipsoid('km')); 
            % ---- ����ÿ����Ԫ ----
            parfor j = 1:npts
                lon_j = target_non_cloud_ist(j,1);
                lat_j = target_non_cloud_ist(j,2);
                Ts_j  = target_non_cloud_ist(j,3);   % K
                sza_j = target_solar_zenith(j,3);    % deg

                % ��� AWI ��ֵ
                awi  = const_AWI.Value;    % [lon, lat, sit]
                ang  = distance(lat_j, lon_j, awi(:,2), awi(:,1), constEll.Value);
                [~, j_inx]     = min(deg2km(ang));
                reference_sit  = awi(j_inx,3);

                % ���ߣ�ʹ�������ϣ������ò����������滻��һ�У�
                nswf_fit = net_shortwave_flux_fit(reference_sit, target_month);
                tmp3     = ist2sit([lon_j, lat_j, Ts_j], const_F_Ta.Value, const_F_U2.Value, nswf_fit);
                if isnan(tmp3), tmp3 = -1; end

                % ɨ�������
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
                    [~,best_idx] = min(abs(svec - cvec));   %% FIX: ��һ argmin�������Ԫ�ظ�ֵ
                    tmp4 = cvec(best_idx);
                else
                    tmp4 = NaN;
                end

                sit_fit(j)   = tmp3;   % ��Ƭд��
                sit_param(j) = tmp4;
            end

            target_non_cloud_sit = [target_non_cloud_ist, sit_fit, sit_param];
            % ---- ���-----
            output_path = 'W:\ljx_aux\MODIS\output_sit';
             if ~exist(output_path,'dir'); mkdir(output_path); end
            output_name = [timestamp,'.mat'];
            save(fullfile(output_path,output_name),'target_non_cloud_sit');
        end
    end
end

function era_output = era_data_exact(era_path,target_year,era_mat_path)
                                                                                                        
 %%  era5 ��������
    % era5 ��ͬ�߶�����Ϣ
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
    proj = projcrs(3413);  % WGS84��������ͶӰ
    target_ist_file = fullfile(target_ist_path,target_ist_name);
    target_coord_file = fullfile(ist_coord_path,target_coord_name);
    target_cloud_file = fullfile(cloud_path,target_cloud_name);

    % ��ȡist�����¶� 1km�ֱ���, 
    % ע��mod_ist�ĸ���������2030*1354/2030*1350��������modis�Ĳ�ͬλ�þ����ģ��ٺ����Ҫ��1354���н���ɾ���������
    mod_ist = double(hdfread(target_ist_file,'MOD_Swath_Sea_Ice/Data Fields/Ice_Surface_Temperature'))*0.01;
    % �����Ǵֲ�λ�� 5km
    % �ֲھ�γ�ȵĸ��������� 406*271/406*270
    ist_coarse_lon = double(hdfread(target_ist_file,'MOD_Swath_Sea_Ice/Geolocation Fields/Longitude'));
    ist_coarse_lat = double(hdfread(target_ist_file,'MOD_Swath_Sea_Ice/Geolocation Fields/Latitude'));

    % ��ȡist��Ӧ�ľ�γ�� ��1km�ֱ���
    % ������������Ƕ�Ӧmod_ist��2030*1354/2030*1350
    ist_lat = double(hdfread(target_coord_file,'MODIS_Swath_Type_GEO/Geolocation Fields/Latitude'));
    ist_lon = double(hdfread(target_coord_file,'MODIS_Swath_Type_GEO/Geolocation Fields/Longitude'));

    % ��ȡ����Ĥ�Ʋ㸲������
    % �����������ͬ����Ӧmod_ist��2030*1354/2030*1350,����ע��ist4cloud��6*2030*1354/6*2030*1350����ά����
    ist4cloud = double(hdfread(target_cloud_file,'mod35/Data Fields/Cloud_Mask'));
    [~,ist_c] = size(mod_ist);
    if ist_c == 1354
        % �޳�����ʶ��� ��������
        mod_ist(:,end-3:end)=[];
        ist_lat(:,end-3:end)=[];
        ist_lon(:,end-3:end)=[];
        ist4cloud(:,:,end-3:end)=[];
        ist_coarse_lon(:,end)=[];
        ist_coarse_lat(:,end)=[];
    end

    [myd_x, myd_y] = projfwd(proj, ist_lat, ist_lon);  % ��myd03_lon/myd03_latͬ�ߴ�

    ist_output = {mod_ist,ist_coarse_lon,ist_coarse_lat,ist_lat,ist_lon,ist4cloud,myd_x,myd_y,target_cloud_file};
end

function expanded_solar_zenith = solar_zenith_exact(target_cloud_file,mod_ist)
    % ע�⣬��ȡ�ó���̫���춥�ǵķֱ����Ǵֲڵģ���5KM*5KM��
    % ����ע�⣬coarse_coord����״��һ����ist4solar_zenith��һ�£�����漰��һ��MODISӰ��Ĵ�С�������2030*1350����2030*1354�Ĵ�С
    % �����춥�ǵ����ݶ���406*270��size������������ж��Ǿ����޳���������������
    % ע�� �춥�ǵĶ������Լ��춥�����̫���ļнǣ���ˣ����춥�Ǵ���90�ȣ���˵�������ϡ�
    ist4solar_zenith = double(hdfread(target_cloud_file,'mod35/Data Fields/Solar_Zenith'))*0.01;
    [out_x, out_y] = size(mod_ist);
    [rz, cz] = size(ist4solar_zenith);
    bx = floor(out_x / rz);
    by = floor(out_y / cz);

    expanded_solar_zenith = kron(ist4solar_zenith, ones(bx, by));
    expanded_solar_zenith = expanded_solar_zenith(1:out_x, 1:out_y);
end

function [surf_clear_ist,surf_non_cloud_ist] = ist_no_cloud_exact(mod_ist,ist4cloud,ist_lon,ist_lat,ist_r,ist_c,era_cloud_lon,era_cloud_lat, era_cloud_hcc,era_cloud_mcc,target_doy,target_hour,myd_x,myd_y) 
    %% �ú���������MODIS������Ĥ�����Լ�era5���Ʋ����� ������ȡ���ƴ��ı����¶� 
    % ��һ��8bit 
    % ע����bits�п�ʼ�Ǵ����������ģ���Ŵ�0��ʼ,�����ڸ���ȡ�����������ǵ�8λ
    cloud1 = ist4cloud(1,:,:);
    cloud1_2d = reshape(cloud1,ist_r,ist_c);
    cloud1_8bit = dec2bin(cloud1_2d,8);
    %% ������ֻ����myd35����Ĥ���ݽ����޳��ı����¶�����
    % �����ж�bit0(cloud1_8bit�ĵ�8λ)���Ƿ�Ϊ1 ��������������һ���� 1��������ȷ�������ģ��������1 ��ôɾ����pixel
    % �����ж�bit1-2(cloud1_8bit�ĵ�6-7λ),'00'�����ƣ�'11'������0��'01','10'����ȷ������ĸ��ʷֱ���0,0.33,0.67,1
    % bit3(cloud1_8bit�ĵ�5λ) ������ҹ�仹�ǰ��� '0'����night, '1'�������
    % bit6-7(cloud1_8bit�ĵ�1-2λ) �����Ƿ���½�ػ���ˮ��  '00'water,'01'coastal,'10'desert,'11'Land
    % ��mask_clear_inx���ȷ�������ľ��� ˮ���Ϸ�ҹ����յ�ȷ������
    mask_clear_inx = ismember(cloud1_8bit(:,8),'1') & ismember(cloud1_8bit(:,6:7),'11','rows')... 
                                                   & ismember(cloud1_8bit(:,5),'0') ...
                                                   & ismember(cloud1_8bit(:,1:2),'00','rows');
    mask_clear_2d = reshape(mask_clear_inx,ist_r,ist_c);
    surf_clear_ist = [ist_lon(mask_clear_2d),ist_lat(mask_clear_2d),mod_ist(mask_clear_2d)];
    
    lat_inx = surf_clear_ist(:,2)<66 | surf_clear_ist(:,2)>85;
    surf_clear_ist(lat_inx,:)=[];
    ist_inx= surf_clear_ist(:,3)<243;
    surf_clear_ist(ist_inx,:)=[];
    
    %% �����Ǹ�������Ĥ�Լ�ERA5�ٷ��������еĸ����Ʋ����ݽ����޳�
    % ��mask_cloud_inx���ȷ�������ľ��� ֻҪ���Ƶĵط�
    % ����10*10��block�����Ʋ㣬����10������£���Ҫ�޳�
    % �ò��������� Makynen 2017�������޳��Ʋ��Լ����ݺ�����ȵĲ���
    
    % ���ȶ�10*10��block�����жϣ�����>10% ���ж�Ϊ��
    block_size = 10;
    mask_cloud_inx = ismember(cloud1_8bit(:,8),'1') & ismember(cloud1_8bit(:,6:7),'00','rows');
    mask_cloud_2d = double(reshape(mask_cloud_inx,ist_r,ist_c));
    mask_cloud_nan = mask_cloud_2d;
    
    nr = floor(ist_r/block_size)*block_size;
    nc = floor(ist_c/block_size)*block_size;

    M = mask_cloud_nan(1:nr, 1:nc);    % �óɿ�������С
    % ����ÿ�� 10x10 �����Ԫ��������ԭ�߼�һ�£�
    S = reshape(M, block_size, nr/block_size, block_size, nc/block_size);
    S = squeeze(sum(sum(S,1),3));      % (nr/block) x (nc/block)
    bad = (S > 10);                     % ԭ��׼��>10 ��Ϊ��

    % �ѻ��顰���ͻء���Ԫ�ֱ���
    bad_rep = kron(bad, true(block_size));
    bad_rep = logical(bad_rep);
    bad_rep = bad_rep(1:nr, 1:nc);               % ���գ��õ��� tmp ͬ�ߴ�
    % д�ص� mask_cloud_nan
    tmp = mask_cloud_nan(1:nr,1:nc);
    tmp(bad_rep) = NaN;
    mask_cloud_nan(1:nr,1:nc) = tmp;

    % ���±߽磨������������/�У������ԭ�߼�����Ϊ NaN
    if nc < ist_c, mask_cloud_nan(:,nc+1:end) = NaN; end
    if nr < ist_r, mask_cloud_nan(nr+1:end,:) = NaN; end

    d2h = (target_doy-1)*24 + target_hour;

    % ���ж�½�ص�����£�Ҳ����Ϊnan  �ñ��ʽΪ���Ǻ�ˮ
    mask_land_inx = ~ismember(cloud1_8bit(:,1:2),'00','rows');
    mask_land_2d = reshape(mask_land_inx,ist_r,ist_c);
    % ��½�ص�ֵ����Ϊnan 
    mask_cloud_nan(mask_land_2d) = nan;
    % ���ж�С��10������µĸ��ƺ����Ƶĸ������
    % ����ERA5
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

    % --- �޳�С������� ---
    min_clear_size = 3;
    min_pixels = min_clear_size^2;
    clear_mask = ~isnan(mask_cloud_nan);
    mask_cloud_nan(~bwareaopen(clear_mask, min_pixels, 8)) = nan;
    non_cloud_inx = ~isnan(mask_cloud_nan);
    surf_non_cloud_ist = [ist_lon(non_cloud_inx),ist_lat(non_cloud_inx),mod_ist(non_cloud_inx)];
    % ������γ�ȵ���ѡ�����Լ�modis�����¶ȵ���Чֵ��243-273ֵ
    indx = surf_non_cloud_ist(:,2)<66 | surf_non_cloud_ist(:,2)>85 | surf_non_cloud_ist(:,1)<-170 | surf_non_cloud_ist(:,1)>-120 |surf_non_cloud_ist(:,3)<243 | surf_non_cloud_ist(:,3)>273;
    surf_non_cloud_ist(indx,:)=[];
%     if isempty(surf_non_cloud_ist)
%         ...
%     end   
end

function mod2sit = ist2sit(surf_aux,F_Ta,F_U2,nswf)
    % ע�������surf_auxֻ��һ������
    a = 5.67e-8;        % ������������
    p = 1.3;                % �����ܶ�
    C = 1.44e-3;            % ������Ǳ�������崫��ϵ��
    Li = 2.86e6;
    epsF = 5;     % ��ֵ
    
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
    dtm = dt - days(1);  % ǰһ��
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
