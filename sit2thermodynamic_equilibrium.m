function sit2thermodynamic_equilibrium(result_doy,ist_parent_path,ist_coord_parent_path,cloud_parent_path,sensor,night,target_year,era_path,era_mat_path,str_ib_date)
    era_output = era_data_exact(era_path,target_year,era_mat_path);
    [era_cloud_lat, era_cloud_lon, era_cloud_hcc, era_cloud_mcc, ...
     era_lon, era_lat, era_wind_u, era_wind_v, era_t2m] = era_output{:};
    [era_r,era_c] = size(era_lon);
    
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
    ist_hdf_list = dir([target_ist_path,'\*.','hdf']);
    % �¶�1KM��������
    ist_coord_list = dir([ist_coord_path,'\*.','hdf']);
    %modis����Ĥ����
    cloud_hdf_list = dir([cloud_path,'\*.','hdf']);
    
    % ���ݵó���������Ϣ�������Զ�Ӧ��IB�������ڽ�����ȡ
    % ÿ���������£����ж��IB����
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
            % ת��Ϊ2m�ķ���
            wind_2m_speed = wind_10m_speed/1.27;

            mod_aux = sortrows(double([reshape(era_lon,[era_r*era_c 1]),reshape(era_lat,[era_r*era_c 1]),... 
                                    reshape(target_air_t2m,[era_r*era_c 1]), reshape(wind_2m_speed, [era_r*era_c 1])]));
            
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
                                                        
            sort_non_cloud_ist = roundn(sortrows(surf_non_cloud_ist),-4);                                                                                       
            % �ȶ�ҹ�����ݽ��д����ó�ҹ�亣����ȣ���Ϊ�����ռ亣����ȵĳ�ʼֵ,�����趨��Χ�����е�������
           %% ����ҹ�����
%             nswf = 0;
%             sit_night  = ist2sit(surf_non_cloud_ist,mod_aux,nswf);
%             sit_inx = sit_night < 0 | sit_night > 1;
%             sit_night(sit_inx,:)=[];

            expanded_solar_zenith = solar_zenith_exact(target_cloud_file,mod_ist);
           %% ������� �ó�Ӧ�õ�̫���춥�ǵ����ݾ���2030*1350
            inx_zenith = ist_lon < -120 & ist_lon > -170 & ist_lat>66 & ist_lat<80;
            sort_solar_zenith = roundn(sortrows([ist_lon(inx_zenith),ist_lat(inx_zenith),expanded_solar_zenith(inx_zenith)]),-4);
            
            ist2zenith_inx = ismember(sort_solar_zenith(:,1:2),sort_non_cloud_ist(:,1:2),'rows');
            % �ó���Ч�����¶��µ�̫���춥��
            target_solar_zenith = sort_solar_zenith(ist2zenith_inx,:);
            % �ó�����Ч̫���춥���µı����¶�
            zenith2ist_inx = ismember(sort_non_cloud_ist(:,1:2),target_solar_zenith(:,1:2),'rows');
            target_non_cloud_ist = sort_non_cloud_ist(zenith2ist_inx,:);
            % target_solar_zenith��target_non_cloud_ist ���߾�γ�Ⱥ����ж�һ��
            % ����ļ�Ϊtarget_non_cloud_sit,��������϶̲����亣����ȣ�������Ϊ�������̲����亣�����
            target_non_cloud_sit = target_non_cloud_ist;
            % ��ȡ�����ǰһ��� awi �������,��Ϊ��ʼ�������ֵ,
            % ���AWI�������n*3���ݾ���
            sit_awi = awi_exact_sit(str_ib_date);
            % ��sit_awi������ȶ�Ӧ��ÿһ����Ч�����¶ȣ�����ֻ�ǳ�ʼֵ�������������ڽ������и�ֵ
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
    [output_size_x, output_size_y] = size(mod_ist);
    expanded_solar_zenith = NaN(output_size_x, output_size_y);
    for i = 1:406
        for j = 1:270
            % ������춥��ֵ����������еĶ�Ӧλ��
            row_start = (i - 1) * 5 + 1;  % �е���ʼ����
            row_end = i * 5;  % �еĽ�������
            col_start = (j - 1) * 5 + 1;  % �е���ʼ����
            col_end = j * 5;  % �еĽ�������

            % ���춥�����ݸ�ֵ����Ӧ��5x5KM���� 
            expanded_solar_zenith(row_start:row_end, col_start:col_end) = ist4solar_zenith(i, j);
        end
    end
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
    num_rows = floor(ist_r / block_size);  % ���ֵ�����
    num_cols = floor(ist_c / block_size);  % ���ֵ�����

    d2h = (target_doy-1)*24 + target_hour;

    grid_count = zeros(num_rows, num_cols); % ���ڼ�¼ÿ���������ڵ��Ʋ㸲����Ԫ��
    for i = 1:num_rows
        for j = 1:num_cols
            % ��ǰ������ķ�Χ
            row_min = (i-1) * block_size + 1;
            row_max = i * block_size;
            col_min = (j-1) * block_size + 1;
            col_max = j * block_size;

            % ���ҵ�ǰ�������ڵ��Ʋ㸲����Ԫ
            in_grid = mask_cloud_nan(row_min:row_max, col_min:col_max);
            cloud_block_coverage = sum(in_grid(:));
            grid_count(i, j) = cloud_block_coverage;
            % cloud_block_coverage >10 ��˵�����Ƶ�
            if cloud_block_coverage > 10
                mask_cloud_nan(row_min:row_max, col_min:col_max) = nan;
            end
        end
    end   
    % ���ж�½�ص�����£�Ҳ����Ϊnan  �ñ��ʽΪ���Ǻ�ˮ
    mask_land_inx = ~ismember(cloud1_8bit(:,1:2),'00','rows');
    mask_land_2d = reshape(mask_land_inx,ist_r,ist_c);
    % ��½�ص�ֵ����Ϊnan 
    mask_cloud_nan(mask_land_2d) = nan;
    % ���ж�С��10������µĸ��ƺ����Ƶĸ������
    % ����ERA5
    [era_lon_grid, era_lat_grid] = meshgrid(era_cloud_lon, era_cloud_lat); 
    % ����������ά�ȣ�1440��101��time �� 101��1440��time��
    target_cloud_hcc = permute(era_cloud_hcc(:,:,d2h), [2,1]);
    target_cloud_mcc = permute(era_cloud_mcc(:,:,d2h), [2,1]);

    proj = projcrs(3413);  % WGS84��������ͶӰ
    [era_x, era_y] = projfwd(proj, era_lat_grid, era_lon_grid);  % ����101��1440����

    % ��ֵ��MODIS��Ԫ
    F_hcc = scatteredInterpolant(era_x(:), era_y(:), target_cloud_hcc(:), 'linear', 'none');
    hcc_interp = F_hcc(myd_x, myd_y);  % myd_x��myd_y����ǰͶӰ
    M_hcc = scatteredInterpolant(era_x(:), era_y(:), target_cloud_mcc(:), 'linear', 'none');
    mcc_interp = M_hcc(myd_x, myd_y);  % myd_x��myd_y����ǰͶӰ
    inx_cc = hcc_interp > 0.33 | mcc_interp > 0.51;
    mask_cloud_nan(inx_cc) = nan;
    % ������û�д�������ݸ�ֵnan
    mask_cloud_nan(:,num_cols*10+1:end) = nan;
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

function mod2sit = ist2sit(surf_aux,mod_aux,nswf)
    % ע�������surf_auxֻ��һ������
    a = 5.67e-8;        % ������������
    p = 1.3;                % �����ܶ�
    C = 1.44e-3;            % ������Ǳ�������崫��ϵ��
    Li = 2.86e6;
    
    % ���޳���ı����¶�ת��Ϊ�������
    R = 6371; % ����뾶��km��
    % ת��Ϊ��ά����
    x = R .* cosd(mod_aux(:,2)) .* cosd(mod_aux(:,1));
    y = R .* cosd(mod_aux(:,2)) .* sind(mod_aux(:,1));
    z = R .* sind(mod_aux(:,2));
    era_coords = [x, y, z];
    tree = KDTreeSearcher(era_coords); % ����KDTree

    surf_lon = surf_aux(1);
    surf_lat = surf_aux(2);
    x_surf = R .* cosd(surf_lat) .* cosd(surf_lon);
    y_surf = R .* cosd(surf_lat) .* sind(surf_lon);
    z_surf = R .* sind(surf_lat);
    surf_coords = [x_surf, y_surf, z_surf];
    [idx, ~] = knnsearch(tree, surf_coords, 'K', 1); % ������ѯ�����
    surf_aux(4:5) = mod_aux(idx, 3:4); % ��ȡ�¶Ⱥͷ���
    inx = surf_aux(4)>273.15-5;
    surf_aux(inx,:)=[];

    lwd_u = 0.99*a*surf_aux(3).^4;   % ���� ���� 
%         e = 0.7855*(1+0.2232*0.4^0.75);
    lwd_d = 0.7855*a*surf_aux(4).^4; % ��������     
    % ����ͨ��
    Fs = p*1004*C*surf_aux(5)*(surf_aux(4) - surf_aux(3));
    % Ǳ��ͨ�� Fe
%         Fe_MXY = 0.622*p*Li*C*sort_surf_T(:,5).*SVP_MXY(1000,sort_surf_T(:,6)-273.15,sort_surf_T(:,3)-273.15)./1000; 
    % ��һ�� Ǳ��ͨ�����㷽ʽ
%         Fe = 0.622*p*2.49e3*C*S_Loca(:,5).*(0.9*SVP(S_Loca(:,4))-SVP(S_Loca(:,3)))./1000;
    % �����ּ���Ǳ��ͨ�����㷽ʽ
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
