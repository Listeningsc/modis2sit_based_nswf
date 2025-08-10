% author:LJX 
% version: I
% date:  2024/08/16
% statement :    查看L2数据和云掩膜数据是否对应
%              I 没有考虑太阳短波辐射，辅助数据则是距离该点处17.5km的AUX处数据
%                 该平衡方程各个参数都是根据Maykut2013，和AVHRR 1996的文章得出

era_mat_path = 'Z:\ERA5_aux\nc2mat\';
era_path = 'Z:\ERA5_aux';
ib_parent_path = 'W:\ljx_aux\IceBridge\IDCSI4\ib_files\';
% interest_region中想要看的是博福特海域附近的数据
ist_parent_path = 'W:\ljx_aux\MODIS\data\Interest_region\L2_IST';
ist_coord_parent_path = 'W:\ljx_aux\MODIS\data\Interest_region\L2_coord';
cloud_parent_path  = 'W:\ljx_aux\MODIS\data\Interest_region\Cloud_4L2';

years = [2011:2015,2017:2019]; % 这是IB数据存在的年份

[~,ib_dates_set] = IB_date(ib_parent_path,years);

tic 
sensor = 'MOD';     % 这里指的是搭载在不同卫星上的MODIS传感器
night = 0;          % 这里判断是白天还是晚上，白天的话需要计算短波辐射。

% 选取特定的表面温度hdf
for sub_ib_dates = ib_dates_set
    % 当年的IB数据日期
    ib_date_mat = cell2mat(sub_ib_dates);
    if isempty(ib_date_mat)
        continue
    else
        for target_ib_date = ib_date_mat
            % 当年内某一天的IB数据
            str_ib_date = num2str(target_ib_date);
            target_year = str_ib_date(1:4);
            
            datetime_date = datetime(str_ib_date, 'InputFormat', 'yyyyMMdd');  % 转换为 datetime
            target_ib_day = day(datetime_date,'dayofyear');
            result_doy = cell2mat(strcat(target_year,compose('%03d',target_ib_day)));
            
            [~,sit_day] = sit2thermodynamic_equilibrium(result_doy,ist_parent_path,ist_coord_parent_path,cloud_parent_path,...
                                                               sensor,night,target_year,era_path,era_mat_path,str_ib_date);
            fprintf('success in %d',target_ib_date);
                                                                                    
        end
    end
end
toc
