target_path = [pwd,'\'];
hdf_list = dir([target_path,'*.','hdf']);
myd_info = hdfinfo(hdf_list(1).name);

myd35_lat = hdfread(hdf_list(1).name,'Geolocation Fields/Latitude');
myd35_lon = hdfread(hdf_list(1).name,'Geolocation Fields/Longitude');
myd35_cloud = hdfread(hdf_list(1).name,'Data Fields/Cloud_Mask');
cloud1 = myd35_cloud(1,:,:);
cloud1_2d = reshape(cloud1,2040,1354);

cloud2 = myd35_cloud(2,:,:);
cloud2_2d = reshape(cloud2,2040,1354);

cloud3 = myd35_cloud(3,:,:);
cloud3_2d = reshape(cloud3,2040,1354);

cloud4 = myd35_cloud(4,:,:);
cloud4_2d = reshape(cloud4,2040,1354);

cloud5 = myd35_cloud(5,:,:);
cloud5_2d = reshape(cloud5,2040,1354);

cloud6 = myd35_cloud(6,:,:);
cloud6_2d = reshape(cloud6,2040,1354);