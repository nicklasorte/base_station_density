clear;
clc;
close all force;
close all;
app=NaN(1);  %%%%%%%%%This is to allow for Matlab Application integration.
format shortG
%format longG
top_start_clock=clock;
folder1='C:\Users\nlasorte\OneDrive - National Telecommunications and Information Administration\MATLAB2024\7GHz FSS Neighborhoods';
cd(folder1)
addpath(folder1)
addpath('C:\Users\nlasorte\OneDrive - National Telecommunications and Information Administration\MATLAB2024\Basic_Functions')
pause(0.1)
%load('us_cont.mat','us_cont')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'Increase the density of the Randomized Real: x2, x3, x4'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Read in randomized real data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Base Station Deployment
%%%%Deployment Lat/Lon
tf_repull_rand=0%1%0
excel_filename_rand='Rand_Real_2025_1Sector_idx.xlsx' %%%%%(This is the 1 sector and nationwide idx)
mat_filename_str_rand=strcat('rand_real_2025_one_sector_idx.mat')
tic;
[cell_rand_real]=load_full_excel_rev1(app,mat_filename_str_rand,excel_filename_rand,tf_repull_rand);
toc;
rand_real_2025=cell2mat(cell_rand_real([2:end],:)); %%%%%%%%1)Lat, 2)Lon, 3)Antenna Height 4)Azimuth 5)IDX
base_station_latlonheight=rand_real_2025;  %%1)Lat, 2)Lon, 3)Height meters, 4)Azimuth 5)Idx 6)EIRP IDX 7)Clutter IDX, 1==Urban, 2==Suburban, 3==Rural
base_station_latlonheight(1:10,:)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
array_latlon=base_station_latlonheight(:,[1,2,3]);  %%%%1) Latitude, 2) Longitude, 3)Height

%%%Just find the knn for each, find the distance between, divide it by 2,
%%%%Randomize within that circle
knn_idx=knnsearch(fliplr(array_latlon),fliplr(array_latlon),'K',3);

%%%Find the distance between array_latlon and
%%%array_latlon(knn_idx(:,2),:)
knn_latlon=array_latlon(knn_idx(:,2),:);
pair_dist_km=deg2km(distance(array_latlon(:,1),array_latlon(:,2),knn_latlon(:,1),knn_latlon(:,2)));

%%%%%Find the zero distance and put in the next
zero_idx=find(pair_dist_km==0)
knn_latlon(zero_idx,:)=array_latlon(knn_idx(zero_idx,3),:);

pair_dist_km=deg2km(distance(array_latlon(:,1),array_latlon(:,2),knn_latlon(:,1),knn_latlon(:,2)));
min(pair_dist_km)

%%%%%%Now split the difference for the radius
array_latlonrad=horzcat(array_latlon(:,[1,2]),pair_dist_km/3);

% tic;
% rand_pts=randomPointsForCircles(array_latlonrad(:,1), array_latlonrad(:,2), array_latlonrad(:,3), 1);
% toc;
% rand_latlon=vertcat(rand_pts{:});
% size(rand_latlon)
% size(array_latlon)
%%%%%%%%%x2 to x4
multiplier=4
cell_rand_sample=cell(multiplier,1);
cell_rand_sample{1}=array_latlon;
for i=1:1:multiplier-1
    %%%%%
    tic;
    [latP, lonP, id]=randCircumference_scircle1(array_latlonrad(:,1), array_latlonrad(:,2), array_latlonrad(:,3),i,50);
    toc;
    rand_latlon=horzcat(latP, lonP,array_latlon(id,3));
    if i==1
        rand_latlon1=rand_latlon;
    end
    cell_rand_sample{i+1}=vertcat(array_latlon,rand_latlon);
    %     rand_latlon(1:10,:)
    % unique(rand_latlon(:,3))
    % size(rand_latlon)
end
cell_rand_sample

save(strcat('cell_rand_sample_',num2str(multiplier),'.mat'),'cell_rand_sample')


%%%%%Write to excel
filename = strcat('cell_rand_sample_',num2str(multiplier),'.xlsx');
num_rows = numel(cell_rand_sample);
for k = 1:1:num_rows
    sheetName = "Sheet" + k;   % Create sheet names: Sheet1, Sheet2, ...
    tic;
    %writematrix(cell_rand_sample{k}, filename,Sheet=sheetName);
    toc;
end

% Write KML file (NO labels)
kml_latlon=cell_rand_sample{2};
[num_pts,~]=size(kml_latlon)
names = strings(num_pts,1);  % empty labels
tic;
%kmlwritepoint(strcat('UpsampleRandReal2025_x2.kml'), kml_latlon(:,1),kml_latlon(:,2),Name=names,IconScale=1);
toc; %%%%36 seconds for x2




%%%%%%%%%%%%%%%%%Show an example of Nats Park
%base_polygon=horzcat(39.617,-74.82)
base_polygon=horzcat(38.873,-77.0074)
sim_radius_km=10;
[sim_bound]=calc_circle_bound(app,base_polygon,sim_radius_km);

%%%%%%Find the points inside

[inside_idx1]=find_points_inside_contour_two_step(app,sim_bound,array_latlon);
[inside_idx2]=find_points_inside_contour_two_step(app,sim_bound,rand_latlon1);
 
f1=figure;
geoplot(sim_bound(:,1),sim_bound(:,2),'-k','LineWidth',3)
hold on;
geoplot(array_latlon(inside_idx1,1),array_latlon(inside_idx1,2),'ob','LineWidth',2)
geoplot(rand_latlon1(inside_idx2,1),rand_latlon1(inside_idx2,2),'xr','LineWidth',2)
geobasemap streets-light%landcover
pause(0.1)
f1.Position = [100 100 1200 900];
pause(1)
filename1=strcat('NatsPark_Density.png');
retry_save=1;
while(retry_save==1)
    try
        saveas(gcf,char(filename1))
        retry_save=0;
    catch
        retry_save=1;
        pause(1)
    end
end

