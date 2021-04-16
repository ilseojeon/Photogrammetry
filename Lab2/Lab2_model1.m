close all;
% clear all;
clc;
format long g
% Ground Coverage Visualization with Photogrammetric Model
% Ilseo Jeon
% University of Seoul
% Adv. Photogrammetry Lab2 in Master Course
% 2020. 04. 10

%% Read a .txt file including Exif tags
fid = fopen( 'mavic2_eoini_tm.txt', 'r');

% % Easting, Northing, Altitude, Roll, Pitch, Yaw 
% % (body) roll pitch yaw (gimbal) roll pitch yaw (camera) roll pitch yaw
imgexif = textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'CommentStyle','#','Delimiter','\t');
fclose(fid);
no_data = size(imgexif{2},1); % number of files

%% Pos/Ot Parameters Initialization
C_WB = zeros(no_data,3); % lon, lat, alt
A_WB = zeros(no_data,3); % body roll pitch yaw
A_WG = zeros(no_data,3); % gimbal roll pitch yaw

mat_WB = cell(no_data,1); % BCS(body coordinate system) ot in WCS(world coordinate system)
mat_WG = cell(no_data,1); % GCS(gimbal coordinate system) ot in WCS(world coordinate system)
mat_WC = cell(no_data,1); % CCS(camera coordinate system) ot in WCS(world coordniate system)

     % X Y Z   definition of NED coordinate system
NED = [0 1 0; ...
       1 0 0; ...
       0 0 -1];
    % input data for each variable
for n = 1:no_data
    C_WB(n,:) = [imgexif{2}(n) imgexif{3}(n) imgexif{4}(n)];  
    C_WB(n,1:2) = ll2utm(C_WB(n,2), C_WB(n,1)); % WGS84 to TM
    A_WB(n,:) = [imgexif{8}(n) imgexif{9}(n) imgexif{10}(n)];
    A_WG(n,:) = [imgexif{11}(n) imgexif{12}(n) imgexif{13}(n)]; 
    
    mat_WB{n} = NED*A2R_YPR(A_WB(n,:)*pi/180); %body rpy to matrix; represents bcs in WCS
    mat_WG{n} = NED*A2R_YPR(A_WG(n,:)*pi/180); %gimbal rpy to matrix; represents gcs in WCS
%     R_CW{n} = [0 0 1; 0 1 0; -1 0 0]*[1 0 0; 0 0 1; 0 -1 0]*mat_WG{n}; %  pitch 90 roll -90
    mat_WC{n} = [0 0 -1; 0 1 0; 1 0 0]*[1 0 0; 0 0 1; 0 -1 0]*mat_WG{n}; %  represents ccs in wcs % pitch -90 roll -90 
end

%% Image Initialization
%Drone Model: Mavic Pro
p_size = (6.3/4000)/1000; % pixel size; 6.3 means image width (mm) -> unit in (m)
f= (4.73)/1000; % focal length (mm) -> unit in (m)
f_mat = [0 0 f]; % focal length representation with matrix
i_width = 4000;
i_height = 3000;

ul_ics = [0 0 0];
dl_ics = [0 i_height 0];
ur_ics = [i_width 0 0];
dr_ics = [i_width i_height 0];

img_origin=[i_width/2 i_height/2 0];

% 영상의 3차원 지상좌표에서의 위치
% 영상 바운더리 4개의 위치에서 x,y,z로 3개의 값을 가짐
% (up left)   (up right)
% |ㅡㅡㅡㅡㅡㅡㅡㅡ|
% |     X     |
% |ㅡㅡㅡㅡㅡㅡㅡㅡ|
% (down left) (down right)

% ICS -> CCS
ul_pcs = zeros(no_data, 3);
dl_pcs = zeros(no_data, 3);
ur_pcs = zeros(no_data, 3);
dr_pcs = zeros(no_data, 3);
R_IC = [1 0 0; 0 -1 0; 0 0 1];
%PCS -> GCS
%Image in GCS
ul_gcs = zeros(no_data, 3);
dl_gcs = zeros(no_data, 3);
ur_gcs = zeros(no_data, 3);
dr_gcs = zeros(no_data, 3);

for n=1:no_data
    ul_pcs(n,:) = (R_IC*((ul_ics-img_origin)*p_size-f_mat)')';
    dl_pcs(n,:) = (R_IC*((dl_ics-img_origin)*p_size-f_mat)')';
    ur_pcs(n,:) = (R_IC*((ur_ics-img_origin)*p_size-f_mat)')';
    dr_pcs(n,:) = (R_IC*((dr_ics-img_origin)*p_size-f_mat)')';

    ul_gcs(n,:) = (mat_WC{n}*ul_pcs(n,:)'+C_WB(n,:)')';
    dl_gcs(n,:) = (mat_WC{n}*dl_pcs(n,:)'+C_WB(n,:)')';
    ur_gcs(n,:) = (mat_WC{n}*ur_pcs(n,:)'+C_WB(n,:)')';
    dr_gcs(n,:) = (mat_WC{n}*dr_pcs(n,:)'+C_WB(n,:)')';
end

%% Image projection onto ground: setting scale with the collinearity model
dem = 0; % unit: m
scale = zeros(no_data,1);

gp_u_l = zeros(no_data, 3);
gp_d_l = zeros(no_data, 3);
gp_u_r = zeros(no_data, 3);
gp_d_r = zeros(no_data, 3);

for n=1:no_data
    scale(n,1) = (dem - C_WB(n,3))/(-f);

    gp_u_l(n,1:2) = scale(n,1)*ul_pcs(n,1:2);
    gp_u_l(n,:)=(mat_WC{n}*gp_u_l(n,:)'+C_WB(n,:)')';
    gp_u_l(n,3) = dem;
    gp_d_l(n,1:2) = scale(n,1)*dl_pcs(n,1:2);
    gp_d_l(n,:)=(mat_WC{n}*gp_d_l(n,:)'+C_WB(n,:)')';
    gp_d_l(n,3) = dem;
    gp_u_r(n,1:2) = scale(n,1)*ur_pcs(n,1:2);
    gp_u_r(n,:)=(mat_WC{n}*gp_u_r(n,:)'+C_WB(n,:)')';
    gp_u_r(n,3) = dem;
    gp_d_r(n,1:2) = scale(n,1)*dr_pcs(n,1:2);
    gp_d_r(n,:)=(mat_WC{n}*gp_d_r(n,:)'+C_WB(n,:)')';
    gp_d_r(n,3) = dem;
end

%% Visualization
%Initialization for Visualization
img_x=cell(no_data, 1);
img_y=cell(no_data, 1);
img_z=cell(no_data, 1);
imggp_x=cell(no_data, 1);
imggp_y=cell(no_data, 1);
imggp_z=cell(no_data, 1);

figure ()
hold on
offset = 10;
cs_scale=5; % coordinate system scale
C_WB_tmp=zeros(no_data,3);
%%NED Coordinate System Visualization
vis_coord_system(C_WB(1,:)'+20, NED', 30, sprintf("NED"),'r');

%%Flight, Gimbal, Camera Coordinate System Visualization
for n = 1:no_data
    %% FCS Visualization
    C_WB_tmp(n,:)=C_WB(n,:);
    C_WB_tmp(n,3)=C_WB(n,3)+offset; 
    vis_coord_system(C_WB_tmp(n,:)', mat_WB{n}', cs_scale, sprintf('body %d',n),'r');
    %% GCS Visualization
    C_WB_tmp(n,3)=C_WB_tmp(n,3)+offset;
    vis_coord_system(C_WB_tmp(n,:)', mat_WG{n}, cs_scale, sprintf('gimbal %d',n),'r');
    %% CCS Visualization
    vis_coord_system(C_WB(n,:)', mat_WC{n}', cs_scale, sprintf('camera %d',n),'r');
    %% Image Frame Visualization
    hold on
    img_x{n} = [ul_gcs(n,1) ur_gcs(n,1) dr_gcs(n,1) dl_gcs(n,1) ul_gcs(n,1)]';
    img_y{n} = [ul_gcs(n,2) ur_gcs(n,2) dr_gcs(n,2) dl_gcs(n,2) ul_gcs(n,2)]';
    img_z{n} = [ul_gcs(n,3) ur_gcs(n,3) dr_gcs(n,3) dl_gcs(n,3) ul_gcs(n,3)]';    
    plot3(img_x{n}, img_y{n}, img_z{n}, 'b-', 'LineWidth', 2);
    %% Ground Coverage Visualization
    hold on
    imggp_x{n} = [gp_u_l(n,1) gp_u_r(n,1) gp_d_r(n,1) gp_d_l(n,1) gp_u_l(n,1)]';
    imggp_y{n} = [gp_u_l(n,2) gp_u_r(n,2) gp_d_r(n,2) gp_d_l(n,2) gp_u_l(n,2)]';
    imggp_z{n} = [gp_u_l(n,3) gp_u_r(n,3) gp_d_r(n,3) gp_d_l(n,3) gp_u_l(n,3)]';    
    plot3(imggp_x{n}, imggp_y{n}, imggp_z{n}, 'g-', 'LineWidth', 2);
    %% lines from camera center to each image corner point
    plot3([C_WB(n,1), gp_u_l(n,1)], [C_WB(n,2),gp_u_l(n,2)], [C_WB(n,3),gp_u_l(n,3)], 'k');
    plot3([C_WB(n,1), gp_d_l(n,1)], [C_WB(n,2),gp_d_l(n,2)], [C_WB(n,3),gp_d_l(n,3)], 'k');
    plot3([C_WB(n,1), gp_u_r(n,1)], [C_WB(n,2),gp_u_r(n,2)], [C_WB(n,3),gp_u_r(n,3)], 'k');
    plot3([C_WB(n,1), gp_d_r(n,1)], [C_WB(n,2),gp_d_r(n,2)], [C_WB(n,3),gp_d_r(n,3)], 'k');
end
view(3)
grid on, axis equal
title('Ground Coverage with Photorgrammetric Model') 
xlabel('Easting'), ylabel('Northing'), zlabel('Z')