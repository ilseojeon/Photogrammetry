clear all
close all
clc
format long g

% Lab 4 'Single Photo Resection'
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 1. Initialization of GCP, image
% 2. Measuring IPs
% 3. Calculatiing initial approximates of camera EOs
%       Input: GCP, IP, IO
% 4. SPR
% 5. Analysis of statistical results
% 6. Visualization of everything that we need for showing results
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Ilseo Jeon
% Dept. of Geoinformatics, The University of Seoul
% Adv. Photogrammetry Lab4 in Master Course
% 2020. 05. 8

%% 1. Initialization of GCP, an image
im = "DJI_0021.jpg";
GCP = [205000.384 553767.153 62.837;...
       205000.407 553762.250 62.952;...
       204998.115 553767.156 62.828;
       204998.120 553762.259 62.941;...
       204990.079 553780.005 69.556;...
       204990.075 553780.013 66.118;...
       205009.231 553779.912 66.118]';

% % Initialization of camera IOs
[no_data, no_var] = size(GCP);
resol = [5427 3648];
u0 = resol(1)/2;
v0 = resol(2)/2;

p_size =  0.00153452/1000; %unit: m/pixel
f= (8.8)/1000; %unit : m

ku = 1/p_size; % unit : pixel/m
kv = ku;
fu = f*ku;
fv = f*kv;

IO = [u0 v0 fu p_size]';
%% 2. Measuring IPs
% [IP_ccs, IP_ics]=M_IP(im, IO, resol);
IP_ics = [2929.19407894737,2780.10197368421;
        2943.65460526316,3372.79276315789;
        2571.67434210526,2771.38486842105;
        2501.30592105263,3361.12171052632;
        1741.44407894737,1014.39144736842;
        1783.43750000000,1394.25986842105;
        3911.72039473684,1425.00328947368]';
IP_ics(2,:) = 3000 - IP_ics(2,:);
% IP_ccs = [0.000332369965131579,-0.00146688502144737,-0.00880000000000000;
%         0.000353691716710527,-0.00237992442144737,-0.00880000000000000;
%         -0.000217765550394737,-0.00145365988197368,-0.00880000000000000;
%         -0.000324858893552631,-0.00236050062881579,-0.00880000000000000;
%         -0.00149178058934211,0.00124490458881579,-0.00880000000000000;
%         -0.00142644234302632,0.000656900754078947,-0.00880000000000000;
%         0.00184158048065789,0.000615297090131579,-0.00880000000000000];

%% 3. Calculatiing initial approximates of camera EOs
%       Input: GCP, IP, IO
%       output: projection centre coordinations of camera x, y, z (in WCS)
%       , orientations om, ph, kp
% ini_EO=ini(GCP, IO, IP_ccs);
ini_EO=[204998.569946372,553775.411712265,114.221961241597,0,0,-0.0584551365782694]';

%% 4. SPR
delta = 1e-9; %tolerance
[EstiEO, deo, et, vct, veo, f0]=SPR(IP_ics, GCP, IO, ini_EO, delta);
coef=corrcoef(veo);
%% 5. Analysis
% on Lab4 document
figure ()
hold on
plot3(veo(1,:), veo(2,:), veo(3,:), 'r*');
view(3)

figure()
hold on
plot(1:15, deo(1,:), 'r-', 'LineWidth', 2);
plot(1:15, deo(2,:), 'g-', 'LineWidth', 2);
plot(1:15, deo(3,:), 'b-', 'LineWidth', 2);

figure()
hold on
plot(1:15, deo(4,:), 'r-', 'LineWidth', 2);
plot(1:15, deo(5,:), 'g-', 'LineWidth', 2);
plot(1:15, deo(6,:), 'b-', 'LineWidth', 2);

%% 6. Visualization
[no, no_it]=size(EstiEO);
img_x=cell(no_it,1);
img_y=cell(no_it,1);
img_z=cell(no_it,1);
scale=1000;

figure ()
hold on
% GCP point
plot3(GCP(1,:),GCP(2,:),GCP(3,:),'b*');
IO(4)=IO(4)*scale; %scale
IO(3)=-IO(3); % for scale
[gcs_ul, gcs_dl, gcs_dr, gcs_ur, R]=Cal_IF(GCP, EstiEO, IO, resol);

for n=1:no_it
    % ccs in gcs
%     vis_coord_system(EstiEO(1:3,n)'+5, NED', 30, sprintf("NED"),'r');
    vis_coord_system(EstiEO(1:3,n), R{n}', 5, sprintf('camera %d',n),'r');
    % image frame
    hold on
    img_x{n} = [gcs_ul(1,n) gcs_ur(1,n) gcs_dr(1,n) gcs_dl(1,n) gcs_ul(1,n)]';
    img_y{n} = [gcs_ul(2,n) gcs_ur(2,n) gcs_dr(2,n) gcs_dl(2,n) gcs_ul(2,n)]';
    img_z{n} = [gcs_ul(3,n) gcs_ur(3,n) gcs_dr(3,n) gcs_dl(3,n) gcs_ul(3,n)]';    
    plot3(img_x{n}, img_y{n}, img_z{n}, 'b-', 'LineWidth', 2);

    hold on
    % lines from camera center to each image corner point
    
    plot3([EstiEO(1,n),gcs_ul(1,n)], [EstiEO(2,n),gcs_ul(2,n)], [EstiEO(3,n),gcs_ul(3,n)], 'k', 'LineWidth', 2);
    plot3([EstiEO(1,n),gcs_dl(1,n)], [EstiEO(2,n),gcs_dl(2,n)], [EstiEO(3,n),gcs_dl(3,n)], 'k', 'LineWidth', 3);
    plot3([EstiEO(1,n),gcs_ur(1,n)], [EstiEO(2,n),gcs_ur(2,n)], [EstiEO(3,n),gcs_ur(3,n)], 'k', 'LineWidth', 3);
    plot3([EstiEO(1,n),gcs_dr(1,n)], [EstiEO(2,n),gcs_dr(2,n)], [EstiEO(3,n),gcs_dr(3,n)], 'k', 'LineWidth', 3);
end

for n=1:7
    
end
view(3)
grid on, axis equal
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
title('SPR Visualization') 
xlabel('Easting'), ylabel('Northing'), zlabel('Z')