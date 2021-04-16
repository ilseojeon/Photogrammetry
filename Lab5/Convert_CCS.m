% Converting IF in ICS into CCS
% Ilseo Jeon
% Dept. of Geoinformatics, The University of Seoul
% Adv. Photogrammetry Lab4 in Master Course
% 2020. 05. 28

% ccs; camera coordinate system
% ics; image coordinate system

clear all
close all
clc

R_IC=[1 0; 0 -1];
img_origin =[4000/2 3000/2 0]';

% unit mm
p_size=0.00153452;

% read txt
fid = fopen( 'TP_ics.txt','r');
fgets(fid); fgets(fid);
TP_ics = fscanf(fid, '%f', [5 inf] )';
fclose(fid);
no_TP = size(TP_ics,1);

% Read the image
im1 = imread( 'DJI_0007.jpg', 'jpg');
im2 = imread( 'DJI_0008.jpg', 'jpg');

% Visualize the image point over image 1
figure(1)
imshow(im1);
hold on
axis on

figure(2)
imshow(im2);
hold on
axis on

for m = 1:2
    figure(m)
    for n = 1:no_TP
        plot(TP_ics(n,m*2), TP_ics(n,m*2+1), 'rx', 'MarkerSize', 20, 'LineWidth', 2);
        h = text(TP_ics(n,m*2)+5, TP_ics(n,m*2+1), sprintf('%d', n) );
        set(h, 'Color', 'r');
    end
end

%Converting ICS into CCS
TP_ccs = zeros(size(TP_ics));

TP_ccs(:, 1) = TP_ics(:,1);    %Index
for n=1:no_TP
    TP_ccs(n, 2:3) = R_IC*[TP_ics(n,2)-img_origin(1) TP_ics(n,3)-img_origin(2)]'*p_size;
    TP_ccs(n, 4:5) = R_IC*[TP_ics(n,4)-img_origin(1) TP_ics(n,5)-img_origin(2)]'*p_size;
end

figure
axis equal
for i = 1:no_TP
    hold on
    plot(TP_ccs(i, 2), TP_ccs(i, 3), 'rx', 'MarkerSize', 20, 'LineWidth', 2);
    xlim([-3.15, 3.15])
    ylim([-2.3625, 2.3625])
    
end

figure
for i = 1:no_TP
    hold on
    plot(TP_ccs(i, 4), TP_ccs(i, 5), 'rx', 'MarkerSize', 20, 'LineWidth', 2);
    xlim([-3.15, 3.15])
    ylim([-2.3625, 2.3625])
    
end

fid = fopen ('TP_ccs.txt', 'w');
fprintf(fid, '%% TP Coordinates in FCS\r\n' );
fprintf(fid, '%d\r\n', size(TP_ccs,1) );
fprintf(fid, '%d\t%.6f\t%.6f\t%.6f\t%.6f\r\n', TP_ccs' );
fclose(fid);   

