% Simulate the test data for AT
% Impyeong Lee
% 2009. 1. 12

clear all   % clear all the variable previsoully defined
close all   % close all the opend figure window
clc         % clear the command prompt window

path = 'C1\\';

% Interior Orientation
x0 = 0.001;     % x-coordinates of the priniciple point expressed in Fiducial Coordinate System (FCS)
y0 = -0.001;    % y-coordinates of the priniciple point expressed in Fiducial Coordinate System (FCS)
c = 17;    % the calibrated focal lenth (probably, in case of a wide angle photogrammetric cammera)

IO_t = [x0 y0 c];

% Photo size
pixel_size = [3.45 3.45]; % um
pixel_cnt = [2456 2058];
psz = pixel_cnt .* pixel_size * 1e-3; % mm

% Project Area
pa = [0 200 0 120];   % The range of the project area, [xmin xmax ymin ymax], expressed in absolute coordniate system

% Nominal terrain elevation
Za = 0;     %   average terrain elevation

% Nominal scale
plat_h = 200; % 200,400,600
scale = c/plat_h;
nsc = plat_h/c;
% actually, this means 5*1000, 1:5000
% because all the camera or photo related coordinates are expressed in mm,
% but all the ground related coordinates in meter.

% Nominal altitude
Zh = nsc * c + Za;
% the nominal altitude of the image acquisition in absolute height.

% no. strips and no. photos per strip
plat_vel = 36; % km
frame_rate = 0.5;
im_gnd_res = pixel_size .*nsc * 1e-3 ;
im_gnd_size = pixel_cnt .* im_gnd_res; 
dst_bwn_im = plat_vel/ 3600 / frame_rate* 1e3 ;
ol = (im_gnd_size(1) - dst_bwn_im) / im_gnd_size(1);       % overlap = 97.83%
sl = 0.3;       % sidelap = 20%
% photos per strips
ps = ceil( ( (pa(2) - pa(1) ) / nsc / psz(1) - 1) / (1 - ol) ) + 3;
% no. strips
ns = ceil(( (pa(4) - pa(3) ) / nsc / psz(2) - 1 ) / ( 1 - sl)) + 1;

% no. images
no_IM = ps * ns;

% Determine the position of the first photograph, locating at the top-left
% postion of the project area
EO_f(1) = pa(1) + psz(1) * nsc * (ol - 0.5);
EO_f(2) = pa(4) - psz(2) * nsc * 0.5;

std_plat_pos = 1;
std_plat_att = 1 * pi / 180;

% EOs of each photograph
for n = 1:ps,
    for m = 1:ns,
        id = (n-1) * ns + m;
        EO_t(id,1) = (n-1) * (1-ol) * psz(1) * nsc + EO_f(1) + randn(1) * std_plat_pos;
        EO_t(id,2) = - (m-1) * (1-sl) * psz(2) * nsc + EO_f(2) + randn(1) * std_plat_pos;
        EO_t(id,3) = Zh + randn(1) * std_plat_pos;
        EO_t(id,4:6) = [0 0 0] + randn(1,3) * std_plat_att;
    end
end

% Image corner points
bipf = [-psz(1) psz(2); psz(1) psz(2); psz(1) -psz(2); -psz(1) -psz(2)] / 2;
bip(:,1) = bipf(:,1) - x0;  % x: converting to PCS
bip(:,2) = bipf(:,2) - y0;  % y: converting to PCS

% Compute the image boundary on the ground
figure
hold on
for n = 1:no_IM,
    R = Rot3D( EO_t(n,4:6) );
    Pc = EO_t(n,1:3)';
    for nbip = 1:size(bip,1),
        p = [bip(nbip,1) bip(nbip,2) -c]';  % the vector from the PC to an image corner point expressed in PCS
        s(n, nbip) =  (Za - EO_t(n,3)) / ( R(:,3)' * p ); % determine the scale
        bgp{n}(nbip,:) = s(n, nbip) * R' * p + Pc;  % compute the ground points using the scale
    end
    h = patch(bgp{n}(:,1), bgp{n}(:,2), 'b')
    set( h, 'FaceColor', 'b', 'FaceAlpha', 0.05, 'EdgeColor', 'b' );
    h = text ( bgp{n}(1,1),  bgp{n}(1,2), sprintf('%d', n) );
    set( h, 'FontWeight', 'bold', 'Color', 'b');
end
grid on
axis equal
hold off

% Determine the entire ground coverage
tmp = cell2mat(bgp');
egp = [min(tmp(:,1)) max(tmp(:,1)) min(tmp(:,2)) max(tmp(:,2))];

% Generate the ground points
gdst = nsc * psz .* 0.24;     % interval between two adjacent ground points
gofs = nsc * psz .* [0.30 0.14];     % the position of the first ground point (top-left)
ngpx = ceil ( ( egp(2) - egp(1) - 2 * gofs(1) ) / gdst(1) );
gdst(1) = ( egp(2) - egp(1) - 2 * gofs(1) ) / (ngpx-1);
ngpy = ceil ( ( egp(4) - egp(3) - 2 * gofs(2) ) / gdst(2) );
gofs(2) = (  egp(4) - egp(3) - (ngpy-1) * gdst(2) ) / 2;
for m = 1:ngpx,
    for n = 1:ngpy,
        idx = (m-1) * ngpy + n;
        GP_t(idx,:) = [(m-1) * gdst(1) + egp(1)+gofs(1) egp(4) - (n-1) * gdst(2) - gofs(2) ];
    end
end
no_GP = size(GP_t,1);
% GP_t(:,3) = 40 * rand(no_GP,1) - 20;
GP_t(:,3) = 10 * rand(no_GP,1) - 5;

% plot the ground points
figure
plot3( GP_t(:,1), GP_t(:,2), GP_t(:,3), 'ro')
axis equal
xlabel('x'), ylabel('y'), zlabel('z')
title('Plot of ground points');
grid on
view(2)

% Generate the true image points
m = 0;
for ni = 1:no_IM,
    R = Rot3D( EO_t(ni,4:6) );    % compute the rotational matrix
    for np = 1:no_GP,
        GC = GP_t(np,:)' - EO_t(ni,1:3)';   % Ground point - Camera position
        ND = R * GC;    % ND: the numerator and denominator of the the collinearity equation                    
        pos = [x0 y0]' - c / ND(3) * ND(1:2,1);     % the collinearity equation
        if pos(1) >= -psz(1)/2 && pos(1) <= psz(1)/2 && pos(2) >= -psz(2)/2 && pos(2) <= psz(2)/2,
            % Check if the image point locate within the range of the image
            m = m + 1;
            IP_t(m,1) = ni;
            IP_t(m,2) = np;
            IP_t(m,3:4) = pos';
        end
    end
end
no_ip = size(IP_t,1);

% Generate the noisy image points
std_ip = 0.00345; % 5 um
IP_m = [IP_t(:,1:2) IP_t(:,3:4) + randn(size(IP_t,1),2) * std_ip]; % Add the random noises associated with the measurement of the image point

% Generate the noisy ground points for ground control points
std_gcp = 0.05; % 5 cm
GP_m = GP_t + randn(no_GP,3) * std_gcp;

% No error of EO using GPS/INS
%std_gps = 0.001; % 1 mm
%std_ins = 5e-6;  % 1 mm / 200 m

%% seperated GPS/IMU
std_gps = 0.3; % 30 cm
std_ins = 0.1 * pi / 180;  % 0.1 deg

% Direct Observations of EO using GPS/INS
%std_gps = 2.25; % mems 2 ~ 2.5 m
%std_ins = 2 * pi / 180;  % mems 2 deg

EO_m(:,1:3) = EO_t(:,1:3) + randn(size(EO_t,1),3) * std_gps;
EO_m(:,4:6) = EO_t(:,4:6) + randn(size(EO_t,1),3) * std_ins;

% Compute the statistics of measured IP difference
for n = 3:4,
    IP_m_stat(:,n-2) = comp_stat(IP_m(:,n)-IP_t(:,n))*1e3;
end

% Compute the statistics of measured EO difference
for n = 1:6,
    EO_m_stat(:,n) = comp_stat(EO_m(:,n)-EO_t(:,n));
end
EO_m_stat(:,1:3) = EO_m_stat(:,1:3) * 1e3;
EO_m_stat(:,4:6) = EO_m_stat(:,4:6) * 180 / pi;

% Compute the statistics of measure GP difference
for n = 1:3,
    GP_m_stat(:,n)  = comp_stat(GP_m(:,n)-GP_t(:,n))*1e3;
end

% Store AT configuration
mkdir ( path );
fid = fopen(strcat(path, 'Sim_summary.txt'), 'w');
fprintf(fid, 'Interior Orientation (xo, yo, c) [mm]: %g, %g, %g\r\n', IO_t );
fprintf(fid, 'No. pixels: %d x %d\r\n', pixel_cnt );
fprintf(fid, 'Pixel size [um]: %g x %g\r\n', pixel_size );
fprintf(fid, 'Image size [mm]: %g x %g\r\n', psz );
fprintf(fid, 'Average terrain height [m]: %g\r\n', Za );
fprintf(fid, 'Platform height from the terrain [m]: %g\r\n', plat_h );
fprintf(fid, 'Image scale : %g\r\n', nsc * 1e3 );
fprintf(fid, 'Image ground resolution [m]: %g x %g\r\n', im_gnd_res );
fprintf(fid, 'Image ground size [m]: %g x %g\r\n', im_gnd_size );
fprintf(fid, 'Project Area (x_min, x_max, y_min, y_max) [m]: %g ~ %g, %g ~ %g\r\n', pa );
fprintf(fid, 'Platform velocity [km/h]: %g\r\n', plat_vel );
fprintf(fid, 'Frame rate [image/s]: %d\r\n', frame_rate );
fprintf(fid, 'Distance between images [m]: %g\r\n', dst_bwn_im );
fprintf(fid, 'Overlap ratio [%%]: %g\r\n', ol );
fprintf(fid, 'Sidelap ratio [%%]: %g\r\n', sl );
fprintf(fid, 'No. photos per strip: %d\r\n', ps );
fprintf(fid, 'No. strips: %d\r\n', ns );
fprintf(fid, 'No. images: %d\r\n', no_IM );
fprintf(fid, 'No. ground points: %d\r\n', no_GP );
fprintf(fid, 'No. image points: %d\r\n', no_ip );
fprintf(fid, 'Standard deviation of image point measurement error [mm]: %g\r\n', std_ip );
fprintf(fid, 'Standard deviation of ground point measurement error [m]: %g\r\n', std_gcp );
fprintf(fid, 'Standard deviation of GPS measurement errors [m]: %g\r\n', std_gps );
fprintf(fid, 'Standard deviation of INS measurement errors [rad]: %g\r\n', std_ins );
fprintf(fid, 'Standard deviation of INS measurement errors [deg]: %g\r\n', std_ins*180/pi );
fprintf ( fid, 'Statistics for IP_m - IP_t\r\n' );
fprintf ( fid, '        %11s\t%11s\r\n', 'x[um]', 'y[um]' );
fprintf ( fid, 'Minimum:%11.3f\t%11.3f\r\n', IP_m_stat(1,:) );
fprintf ( fid, 'Maximum:%11.3f\t%11.3f\r\n', IP_m_stat(2,:) );
fprintf ( fid, 'Average:%11.3f\t%11.3f\r\n', IP_m_stat(3,:) );
fprintf ( fid, 'Std_dv.:%11.3f\t%11.3f\r\n', IP_m_stat(4,:) );
fprintf ( fid, 'RMS    :%11.3f\t%11.3f\r\n', IP_m_stat(5,:) );
fprintf ( fid, 'Statistics for EO_m - EO_t\r\n' );
fprintf ( fid, '        %11s\t%11s\t%11s\t%11s\t%11s\t%11s\r\n', 'Xc[mm]', 'Yc[mm]', 'Zc[mm]', 'Omega[deg]', 'Phi[deg]', 'Kappa[deg]' );
fprintf ( fid, 'Minimum:%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_m_stat(1,:) );
fprintf ( fid, 'Maximum:%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_m_stat(2,:) );
fprintf ( fid, 'Average:%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_m_stat(3,:) );
fprintf ( fid, 'Std_dv.:%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_m_stat(4,:) );
fprintf ( fid, 'RMS    :%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_m_stat(5,:) );
fprintf ( fid, 'Statistics for GP_m - GP_t\r\n' );
fprintf ( fid, '        %11s\t%11s\t%11s\r\n', 'X[mm]', 'Y[mm]', 'Z[mm]' );
fprintf ( fid, 'Minimum:%11.3f\t%11.3f\t%11.3f\r\n', GP_m_stat(1,:) );
fprintf ( fid, 'Maximum:%11.3f\t%11.3f\t%11.3f\r\n', GP_m_stat(2,:) );
fprintf ( fid, 'Average:%11.3f\t%11.3f\t%11.3f\r\n', GP_m_stat(3,:) );
fprintf ( fid, 'Std_dv.:%11.3f\t%11.3f\t%11.3f\r\n', GP_m_stat(4,:) );
fprintf ( fid, 'RMS    :%11.3f\t%11.3f\t%11.3f\r\n', GP_m_stat(5,:) );
fclose(fid);

% Store the true IO
fid = fopen(strcat(path, 'IO_t.txt'), 'w');
fprintf(fid, '%g\r\n', IO_t );
fclose(fid);

% Store the true EO
fid = fopen(strcat(path, 'EO_t.txt'), 'w');
fprintf(fid, '%d\t%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', [(1:no_IM)' EO_t]' );
fclose(fid);

% Store the noisy EO
fid = fopen(strcat(path, 'EO_m.txt'), 'w');
fprintf(fid, '%d\t%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', [(1:no_IM)' EO_m]' );
fclose(fid);

% Store the true ground points
fid = fopen(strcat(path, 'GP_t.txt'), 'w');
fprintf(fid, '%d\t%11.3f\t%11.3f\t%11.3f\r\n', [(1:no_GP)' GP_t]' );
fclose(fid);

% Store the noisy ground points
fid = fopen(strcat(path, 'GP_m.txt'), 'w');
fprintf(fid, '%d\t%11.3f\t%11.3f\t%11.3f\r\n', [(1:no_GP)' GP_m]' );
fclose(fid);

% Store the true image points
fid = fopen(strcat(path, 'IP_t.txt'), 'w');
fprintf(fid, '%d\t%d\t%11.4f\t%11.4f\r\n', IP_t' );
fclose(fid);

% Store the noisy image points
fid = fopen(strcat(path, 'IP_m.txt'), 'w');
fprintf(fid, '%d\t%d\t%11.4f\t%11.4f\r\n', IP_m' );
fclose(fid);

% Plot the ground coverage
figure
hold on
%h = patch ( [pa(1); pa(1); pa(2); pa(2)], [pa(3); pa(4); pa(4); pa(3)], 'r' );
%set( h, 'FaceAlpha', 0.3, 'EdgeColor', 'r' );
for n = 1:no_IM,
    h = patch( bgp{n}(:,1), bgp{n}(:,2), 'b' );
    set( h, 'FaceColor', 'b', 'FaceAlpha', 0.05, 'EdgeColor', 'b' );
    h = text ( bgp{n}(1,1), bgp{n}(1,2), sprintf('%d', n) );
    set( h, 'Color', 'b' );
end
%h = plot ( [egp(1); egp(1); egp(2); egp(2); egp(1)], [egp(3); egp(4); egp(4); egp(3); egp(3)], 'm' );
%set(h, 'LineWidth', 2);
for n = 1:no_GP,
    h = plot ( GP_t(n,1), GP_t(n,2), 'ro' );
    set ( h, 'MarkerFaceColor', 'r' );
    h = text ( GP_t(n,1), GP_t(n,2), sprintf('%d', n) );
end
grid on
xlabel('x')
ylabel('y')
axis equal
