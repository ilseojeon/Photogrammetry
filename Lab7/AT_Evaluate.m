% Evaluate the Estimation Results
% Impyeong Lee
% 2009. 1. 12

% Initialize the work space
clear all
close all
clc

path = 'KeKg\\';

% Read true values for EO
fid = fopen( strcat(path,'EO_t.txt'),'r');
EO_t = fscanf(fid, '%f', [7 inf] )';
fclose(fid);

% Read true values for EO
fid = fopen( strcat(path,'EO_i.txt'),'r');
EO_i = fscanf(fid, '%f', [7 inf] )';
fclose(fid);

% Read estimates for EO
fid = fopen( strcat(path,'EO_e.txt'),'r');
EO_e = fscanf(fid, '%f', [7 inf] )';
fclose(fid);

% Read initial approx. for GP
fid = fopen( strcat(path,'GP_i.txt'),'r');
GP_i = fscanf(fid, '%f', [4 inf] )';
fclose(fid);

% Read true values for GP
fid = fopen( strcat(path,'GP_t.txt'),'r');
GP_t = fscanf(fid, '%f', [4 inf] )';
fclose(fid);

% Read estimates for GP
fid = fopen( strcat(path,'GP_e.txt'),'r');
GP_e = fscanf(fid, '%f', [4 inf] )';
fclose(fid);

m = 0;
for n = 1:size(EO_e),
    id = find (EO_t(:,1)==EO_i(n,1) );
    if length(id)==1
        m = m + 1;
        EO_i_df(m,:) = [EO_i(n,1) EO_i(n,2:7) - EO_t(id,2:7)];
    end
end
EO_i_df(:,2:4) = EO_i_df(:,2:4) * 1e3;
EO_i_df(:,5:7) = EO_i_df(:,5:7) * 180 / pi;

m = 0;
for n = 1:size(EO_e),
    id = find (EO_t(:,1)==EO_e(n,1) );
    if length(id)==1
        m = m + 1;
        EO_e_df(m,:) = [EO_e(n,1) EO_e(n,2:7) - EO_t(id,2:7)];
    end
end
EO_e_df(:,2:4) = EO_e_df(:,2:4) * 1e3;
EO_e_df(:,5:7) = EO_e_df(:,5:7) * 180 / pi;

m = 0;
for n = 1:size(GP_i),
    id = find (GP_t(:,1)==GP_i(n,1) );
    if length(id)==1
        m = m + 1;
        GP_i_df(n,:) = [GP_i(n,1) GP_i(n,2:4) - GP_t(id,2:4)];
    end
end
GP_i_df(:,2:4) = GP_i_df(:,2:4) * 1e3;

m = 0;
for n = 1:size(GP_e),
    id = find (GP_t(:,1)==GP_e(n,1) );
    if length(id)==1
        m = m + 1;
        GP_e_df(n,:) = [GP_e(n,1) GP_e(n,2:4) - GP_t(id,2:4)];
    end
end
GP_e_df(:,2:4) = GP_e_df(:,2:4) * 1e3;

% Compute the statistics of initial EO difference
for n = 1:6,
    EO_i_stat(:,n) = comp_stat(EO_i_df(:,n+1));
end

% Compute the statistics of estimated EO difference
for n = 1:6,
    EO_e_stat(:,n)  = comp_stat(EO_e_df(:,n+1));
end

% Compute the statistics of initial GP difference
for n = 1:3,
    GP_i_stat(:,n)  = comp_stat(GP_i_df(:,n+1));
end

% Compute the statistics of estimated GP difference
for n = 1:3,
    GP_e_stat(:,n) = comp_stat(GP_e_df(:,n+1));
end

% Store the evaluation results
fid = fopen( strcat(path, 'Eva_summary.txt'), 'w');
fprintf ( fid, 'EO_i - Before Estimation\r\n' );
fprintf ( fid, '        %11s\t%11s\t%11s\t%11s\t%11s\t%11s\r\n', 'Xc[mm]', 'Yc[mm]', 'Zc[mm]', 'Omega[deg]', 'Phi[deg]', 'Kappa[deg]' );
fprintf ( fid, 'Minimum:%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_i_stat(1,:) );
fprintf ( fid, 'Maximum:%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_i_stat(2,:) );
fprintf ( fid, 'Average:%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_i_stat(3,:) );
fprintf ( fid, 'Std_dv.:%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_i_stat(4,:) );
fprintf ( fid, 'RMS    :%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_i_stat(5,:) );
fprintf ( fid, '\r\nEO_e - After Estimation\r\n' );
fprintf ( fid, '        %11s\t%11s\t%11s\t%11s\t%11s\t%11s\r\n', 'Xc[mm]', 'Yc[mm]', 'Zc[mm]', 'Omega[deg]', 'Phi[deg]', 'Kappa[deg]' );
fprintf ( fid, 'Minimum:%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_e_stat(1,:) );
fprintf ( fid, 'Maximum:%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_e_stat(2,:) );
fprintf ( fid, 'Average:%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_e_stat(3,:) );
fprintf ( fid, 'Std_dv.:%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_e_stat(4,:) );
fprintf ( fid, 'RMS    :%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_e_stat(5,:) );
fprintf ( fid, '\r\nGP_i - Before Estimation\r\n' );
fprintf ( fid, '        %11s\t%11s\t%11s\r\n', 'X[mm]', 'Y[mm]', 'Z[mm]' );
fprintf ( fid, 'Minimum:%11.3f\t%11.3f\t%11.3f\r\n', GP_i_stat(1,:) );
fprintf ( fid, 'Maximum:%11.3f\t%11.3f\t%11.3f\r\n', GP_i_stat(2,:) );
fprintf ( fid, 'Average:%11.3f\t%11.3f\t%11.3f\r\n', GP_i_stat(3,:) );
fprintf ( fid, 'Std_dv.:%11.3f\t%11.3f\t%11.3f\r\n', GP_i_stat(4,:) );
fprintf ( fid, 'RMS    :%11.3f\t%11.3f\t%11.3f\r\n', GP_i_stat(5,:) );
fprintf ( fid, '\r\nGP_e - After Estimation\r\n' );
fprintf ( fid, '        %11s\t%11s\t%11s\r\n', 'X[mm]', 'Y[mm]', 'Z[mm]' );
fprintf ( fid, 'Minimum:%11.3f\t%11.3f\t%11.3f\r\n', GP_e_stat(1,:) );
fprintf ( fid, 'Maximum:%11.3f\t%11.3f\t%11.3f\r\n', GP_e_stat(2,:) );
fprintf ( fid, 'Average:%11.3f\t%11.3f\t%11.3f\r\n', GP_e_stat(3,:) );
fprintf ( fid, 'Std_dv.:%11.3f\t%11.3f\t%11.3f\r\n', GP_e_stat(4,:) );
fprintf ( fid, 'RMS    :%11.3f\t%11.3f\t%11.3f\r\n', GP_e_stat(5,:) );
fclose(fid);

hist_x = -4950:50:4950;

% Visualize the EO difference before adjustment
figure
subplot(2,2,1);
axis on
bar(EO_i_df(:,1), EO_i_df(:,2), 0.5);
ylabel('\Delta x [mm]');
xlabel('ID');
title('EO_i: \Delta x');

subplot(2,2,2);
axis on
bar(EO_i_df(:,1), EO_i_df(:,2), 0.5);
ylabel('\Delta y [mm]');
xlabel('ID');
title('EO_i: \Delta y');

subplot(2,2,3);
axis on
bar(EO_i_df(:,1), EO_i_df(:,2), 0.5);
ylabel('\Delta z [mm]');
xlabel('ID');
title('EO_i: \Delta z');

% Visualize the EO difference before adjustment
figure
subplot(2,2,1);
axis on
bar(EO_i_df(:,1), EO_i_df(:,5), 0.5);
ylabel('\Delta \omega [deg]');
xlabel('ID');
title('EO_i: \Delta \omega');

subplot(2,2,2);
axis on
bar(EO_i_df(:,1), EO_i_df(:,6), 0.5);
ylabel('\Delta \phi [deg]');
xlabel('ID');
title('EO_i: \Delta \phi');

subplot(2,2,3);
axis on
bar(EO_i_df(:,1), EO_i_df(:,7), 0.5);
ylabel('\Delta \delta [deg]');
xlabel('ID');
title('EO_i: \Delta \kappa');

% Visualize the EO difference after adjustment
figure
subplot(2,2,1);
axis on
bar(EO_e_df(:,1), EO_e_df(:,2), 0.5);
ylabel('\Delta x [mm]');
xlabel('ID');
title('EO_e: \Delta x');

subplot(2,2,2);
axis on
bar(EO_e_df(:,1), EO_e_df(:,2), 0.5);
ylabel('\Delta y [mm]');
xlabel('ID');
title('EO_e: \Delta y');

subplot(2,2,3);
axis on
bar(EO_e_df(:,1), EO_e_df(:,2), 0.5);
ylabel('\Delta z [mm]');
xlabel('ID');
title('EO_e: \Delta z');

% Visualize the EO difference after adjustment
figure
subplot(2,2,1);
axis on
bar(EO_e_df(:,1), EO_e_df(:,5), 0.5);
ylabel('\Delta \omega [deg]');
xlabel('ID');
title('EO_e: \Delta \omega');

subplot(2,2,2);
axis on
bar(EO_e_df(:,1), EO_e_df(:,6), 0.5);
ylabel('\Delta \phi [deg]');
xlabel('ID');
title('EO_e: \Delta \phi');

subplot(2,2,3);
axis on
bar(EO_e_df(:,1), EO_e_df(:,7), 0.5);
ylabel('\Delta \delta [deg]');
xlabel('ID');
title('EO_e: \Delta \kappa');

% Visualize the ground point difference before adjustment
figure
subplot(2,2,1);
axis on
bar(GP_i_df(:,1), GP_i_df(:,2), 0.5);
ylabel('\Delta x [mm]');
xlabel('ID');
title('GP_e: \Delta x');

subplot(2,2,2);
axis on
bar(GP_i_df(:,1), GP_i_df(:,3), 0.5);
ylabel('\Delta y [mm]');
xlabel('ID');
title('GP_e: \Delta y');

subplot(2,2,3);
axis on
bar(GP_i_df(:,1), GP_i_df(:,4), 0.5);
ylabel('\Delta z [mm]');
xlabel('ID');
title('GP_e: \Delta z');

% Visualize the histogram of the ground point difference  before adjustment
figure
subplot(2,2,1);
axis on
hist(GP_i_df(:,2), hist_x);
xlabel('\Delta x [mm]');
title('GP_e: Histogram \Delta x');

subplot(2,2,2);
axis on
hist(GP_i_df(:,3),  hist_x);
xlabel('\Delta y [mm]');
title('GP_e: Histogram \Delta y');

subplot(2,2,3);
axis on
hist(GP_i_df(:,4),  hist_x);
xlabel('\Delta z [mm]');
title('GP_e: Histogram \Delta z');

subplot(2,2,4);
axis on
hist([GP_i_df(:,2); GP_i_df(:,3); GP_i_df(:,4)], hist_x);
xlabel('\Delta a [mm]');
title('Histogram \Delta a');

% Visualize the ground point difference after adjustment
figure
subplot(2,2,1);
axis on
bar(GP_e_df(:,1), GP_e_df(:,2), 0.5);
ylabel('\Delta x [mm]');
xlabel('ID');
title('GP_e: \Delta x');

subplot(2,2,2);
axis on
bar(GP_e_df(:,1), GP_e_df(:,3), 0.5);
ylabel('\Delta y [mm]');
xlabel('ID');
title('GP_e: \Delta y');

subplot(2,2,3);
axis on
bar(GP_e_df(:,1), GP_e_df(:,4), 0.5);
ylabel('\Delta z [mm]');
xlabel('ID');
title('GP_e: \Delta z');

% Visualize the histogram of the ground point difference after adjustment
figure
subplot(2,2,1);
axis on
hist(GP_e_df(:,2), hist_x);
xlabel('\Delta x [mm]');
title('GP_e: Histogram \Delta x');

subplot(2,2,2);
axis on
hist(GP_e_df(:,3),  hist_x);
xlabel('\Delta y [mm]');
title('GP_e: Histogram \Delta y');

subplot(2,2,3);
axis on
hist(GP_e_df(:,4),  hist_x);
xlabel('\Delta z [mm]');
title('GP_e: Histogram \Delta z');

subplot(2,2,4);
axis on
hist([GP_e_df(:,2); GP_e_df(:,3); GP_e_df(:,4)],  hist_x);
xlabel('\Delta a [mm]');
title('Histogram \Delta a');

% Visualize the EO
figure
plot3(EO_t(:,2), EO_t(:,3), EO_t(:,4), 'r-', EO_i(:,2), EO_i(:,3), EO_i(:,4), 'g-', ...
    EO_e(:,2), EO_e(:,3), EO_e(:,4), 'b-' );
view(2)
axis on
axis equal
grid on
xlabel('x[m]'), ylabel('y[m]'), zlabel('z[m]');



