% Determine the initial approximations for AT
% Impyeong Lee
% 2009. 1. 12

clear all   % clear all the variable previsoully defined
close all   % close all the opend figure window
clc         % clear the command prompt window


path = 'C1\\';

% Read the true IO
fid = fopen( strcat(path,'IO_t.txt'),'r'); %% true values of IO
IO = fscanf(fid, '%f');
fclose(fid);

% Read the image points
fid = fopen( strcat(path,'IP_m.txt'),'r'); %% values of image points including random errors
IP = fscanf(fid, '%f', [4 inf] )'; %% 4 columns
fclose(fid);
no_IP = size(IP,1);

% Read true values for EO
fid = fopen( strcat(path,'EO_t.txt'),'r'); %% true values of EO
EO_t = fscanf(fid, '%f', [7 inf] )'; %% X, Y, Z, o , p, k
fclose(fid);

% Read directly measured EO
fid = fopen( strcat(path,'EO_m.txt'),'r'); %% values of EO including random errors
EO_m = fscanf(fid, '%f', [7 inf] )';
fclose(fid);

% Read true values for GP
fid = fopen( strcat(path,'GP_t.txt'),'r'); %% true values of ground points
GP_t = fscanf(fid, '%f', [4 inf] )'; % x, y, z
fclose(fid);

%% Intermediate Outputs
% Compute the no. images
id_IM = [];
for n = 1:size(IP,1)
    if length ( find( id_IM == IP(n,1) ) ) == 0,
        id_IM = [id_IM; IP(n,1)];
    end
end
no_IM = size(id_IM,1);

% Compute the no. ground points
id_GP = [];
for n = 1:size(IP,1)
    if length ( find( id_GP == IP(n,2) ) ) == 0,
        id_GP = [id_GP; IP(n,2)];
    end
end
id_GP = sort(id_GP);
no_GP = size(id_GP, 1);

% Compute no. images of each GP appearing
cnt_GP = zeros( no_GP, 1);
for n = 1:no_IP,
    idx = find ( id_GP == IP(n,2) );
    cnt_GP(idx) = cnt_GP(idx) + 1; %% 
end
id_GP_inc = id_GP ( find( cnt_GP > 1 ) ); %% find id of ground point which appears more than 1 time per image
id_GP_inc = sort(id_GP_inc);
id_GP_exc = id_GP ( find( cnt_GP <= 1 ) );
no_GP_inc = length(id_GP_inc);

% Select GP and IP for AT
m = 0;
for np = 1:no_IP,
    if length(find(id_GP_exc==IP(np,2))) == 0,
        % check if the point index is in the index of gp exclusive
        m = m + 1;
        IP_inc(m,:) = IP(np,:);
    end
end
no_IP_inc = size(IP_inc,1);

IP_inc = sortrows(IP_inc,[2 1]);

%% Initial Approximations of EOs
% Start a timer
tic;

% Determine the initial approximation of EOs
EO_i = EO_m;

for ni = 1:no_IM,
    R_e{ni} = Rot3D(EO_i(ni,5:7));
end
    
for gi = 1:no_GP_inc,
    idx = find( IP_inc(:,2) == id_GP_inc(gi) );
    no_line = length(idx);
    no_line_pair = no_line * ( no_line - 1 ) / 2;
    A = zeros( no_line_pair * 3, no_line );
    y = zeros( no_line_pair * 3, 1 );
    cnt = 0;
    for n = 1:length(idx)-1,
        I1 = IP_inc(idx(n),1);
        IP1 = [ IP_inc(idx(n),3:4)' - IO(1:2); -IO(3)];
        for m = n+1:length(idx),
            cnt = cnt + 1;
            I2 = IP_inc(idx(m),1);
            IP2 = [ IP_inc(idx(m),3:4)' - IO(1:2); -IO(3)];
            A(3*cnt-2:3*cnt, n) =  R_e{I1}' * IP1;
            A(3*cnt-2:3*cnt, m) =  -R_e{I2}' * IP2;
            y(3*cnt-2:3*cnt, 1) = EO_i(I2,2:4)' - EO_i(I1,2:4)';            
        end
    end
    kt = inv(A'*A)*A'*y;
    et = y - A * kt;
    Om = et' * et;
    vct = Om / ( size(A,1) - rank(A) );
    et_stat = comp_stat(et);

    GP = [];
    for n = 1:length(idx)
        IID = IP_inc(idx(n),1);
        IP = [ IP_inc(idx(n),3:4)' - IO(1:2); -IO(3)];
        GP(:,n) = kt(n) * R_e{IID}' * IP + EO_i(IID,2:4)';
    end
    GP_i(gi,:) = [id_GP_inc(gi) mean(GP,2)'];
    Est_stat(gi,:) = [id_GP_inc(gi) vct et_stat'];
end

Time_Elapsed = toc

% Compute the statistics of measured EO difference
for n = 1:6,
    EO_i_stat(:,n) = comp_stat(EO_i(:,n+1)-EO_t(:,n+1));
end
EO_i_stat(:,1:3) = EO_i_stat(:,1:3) * 1e3;
EO_i_stat(:,4:6) = EO_i_stat(:,4:6) * 180 / pi;

% Compute the statistics of measure GP difference
for n = 1:3,
    GP_i_stat(:,n)  = comp_stat(GP_i(:,n+1)-GP_t(:,n+1));
end
GP_i_stat = GP_i_stat * 1e3;

%% Store Outputs
% Store the initial appox. for EO
fid = fopen(strcat(path, 'EO_i.txt'), 'w');
fprintf(fid, '%d\t%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_i' );
fclose(fid);

% Store the initial appox. for GP
fid = fopen(strcat(path, 'GP_i.txt'), 'w');
fprintf(fid, '%d\t%11.3f\t%11.3f\t%11.3f\r\n', GP_i' );
fclose(fid);

% Store Initial Approximation Summary
fid = fopen(strcat(path, 'IA_summary.txt'), 'w');
fprintf(fid, 'Processing time: %g\r\n', Time_Elapsed );
fprintf ( fid, 'Statistics for EO_i - EO_t\r\n' );
fprintf ( fid, '        %11s\t%11s\t%11s\t%11s\t%11s\t%11s\r\n', 'Xc[mm]', 'Yc[mm]', 'Zc[mm]', 'Omega[deg]', 'Phi[deg]', 'Kappa[deg]' );
fprintf ( fid, 'Minimum:%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_i_stat(1,:) );
fprintf ( fid, 'Maximum:%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_i_stat(2,:) );
fprintf ( fid, 'Average:%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_i_stat(3,:) );
fprintf ( fid, 'Std_dv.:%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_i_stat(4,:) );
fprintf ( fid, 'RMS    :%11.3f\t%11.3f\t%11.3f\t%11.6f\t%11.6f\t%11.6f\r\n', EO_i_stat(5,:) );
fprintf ( fid, 'Statistics for GP_i - GP_t\r\n' );
fprintf ( fid, '        %11s\t%11s\t%11s\r\n', 'X[mm]', 'Y[mm]', 'Z[mm]' );
fprintf ( fid, 'Minimum:%11.3f\t%11.3f\t%11.3f\r\n', GP_i_stat(1,:) );
fprintf ( fid, 'Maximum:%11.3f\t%11.3f\t%11.3f\r\n', GP_i_stat(2,:) );
fprintf ( fid, 'Average:%11.3f\t%11.3f\t%11.3f\r\n', GP_i_stat(3,:) );
fprintf ( fid, 'Std_dv.:%11.3f\t%11.3f\t%11.3f\r\n', GP_i_stat(4,:) );
fprintf ( fid, 'RMS    :%11.3f\t%11.3f\t%11.3f\r\n', GP_i_stat(5,:) );
fprintf(fid, 'Statistics for GP initial appoximation\r\n' );
fprintf(fid, 'GP\t%11s\t%11s\t%11s\t%11s\t%11s\t%11s\r\n', 'vct', 'Min', 'Max', 'Avg', 'Std_dv', 'RMS' );
fprintf(fid, '%d\t%11.3f\t%11.3f\t%11.3f\t%11.3f\t%11.3f\t%11.3f\r\n', Est_stat' );
fclose(fid);

% Visualize the variance component estimate
figure
axis on
bar(Est_stat(:,1), Est_stat(:,2), 0.5);
ylabel('Variance Component Estimate [m]');
xlabel('ID');
title('VCT for GP Init Apporx');
