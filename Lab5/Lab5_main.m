% Lab 5 'Relative Orientation'
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 1. Initialization of a pair of images, Converting CS into CCS
% 2. Computing object points' coordinates corresponding to a pair of image points
%    Input: a list of tie points, IO(cx, cy, focal length), ini_EO
%    Output: object points' coordinates 
% 3. RO
%    Input: a list of tie points, IO(cx,cy, focal length), ini_EO
% 4. Analysis of statistical results
% 5. Visualization of everything that we need for showing results (included
% in prob 2, 3, 4.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all %clear workspace
close all %close figures
clc %clear command line
format long g

%% 1. Initialization of a pair of images, Converting CS into CCS
% Before you run this code, compute MP in ics into in ccs
% Use Measure_TP.m to save a TP_ics.txt which includes a list of tie points
% between a pair of images
% Next, use Convert_CCS.m to save a TP_ccs_rc.txt to get MPs in ccs which
% are corrected parameters

% Read tie points in ccs
fid = fopen( 'TP_ccs.txt','r');
fgets(fid); fgets(fid); %ignore
TP = fscanf(fid, '%f', [5 inf])'; % inf 행 5개 열을 읽는것. (전치해서)
fclose(fid);
no_TP = size(TP,1); % # of tie points

% IO Parameters
resolution=[4000 3000];
p_size =  0.00153452;
IO = [0 0 4.73 p_size]'; % unit : mm

% Initialize Parameters
EO_i = [0 0 0 0 0 0; 0 5 0 0 0 0]; % initial EO
EO{1} = EO_i;
BP(:,1) = [1 4 2 3]';
[BP(1,2:3), BP(2,2:3), BP(3, 2:3), BP(4, 2:3)]=Cal_IB_CCS(IO, resolution); % Image Boundary in CCS

%% 2. Computing GP(MP)
GP_i = Comp_GP(TP, IO, EO_i); %nonchangable: TP, IO, changable: EO_i
GP{1} = GP_i;

%% Visualization of Initial State
% 3D Visualization
figure
subplot(1,2,1)
hold on
Comp_GP_with_Visualization(TP, IO, EO_i, BP);
view(3)
axis equal
grid on
xlabel('X_G'), ylabel('Y_G'), zlabel('Z_G')
hold off % default상태로 돌려줌
title('Initial');
ax1 = axis; %axis에 대한 정보를 저장해서 다음 그림에 사용해줄 수 있도록 함

% 2D Visualization 
subplot(1,2,2)
hold on
Comp_GP_with_Visualization(TP, IO, EO_i, BP);
view(2)
axis equal
grid on
xlabel('X_G'), ylabel('Y_G'), zlabel('Z_G')
hold off
title('Initial');
ax2 = axis;

%% 3. Relative Orientation (plus 5)
[EO, GP, stat] = RO_F(TP, GP, IO, EO);

for k = 1:size(stat{1},2)
    if mod(k, 2) ==0
        figure
        subplot(1,2,1)
        hold on
        Comp_GP_with_Visualization(TP, IO, EO{k+1}, BP);
        view(3)
        axis equal
        axis(ax1)
        grid on
        xlabel('X_G'), ylabel('Y_G'), zlabel('Z_G')
%         ax = gca;
%         ax.Clipping = 'off';
        hold off
        title(sprintf('After %d iterations (vct: %g)', k, stat{3}(k)));

        subplot(1,2,2)
        hold on
        Comp_GP_with_Visualization(TP, IO, EO{k+1}, BP);
        view(2)
        axis equal
        axis(ax2)
        grid on
        xlabel('X_G'), ylabel('Y_G'), zlabel('Z_G')
%         ax = gca;
%         ax.Clipping = 'off';
        hold off
        title(sprintf('After %d iterations (vct: %g)', k, stat{3}(k)));
    end
end

%% 4. Statistical Analyzation
% kt; ksihat
% et; residuals
% vct; variance component estimate
% Dkt; cofactor estimate
% stat = {kt, et, vct, Dkt};
no_it = size(stat{1},2);

% plot kt
% variances of EO1, EO2, MPs
figure ()
for n=1:no_it
    hold on
    % EO1&EO2
    y1=norm(stat{1}{n}(1:2)); % EO2 X^2+Z^2
    y2=norm(stat{1}{n}(3:5)); % EO2 norm(opk)
    
    % MP
    y3=norm(stat{1}{n}(6:23));
    plot(n, y1, 'r*', n, y2, 'b*', n, y3, 'go');
end
legend('Location of EO2','Orientation of EO2', 'Location of Model Points')
title('KsiHat', 'FontSize', 14);

% plot et
figure()
for k=1:no_it
    hold on
    plot(k, norm(stat{2}{k}), 'r-*', 'LineWidth', 2);
end
title('Norm of Residuals', 'FontSize', 14);

% plot vct
figure()
hold on
plot(1:1:no_it, stat{3}(:), 'b-o', 'LineWidth', 2);
title('Variance Component Estimate', 'FontSize', 14);

% plot correlation
% coef=corrcoef(veo); lab 4
figure()
coef=corrcoef(stat{4}{no_it}); % iteration 마지막 Dkt 분석
n = size(coef, 1);
imagesc(coef);
set(gca, 'XTick', 1:n); % center x-axis ticks on bins
set(gca, 'YTick', 1:n); % center y-axis ticks on bins
% set(gca, 'XTickLabel', L); % set x-axis labels
% set(gca, 'YTickLabel', L); % set y-axis labels
title('Correlation Matrix', 'FontSize', 14); % set title
colormap('gray'); % set the colorscheme
colorbar
