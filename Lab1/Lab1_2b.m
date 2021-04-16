clear;
close all;

%% ���� �б�
fid = fopen( 'mavic2_eoini_tm.txt', 'r');
% % ���ʴ�� �����̸� Easting, Northing, Altitude, Yaw, Pitch, Roll
imgexif = textscan(fid, '%s %f %f %f %f %f %f', 'CommentStyle','#','Delimiter','\t');
fclose(fid);
no_data = size(imgexif{2},1);

%% ���� ��ǥ�迡 ǥ���� ī�޶� ��ġ/�ڼ� �Է�
C_BG = zeros(no_data,3);
A_BG = zeros(no_data,3);
R_BG = cell(no_data,1);

for n = 1:no_data
    C_BG(n,:) = [imgexif{2}(n) imgexif{3}(n) imgexif{4}(n)];   
    A_BG(n,:) = [imgexif{5}(n) imgexif{6}(n) imgexif{7}(n)];
    R_BG{n} = A2R_OPK1(A_BG(n,:));
    R_BG{n}(:,3) = -R_BG{n}(:,3);
end

%% ���� ��ǥ�迡 ǥ���� ���� ��ġ
%mavic pro spec
p_size = (0.00241228)/1; %���� : m
f= (6.3)/1; %���� : m
f_mat = [0 0 f];
s_width = 4000;
s_height = 3000;

u_l_ics = [0 0 0];
d_l_ics = [0 s_width 0];
u_r_ics = [s_height 0 0];
d_r_ics = [s_height s_width 0];

cam_origin=[s_height/2 s_width/2 0];

% ������ 3���� ������ǥ������ ��ġ
% ���� �ٿ���� 4���� ��ġ���� x,y,z�� 3���� ���� ����
% (up left)      (up right)
% |�ѤѤѤѤѤѤѤ�|
% |       X       |
% |�ѤѤѤѤѤѤѤ�|
% (down left)    (down right)

% ICS -> PCS
u_l_pcs = zeros(no_data, 3);
d_l_pcs = zeros(no_data, 3);
u_r_pcs = zeros(no_data, 3);
d_r_pcs = zeros(no_data, 3);

for n=1:no_data
    u_l_pcs(n,:) = (u_l_ics-cam_origin)*p_size-f_mat;
    d_l_pcs(n,:) = (d_l_ics-cam_origin)*p_size-f_mat;
    u_r_pcs(n,:) = (u_r_ics-cam_origin)*p_size-f_mat;
    d_r_pcs(n,:) = (d_r_ics-cam_origin)*p_size-f_mat;
 
    u_l_pcs(n,2) = -u_l_pcs(n,2);
    d_l_pcs(n,2) = -d_l_pcs(n,2);
    u_r_pcs(n,2) = -u_r_pcs(n,2);
    d_r_pcs(n,2) = -d_r_pcs(n,2);
end

%PCS -> GCS
u_l_gcs = zeros(no_data, 3);
d_l_gcs = zeros(no_data, 3);
u_r_gcs = zeros(no_data, 3);
d_r_gcs = zeros(no_data, 3);

for n=1:no_data
    u_l_gcs(n,:) = (R_BG{n}'*u_l_pcs(n,:)'+C_BG(n,:)')';
    d_l_gcs(n,:) = (R_BG{n}'*d_l_pcs(n,:)'+C_BG(n,:)')';
    u_r_gcs(n,:) = (R_BG{n}'*u_r_pcs(n,:)'+C_BG(n,:)')';
    d_r_gcs(n,:) = (R_BG{n}'*d_r_pcs(n,:)'+C_BG(n,:)')';
end

pt=zeros(no_data, 3);
cs = {'r', 'g', 'b'};
img_x=cell(no_data, 1);
img_y=cell(no_data, 1);
img_z=cell(no_data, 1);
figure ()
hold on
for n = 1:no_data
    %% ī�޶� ��ǥ�� �׸���
    for k = 1:3
        cstr = cs{k};
        pt(n,:) = C_BG(n,:)' + 3 * R_BG{n}(k,:)';
        plot3([C_BG(n,1), pt(n,1)], [C_BG(n,2), pt(n,2)], [C_BG(n,3), pt(n,3)], 'r-', 'LineWidth', 1, 'Color', cstr);
    end
    %% ���� ��ġ �׸���
    img_x{n} = [u_l_gcs(n,1) u_r_gcs(n,1) d_r_gcs(n,1) d_l_gcs(n,1) u_l_gcs(n,1)]';
    img_y{n} = [u_l_gcs(n,2) u_r_gcs(n,2) d_r_gcs(n,2) d_l_gcs(n,2) u_l_gcs(n,2)]';
    img_z{n} = [u_l_gcs(n,3) u_r_gcs(n,3) d_r_gcs(n,3) d_l_gcs(n,3) u_l_gcs(n,3)]';    
    plot3(img_x{n}, img_y{n}, img_z{n}, 'b-', 'LineWidth', 2);
end
view(3)
grid on, axis equal
title('R(BG)') 
xlabel('X'), ylabel('Y'), zlabel('Z')