% Stage 6. Visualization: calculating image frame
% Ilseo Jeon
% Dept. of Geoinformatics, The University of Seoul
% Adv. Photogrammetry Lab4 in Master Course
% 2020. 05. 8
%  ul             ur
%   |天天天天天天天|
%   |      X      |
%   |天天天天天天天|
%   dl            dr
% gp; ground point
% pcs; camera coordinate system
% ics; image coordinate system

function [gcs_ul, gcs_dl, gcs_dr, gcs_ur, R]=Cal_IF(GCP, EO, IO, resolution)
[no, no_it]=size(EO);
R_IC=[1 0 0; 0 -1 0; 0 0 -1];
ul_ics = [0 0 0]';
dl_ics = [0 resolution(1) 0]';
ur_ics = [resolution(2) 0 0]';
dr_ics = [resolution(2) resolution(1) 0]';

f_mat=[0 0 IO(3)*IO(4)]';
img_origin =[resolution(1)/2 resolution(2)/2 0]';

NED = [0 1 0; ...
       1 0 0; ...
       0 0 -1];
% avg_z = mean(GCP(3,:));
% scale=(mean(GCP(3,:))-EO(3))/(-IO(3));
% vis_coord_system(EO(1:3,no_it)'+5, NED', 5, sprintf("NED"),'r');

pcs_ul=R_IC*((ul_ics-img_origin)*IO(4)-f_mat);
pcs_dl=R_IC*((dl_ics-img_origin)*IO(4)-f_mat);
pcs_dr=R_IC*((dr_ics-img_origin)*IO(4)-f_mat);
pcs_ur=R_IC*((ur_ics-img_origin)*IO(4)-f_mat);

for n=1:no_it
    R{n}=NED*Rot3D(EO(4,n), EO(5,n), EO(6,n))*R_IC;

    gcs_ul(:,n) = R{n}(:,:)'*pcs_ul+EO(1:3,n);
    gcs_dl(:,n) = R{n}(:,:)'*pcs_dl+EO(1:3,n);
    gcs_dr(:,n) = R{n}(:,:)'*pcs_dr+EO(1:3,n);
    gcs_ur(:,n) = R{n}(:,:)'*pcs_ur+EO(1:3,n);
    
%     gp_ul(n,:) = scale*pcs_ul(n,:);
%     gp_ul(n,:)=(R*gp_ul(n,:)'+esti_EO(n,1:3)')';
%     gp_ul(n,3) = avg_z;
%     gp_dl(n,:) = scale*pcs_dl(n,:);
%     gp_dl(n,:)=(R*gp_dl(n,:)'+esti_EO(n,1:3)')';
%     gp_dl(n,3) = avg_z;
%     gp_ur(n,:) = scale*pcs_dr(n,:);
%     gp_ur(n,:)=(R*gp_ur(n,:)'+esti_EO(n,1:3)')';
%     gp_ur(n,3) = avg_z;
%     gp_dr(n,:) = scale*pcs_ur(n,:);
%     gp_dr(n,:)=(R*gp_dr(n,:)'+esti_EO(n,1:3)')';
%     gp_dr(n,3) = avg_z;
end

end