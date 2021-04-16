% Calculating Image Boundary in CCS
% Ilseo Jeon
% Dept. of Geoinformatics, The University of Seoul
% Adv. Photogrammetry Lab4 in Master Course
% 2020. 05. 26
%  ul(3)          ur(2)
%   |天天天天天天天|
%   |      X      |
%   |天天天天天天天|
%   dl(1)         dr(4)
% ccs; camera coordinate system
% ics; image coordinate system

function [ccs_dl, ccs_dr, ccs_ur, ccs_ul]=Cal_IB_CCS(IO, resolution)

R_IC=[1 0; 0 -1];
ul_ics = [0 0]';
dl_ics = [0 resolution(2)]';
ur_ics = [resolution(1) 0]';
dr_ics = [resolution(1) resolution(2)]';

img_origin =[resolution(1)/2 resolution(2)/2]';

ccs_ul=(R_IC*(ul_ics-img_origin)*IO(4))';
ccs_dl=(R_IC*(dl_ics-img_origin)*IO(4))';
ccs_dr=(R_IC*(dr_ics-img_origin)*IO(4))';
ccs_ur=(R_IC*(ur_ics-img_origin)*IO(4))';

end