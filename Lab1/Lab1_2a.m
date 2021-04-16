clear;

fid = fopen( 'mavic2_eoini_tm.txt', 'r');
% % ¬˜∑ ¥Î∑Œ øµªÛ¿Ã∏ß Easting, Northing, Altitude, Yaw, Pitch, Roll
imgexif = textscan(fid, '%s %f %f %f %f %f %f', 'CommentStyle','#');
fclose(fid);

no_data = size(imgexif{2},1);
angles = zeros(no_data, 3);
for n = 1:no_data
   angles(n,:) = [imgexif{5}(n), imgexif{6}(n), imgexif{7}(n)];
end

% angle to rotation matrix
% ypr => pr(-y) 
% R = Rz * Rx * Ry;
A2R_1 = cell(no_data, 1);
for n = 1:no_data
    A2R_1{n} = A2R_RPY(angles(n,:));
end

% ypr => ypr
% R = Rx * Ry * Rz;
A2R_2 = cell(no_data, 1);
for n = 1:no_data
    A2R_2{n} = A2R_OPK1(angles(n,:));
end

% rotation matrix to angle Í≥ÑÏÇ∞
% ?ñâ?†¨?? ?ã§?ùåÍ≥? Í∞ôÏù¥ ?†ï?ùò
% R = Rz * Rx * Ry;
R2A_1 = zeros(no_data, 3);
for n = 1:no_data
    R2A_1_tmp = R2A_RPY(A2R_1{n});
    R2A_1_a = R2A_1_tmp(1,:);
    R2A_1_b = R2A_1_tmp(2,:);
end

% rotation matrix to angle 
% R = Rx * Ry * Rz;
R2A_2 = zeros(no_data, 3);
for n = 1:no_data
    R2A_2_tmp = R2A_OPK1(A2R_1{n});
    R2A_2_a = R2A_2_tmp(1,:);
    R2A_2_b = R2A_2_tmp(2,:);
end
