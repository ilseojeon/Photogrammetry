% 3D rotation matrix
% Impyeong Lee
% University of Seoul
% 2003. 11. 20
% 2004. 11. 28 - [om,ph,kp] -> ra

function R = Rot3D ( ra );

om = ra(1);
ph = ra(2);
kp = ra(3);

R(1,1) = cos(ph) * cos(kp);
R(1,2) = -cos(ph) * sin(kp);
R(1,3) = sin(ph);

R(2,1) = cos(om) * sin(kp) + sin(om) * sin(ph) * cos(kp);
R(2,2) = cos(om) * cos(kp) - sin(om) * sin(ph) * sin(kp);
R(2,3) = -sin(om) * cos(ph);

R(3,1) = sin(om) * sin(kp) - cos(om) * sin(ph) * cos(kp);
R(3,2) = sin(om) * cos(kp) + cos(om) * sin(ph) * sin(kp);
R(3,3) = cos(om) * cos(ph);
