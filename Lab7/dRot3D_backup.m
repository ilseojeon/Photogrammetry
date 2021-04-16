% 3D rotation matrix
% Impyeong Lee
% University of Seoul
% 2003. 11. 20
% 2004. 11. 28 - [om,ph,kp] -> ra

function [dR_om, dR_ph, dR_kp] = dRot3D ( ra );

om = ra(1);
ph = ra(2);
kp = ra(3);

dR_om(1,1) = 0;
dR_om(1,2) = 0;
dR_om(1,3) = 0;

dR_om(2,1) = -sin(om) * sin(kp) + cos(om) * sin(ph) * cos(kp);
dR_om(2,2) = -sin(om) * cos(kp) - cos(om) * sin(ph) * sin(kp);
dR_om(2,3) = -cos(om) * cos(ph);

dR_om(3,1) =  cos(om) * sin(kp) + sin(om) * sin(ph) * cos(kp);
dR_om(3,2) =  cos(om) * cos(kp) - sin(om) * sin(ph) * sin(kp);
dR_om(3,3) = -sin(om) * cos(ph);

dR_ph(1,1) = -sin(ph) * cos(kp);
dR_ph(1,2) =  sin(ph) * sin(kp);
dR_ph(1,3) =  cos(ph);

dR_ph(2,1) =  sin(om) * cos(ph) * cos(kp);
dR_ph(2,2) = -sin(om) * cos(ph) * sin(kp);
dR_ph(2,3) =  sin(om) * sin(ph);

dR_ph(3,1) = -cos(om) * cos(ph) * cos(kp);
dR_ph(3,2) =  cos(om) * cos(ph) * sin(kp);
dR_ph(3,3) = -cos(om) * sin(ph);

dR_kp(1,1) = -cos(ph) * sin(kp);
dR_kp(1,2) = -cos(ph) * cos(kp);
dR_kp(1,3) = 0;

dR_kp(2,1) =  cos(om) * cos(kp) - sin(om) * sin(ph) * sin(kp);
dR_kp(2,2) = -cos(om) * sin(kp) - sin(om) * sin(ph) * cos(kp);
dR_kp(2,3) = 0;

dR_kp(3,1) =  sin(om) * cos(kp) + cos(om) * sin(ph) * sin(kp);
dR_kp(3,2) = -sin(om) * sin(kp) + cos(om) * sin(ph) * cos(kp);
dR_kp(3,3) = 0;
