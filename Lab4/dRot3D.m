% 3D rotation matrix
% Impyeong Lee
% University of Seoul
% 2003. 11. 20

function [dR_om, dR_ph, dR_kp] = dRot3D ( om, ph, kp )

Rx = [1 0 0; 0 cos(om) -sin(om); 0 sin(om) cos(om)];
Ry = [cos(ph) 0 sin(ph); 0 1 0; -sin(ph) 0 cos(ph)];
Rz = [cos(kp) -sin(kp) 0; sin(kp) cos(kp) 0; 0 0 1];

dRx = [0 0 0; 0 -sin(om) -cos(om); 0 cos(om) -sin(om)];
dRy = [-sin(ph) 0 cos(ph); 0 0 0; -cos(ph) 0 -sin(ph)];
dRz = [-sin(kp) -cos(kp) 0; cos(kp) -sin(kp) 0; 0 0 0];

dR_om = dRx * Ry * Rz;
dR_ph = Rx * dRy * Rz;
dR_kp = Rx * Ry * dRz;