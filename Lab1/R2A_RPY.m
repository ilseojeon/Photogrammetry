% 3D rotation matrix
% Impyeong Lee
% University of Seoul
% 2003. 11. 20
% 2004. 11. 28 - [om,ph,kp] -> ra

function ra = R2A_RPY ( R )

% R = Rz * Rx * Ry;

s_om = R(3,2);
c_om = zeros(2,1);
c_om(1) = sqrt(1-s_om^2);
c_om(2) = -sqrt(1-s_om^2);

kp = zeros(2,1);
ph = zeros(2,1);
kp(1) = atan2(-R(1,2), R(2,2));
ph(1) = atan2(-R(3,1), R(3,3));
kp(2) = atan2(R(1,2), -R(2,2));
ph(2) = atan2(R(3,1), -R(3,3));

err = zeros(2,1);
for n = 1:2,
   err(n) = abs(cos(kp(n))*cos(ph(n)) - sin(kp(n))*s_om*sin(ph(n)) - R(1,1));
   err(n) = err(n) + abs(cos(kp(n))*sin(ph(n)) + sin(kp(n))*s_om*cos(ph(n)) - R(1,3));
   err(n) = err(n) + abs(sin(kp(n))*cos(ph(n)) + cos(kp(n))*s_om*sin(ph(n)) - R(2,1));
   err(n) = err(n) + abs(sin(kp(n))*sin(ph(n)) - cos(kp(n))*s_om*cos(ph(n)) - R(2,3));
end
[val,idx] = min(err);
if val > 1e-6,
    fprintf('Error!!');
end

% ra(1) = om(idx);
% ra(2) = atan2(s_ph, c_ph(idx));
% ra(3) = kp(idx);

ra(:,1) = ph;
ra(:,2) = atan2(s_om, c_om);
ra(:,3) = -kp;

