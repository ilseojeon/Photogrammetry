% 3D rotation matrix
% Impyeong Lee
% University of Seoul
% 2003. 11. 20
% 2004. 11. 28 - [om,ph,kp] -> ra

function ra = R2A_OPK1 ( R )

% R = Rx * Ry * Rz;

s_ph = R(1,3); % Kappa
c_ph = zeros(2,1);
c_ph(1) = sqrt(1-s_ph^2); % sqrt(1-kappa^2)
c_ph(2) = -sqrt(1-s_ph^2); % -sqrt(1-kappa^2)

kp = zeros(2,1); 
om = zeros(2,1);
kp(1) = atan2(-R(1,2), R(1,1)); % [ -
om(1) = atan2(-R(2,3), R(3,3));
kp(2) = atan2(R(1,2), -R(1,1));
om(2) = atan2(R(2,3), -R(3,3));

err = zeros(2,1);
for n = 1:2,
   err(n) = abs(cos(om(n))*sin(kp(n)) + sin(om(n))*s_ph*cos(kp(n)) - R(2,1));
   err(n) = err(n) + abs(cos(om(n))*cos(kp(n)) - sin(om(n))*s_ph*sin(kp(n)) - R(2,2));
   err(n) = err(n) + abs(sin(om(n))*sin(kp(n)) - cos(om(n))*s_ph*cos(kp(n)) - R(3,1));
   err(n) = err(n) + abs(sin(om(n))*cos(kp(n)) + cos(om(n))*s_ph*sin(kp(n)) - R(3,2));
end
[val,idx] = min(err);
if val > 1e-6,
    fprintf('Error!!');
end

% ra(1) = om(idx);
% ra(2) = atan2(s_ph, c_ph(idx));
% ra(3) = kp(idx);

ra(:,1) = om;
ra(:,2) = atan2(s_ph, c_ph);
ra(:,3) = kp;

