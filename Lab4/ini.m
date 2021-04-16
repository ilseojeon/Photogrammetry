% Stage 3. Calculatiing initial approximates of camera EOs
% Ilseo Jeon
% Dept. of Geoinformatics, The University of Seoul
% Adv. Photogrammetry Lab4 in Master Course
% 2020. 05. 7

function [ini_EO] = ini(GCP, IO, IP_ccs)

[no_data, no_var] = size(GCP);
Q = eye(no_data*2);

for n = 1:no_data
    %observatioons; in ccs
    L(2*n-1:2*n, 1) = [IP_ccs(n,1) IP_ccs(n,2)]';
    %design matrix 
    A(2*n-1:2*n, :) = [GCP(n,1) -GCP(n,2) 1 0; GCP(n,2) GCP(n,1) 0 1]; 
end
%cos(a) = -kt(1) sin(a) = -kt(2)
%xm=kt(3), ym=kt(4)
[kt, et, vct, Dkt] = LSE (A, L, Q);
avg_z = mean(GCP(:,3)); 
XYc = inv([-kt(1) kt(2); -kt(2) -kt(1)])*[kt(3); kt(4)];
K = atan2(kt(2), kt(1));
h = IO(3)*cos(K)/kt(1);
Zc = h + avg_z;
ini_EO = [XYc(1) XYc(2) Zc 0 0 K];
end