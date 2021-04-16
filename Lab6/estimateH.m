%% estimate Homography
% Ilseo Jeon
% 2020 Adv. Photogrammetry, Geoinformatics, University of Seoul
% calculate Homography from the same projection

function [A, H] = estimateH(img1pts, img2pts)
format long g
no_TP = size(img1pts,1);
A = zeros(no_TP, 9);
for i = 1:no_TP
    A(2*i-1,:) = [-img1pts(i,1) -img1pts(i,2) 1 0 0 0 img1pts(i,1)*img2pts(i,1) img2pts(i,1)*img1pts(i,2) img2pts(i,1)];
    A(2*i,:) = [0 0 0 -img1pts(i,1) -img1pts(i,2) 1 img1pts(i,1)*img2pts(i,2) img2pts(i,2)*img1pts(i,2) img2pts(i,2)];
end

[U, D, V] = svd(A);

eig_val=diag(D);
[min_eig_val, min_eig_val_idx] = min(eig_val);
H = reshape(V(:,min_eig_val_idx), 3, 3)';

end
