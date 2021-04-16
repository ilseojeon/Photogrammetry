%% estimate Fundamental matrix
% Ilseo Jeon
% 2020 Adv. Photogrammetry, Geoinformatics, University of Seoul
% calculate Homography from the same projection

function F = estimateF(img1pts, img2pts)
format long g
for n = 1:size(img1pts, 1)
    A(n,:) = kron(img1pts(n,:), img2pts(n,:));
end

[U, D, V] = svd(A);
eig_val = diag(D);
Fa = reshape(V(:, 9), 3, 3)';

[Ua, Da, Va] = svd(Fa);

F = Ua * diag([Da(1,1), Da(2,2), 0]) * Va';

end