% Lab 6 'Estimate Homography and Fundamental matrix'
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 1. Feature matching with SURF
% 2. Estimate Homography with a pair of images
% 3. Accuracy of H
% 4. Visualize Homography test
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clc
clear all;
close all;

format long g;
%% Feature matching using SURF.m
I1 = imread('rubberwhale1.png');
I2 = imread('rubberwhale2.png');
[img1, img2] = SURF(rgb2gray(I1), rgb2gray(I2), 20, 1);
img1pts = img1.Location;
img2pts = img2.Location;
no_TP = size(img1pts, 1);
[imgwidth, imgheight, dim] = size(I1);

%% Estimate my Homography matrix
% img1pts : img points in image 1
% img2pts : img points in image 2

[A, myH] = estimateH(img1pts, img2pts); % from my code
addpath('ransac_homography')  
[hRANSAC, corr] = findHomography(img1pts', img2pts'); % from finding homography with RANSAC

%% Initialization for test
fid = fopen( 'rubberwhale_test.txt','r');
fgets(fid); fgets(fid); %ignore
HF_test = fscanf(fid, '%f', [5 inf])'; % inf
fclose(fid);
no_HFtest = size(HF_test,1); % # of tie points 
test_img1pts = [HF_test(:,4) HF_test(:,5) ones(no_HFtest,1)]';
test_img2pts = [HF_test(:,2) HF_test(:,3) ones(no_HFtest,1)]';

cal_img2pts_myH = zeros(no_HFtest, 3);
diff_myH = zeros(no_HFtest, 3);

cal_img2pts_hRANSAC = zeros(no_HFtest, 3);
diff_hRANSAC = zeros(no_HFtest, 3);

resi_F = zeros(no_HFtest,1);
%% Homography test
for n=1:no_HFtest
    % subtract truth values of the second image to validate Homography matrix
    cal_img2pts_myH(n,:) = ((myH/myH(end))*test_img1pts(:,n))';
    diff_myH(n,:) = cal_img2pts_myH(n,:)/cal_img2pts_myH(n,end) - test_img2pts(:,n)';
    
    % subtract truth values by yourH (using RANSAC)
    cal_img2pts_hRANSAC(n,:) = (hRANSAC*test_img1pts(:,n))';
    diff_hRANSAC(n,:) = cal_img2pts_hRANSAC(n,:)/cal_img2pts_hRANSAC(n,end) - test_img2pts(:,n)';
end

%% Visualize truth values and calculated points by Homography
figure(3)
imshow(I2); hold on
plot(cal_img2pts_myH(:,1), cal_img2pts_myH(:,2), 'b+', 'MarkerSize', 8);
hold on;
plot(test_img2pts(1,:), test_img2pts(2,:), 'ro', 'MarkerSize', 8);
legend('calculated points by homography','truth values')
hold off;

figure(4)
imshow(I2); hold on
plot(cal_img2pts_hRANSAC(:,1), cal_img2pts_hRANSAC(:,2), 'b+', 'MarkerSize', 8);
hold on;
plot(test_img2pts(1,:), test_img2pts(2,:), 'ro', 'MarkerSize', 8);
legend('calculated points by homography','truth values')

