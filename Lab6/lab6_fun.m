% Lab 6 'Estimate Homography and Fundamental matrix'
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 1. Load .txt file
% 2. Estimate Fundamental matrix with a pair of images
% 3. Accuracy of F
% 4. Visualize epipolar lines from F matrix
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clc
clear all;
close all;

format long g;
%% load .txt file including tie points taken manually
%use this if you use .txt file including mathced points (in case you got
%mathced points taken manually

% fid = fopen( 'stereo-getty.txt','r');
% fgets(fid); fgets(fid); %ignore
% TP = fscanf(fid, '%f', [5 inf])'; % inf 
% fclose(fid);
% no_TP = size(TP,1); % # of tie points
% 
% I1 = imread('stereo-getty-l.jpg'); % left image
% figure(1)
% imshow(I1); hold on
% plot(TP(:,4), TP(:,5), 'r+','MarkerSize', 8);
% I2 = imread('stereo-getty-r.jpg'); % right image
% figure(2)
% imshow(I2); hold on
% plot(TP(:,2), TP(:,3), 'r+','MarkerSize', 8);
% 
% [imgwidth, imgheight, dim] = size(I1);
% img1pts = [TP(:,2) TP(:,3)]; % u,v
% img2pts = [TP(:,4) TP(:,5)]; % u,v

%% Feature matching using SURF.m
I1 = imread('stereo-corridor_l.gif');
I2 = imread('stereo-corridor_r.gif');
[img1, img2] = SURF(I1, I2, 150, 1);
img1pts = img1.Location;
img2pts = img2.Location;
no_TP = size(img1pts, 1);
[imgwidth, imgheight, dim] = size(I1);

%% Estimate my Fundamental matrix
img1pts(:,3)=eye;
img2pts(:,3)=eye;

myF = estimateF(img1pts, img2pts); % ilseo's  
[fRANSAC, inliers] = estimateFundamentalMatrix(img1pts(:,1:2),...
    img2pts(:,1:2),'Method','RANSAC',...
    'NumTrials',3000,'DistanceThreshold',1e-4); % from finding fundamental matrix with RANSAC

%% Initialization for test
fid = fopen( 'stereo-corridor_test.txt','r');
fgets(fid); fgets(fid); %ignore
HF_test = fscanf(fid, '%f', [5 inf])'; % inf
fclose(fid);
no_HFtest = size(HF_test,1); % # of tie points 
test_img1pts = [HF_test(:,2) HF_test(:,3) ones(no_HFtest,1)]';
test_img2pts = [HF_test(:,4) HF_test(:,5) ones(no_HFtest,1)]';

%% Fundamental matrix test
% use symmetric epipolar distance formula see Ch11.5 in Mutiple View Geometry
for n=1:no_HFtest
    numerator1(n,:) = power(test_img2pts(:,n)'*(myF/myF(end))*test_img1pts(:,n),2);
    epl1(n,:) = ((myF/myF(end))*test_img2pts(:,n))';
    epl1(n,:) = epl1(n,:)/epl1(n,3);
    epl2(n,:) = ((myF/myF(end))'*test_img1pts(:,n))';
    epl2(n,:) = epl2(n,:)/epl2(n,3); 
    d_1(n,:) = power(epl1(n,1),2)+power(epl1(n,2),2);
    d_2(n,:) = power(epl2(n,1),2)+power(epl2(n,2),2);
    resi_myF(n,:) = numerator1(n,:)*(1/d_1(n,:) + 1/d_2(n,:));
end
checkmyF = mean(resi_myF);

for n=1:no_HFtest
    numerator2(n) = power(test_img2pts(:,n)'*(fRANSAC/fRANSAC(end))*test_img1pts(:,n),2);
    epl3(n,:) = ((fRANSAC/fRANSAC(end))*test_img2pts(:,n))';
    epl3(n,:) = epl3(n,:)/epl3(n,3);
    epl4(n,:) = ((fRANSAC/fRANSAC(end))'*test_img1pts(:,n))';
    epl4(n,:) = epl4(n,:)/epl4(n,3); 
    d_3(n) = power(epl3(n,1),2)+power(epl3(n,2),2);
    d_4(n) = power(epl4(n,1),2)+power(epl4(n,2),2);
    resi_fRANSAC(n) = numerator2(n)*(1/d_3(n) + 1/d_4(n));
end
checkfRANSAC = mean(resi_fRANSAC);

%% Visualize epipolarlines
% by epilines built-in function
% from myF
figure(3); 
subplot(121);
imshow(I1); 
title('Inliers and Epipolar Lines in First Image (myF)'); hold on;
plot(img1pts(:,1),img1pts(:,2),'go');

epiLines1 = epipolarLine(myF', img1pts(:,1:2));
points = lineToBorderPoints(epiLines1,size(I2));
line(points(:,[1,3])',points(:,[2,4])');

subplot(122); 
imshow(I2);
title('Inliers and Epipolar Lines in Second Image (myF)'); hold on;
plot(img2pts(:,1),img2pts(:,2),'go')

epiLines2 = epipolarLine(myF,img2pts(:,1:2));
points2 = lineToBorderPoints(epiLines2,size(I2));
line(points2(:,[1,3])',points2(:,[2,4])');

% from fRANSAC
figure(4); 
subplot(121);
imshow(I1); 
title('Inliers and Epipolar Lines in First Image (fRANSAC)'); hold on;
plot(img1pts(inliers,1),img1pts(inliers,2),'go');

epiLines1 = epipolarLine(fRANSAC', img1pts(inliers,1:2));
points = lineToBorderPoints(epiLines1,size(I2));
line(points(:,[1,3])',points(:,[2,4])');

subplot(122); 
imshow(I2);
title('Inliers and Epipolar Lines in Second Image (fRANSAC)'); hold on;
plot(img2pts(inliers,1),img2pts(inliers,2),'go')

epiLines2 = epipolarLine(fRANSAC,img2pts(inliers,1:2));
points2 = lineToBorderPoints(epiLines2,size(I2));
line(points2(:,[1,3])',points2(:,[2,4])');

% epipolarlines from test points
figure(5)
subplot(121);
imshow(I1); 
title('Inliers and Epipolar Lines in First Image (test points myF)'); hold on;
plot(test_img1pts(1,:),test_img1pts(2,:),'go');
epiLines3 = epipolarLine(myF', test_img2pts(1:2,:)');
points3 = lineToBorderPoints(epiLines3,size(I2));
line(points3(:,[1,3])',points3(:,[2,4])');

subplot(122);
imshow(I2);
title('Inliers and Epipolar Lines in Second Image (test points myF)'); hold on;
plot(test_img2pts(1,:),test_img2pts(2,:),'go')
epiLines4 = epipolarLine(myF,test_img1pts(1:2,:)');
points4 = lineToBorderPoints(epiLines4,size(I2));
line(points4(:,[1,3])',points4(:,[2,4])');

figure(6)
subplot(121);
imshow(I1); 
title('Inliers and Epipolar Lines in First Image (test points fRANSAC)'); hold on;
plot(test_img1pts(1,:),test_img1pts(2,:),'go');
epiLines3 = epipolarLine(fRANSAC', test_img2pts(1:2,:)');
points3 = lineToBorderPoints(epiLines3,size(I2));
line(points3(:,[1,3])',points3(:,[2,4])');

subplot(122);
imshow(I2);
title('Inliers and Epipolar Lines in Second Image (test points fRANSAC)'); hold on;
plot(test_img2pts(1,:),test_img2pts(2,:),'go')
epiLines4 = epipolarLine(fRANSAC,test_img1pts(1:2,:)');
points4 = lineToBorderPoints(epiLines4,size(I2));
line(points4(:,[1,3])',points4(:,[2,4])');

% %% Visualize epipolar lines on images
% % el1 : epipolar line from the first image; epl1 = Fx2
% % el2 : epipolar line from the second image; epl2 = F'x1
% 
% % this is from myF
% figure(9)
% imshow(I1);
% figure(10)
% imshow(I2);
% 
% epl1 = zeros(no_TP,3);
% epl2 = zeros(no_TP,3);
% for n = 1:no_TP
%     epl1(n,:) = ((myF/myF(end))*img2pts(n,:)')';
%     epl1(n,:) = epl1(n,:)/epl1(n,3); 
%     
%     figure(9)
%     hold on;
%     x = 0:imgwidth;
%     y = (epl1(n,1)/epl1(n,2)) * x + epl1(n,3)/epl1(n,2);
%     plot(x, y, 'LineWidth', 2, 'Color', 'magenta');
%     
%     epl2(n,:) = ((myF/myF(end))'*img1pts(n,:)')';
%     epl2(n,:) = epl2(n,:)/epl2(n,3); 
%     
%     figure(10)
%     hold on;
%     x = 0:imgwidth;
%     y = (epl2(n,1)/epl2(n,2)) * x + epl2(n,3)/epl1(n,2);
%     plot(x, y, 'LineWidth', 2, 'Color', 'magenta');
% end
% hold off;
