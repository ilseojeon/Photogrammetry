%% Feature matching with SURF
% 
% from Lab6, Adv. Photogrammetry, Geoinformatics, Univ. of Seoul
% Ilseo Jeon

%% matching with SURF
function [matchedPoints1, matchedPoints2]=SURF(I1, I2, str, m)

img1pts = detectSURFFeatures(I1, 'MetricThreshold',300);
img2pts = detectSURFFeatures(I2, 'MetricThreshold',300);

img1pts_str = img1pts.selectStrongest(str);
img2pts_str = img2pts.selectStrongest(str);

figure(m);
subplot(1,2,1); imshow(I1); hold on
plot(img1pts_str);
subplot(1,2,2); imshow(I2); hold on
plot(img2pts_str);
hold off

[f1,vpts1] = extractFeatures(I1,img1pts_str);
[f2,vpts2] = extractFeatures(I1,img2pts_str);

indexPairs = matchFeatures(f1,f2, 'Method', 'Approximate','Unique', true, 'MaxRatio', 0.2,'MatchThreshold', 10) ;
matchedPoints1 = vpts1(indexPairs(:,1));
matchedPoints2 = vpts2(indexPairs(:,2));

figure(m+1); showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2);
legend('matched points 1','matched points 2');
hold off;
end
