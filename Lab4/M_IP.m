% Stage 2. Measuring IPs
% Ilseo Jeon
% Dept. of Geoinformatics, The University of Seoul
% Adv. Photogrammetry Lab4 in Master Course
% 2020. 05. 7

% Read the image
function [IP_ccs, IP_ics] = M_IP(img, IO, resolution)
im = imread( img, 'jpg');
R_IC=[1 0 ; 0 -1];
img_origin =[resolution(1)/2 resolution(2)/2];

% Visualize the image point over image 1
figure(1)
imshow(im);
ax = axis;
hold on
axis on

k = 0; %
ws = 50; %
while 1
    figure(1)
    [x1,y1,btn] = ginput(1); %ginput(n) n개 점 좌표 식별 
    if btn == 3 % click three times, break
        break
    end
    k = k + 1;
    axis([x1-ws x1+ws y1-ws y1+ws]);
    [x2,y2,btn] = ginput(1);
    
    IP_ics(k,:) = [k, x2, y2];
    plot(IP_ics(k,2), IP_ics(k,3), 'rx');
    h = text(IP_ics(k,2), IP_ics(k,3), sprintf('%d', k) );
    set(h, 'Color', 'c');
    axis ( ax );

end

for n=1:size(IP_ics(:,1))
    IP_ccs(n,:)=[R_IC*((IP_ics(n,2:3)-img_origin)*IO(4))'; -IO(3)]';
end

fid = fopen ('IP_ics.txt', 'w');
fprintf(fid, '%% IP Coordinates in ICS\r\n' );
fprintf(fid, '%d\r\n', size(IP_ics,1) );
fprintf(fid, '%d\t%.3f\t%.3f\r\n', IP_ics' );
fclose(fid);

fid = fopen ('IP_pcs.txt', 'w');
fprintf(fid, '%% IP Coordinates in ICS\r\n' );
fprintf(fid, '%d\r\n', size(IP_ccs,1) );
fprintf(fid, '%d\t%.3f\t%.3f\r\n', IP_ccs' );
fclose(fid);
end
