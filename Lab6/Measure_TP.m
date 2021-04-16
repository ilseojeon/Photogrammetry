clear all
close all
clc

% Read the image
im1 = imread( 'stereo-corridor_l.gif', 'gif');
im2 = imread( 'stereo-corridor_r.gif', 'gif');

% Visualize the image point over image 1
figure(1)
imshow(im1);
ax(1,:) = axis;
hold on
axis on

figure(2)
imshow(im2);
ax(2,:) = axis;
hold on
axis on

k = 0;
ws = 50;
while 1
    for n = 1:2,
        figure(n)
        [x1,y1,btn] = ginput(1);
        if btn ~= 1,
            break,
        end
        axis([x1-ws x1+ws y1-ws y1+ws]);
        [x2,y2,btn] = ginput(1);
        if btn ~= 1,
            break,
        end
        if n == 1,
            k = k + 1;
            IP_ics(k,1) = k;
        end
        IP_ics(k,n*2:n*2+1) = [x2,y2];
        plot(IP_ics(k,n*2), IP_ics(k,n*2+1), 'rx');
        h = text(IP_ics(k,n*2), IP_ics(k,n*2+1), sprintf('%d', k) );
        set(h, 'Color', 'c');
        axis ( ax(n,:) );
    end
    if btn == 3,
        break;
    end
end

fid = fopen ('stereo-corridor_test.txt', 'w');
fprintf(fid, '%% TP Coordinates in ICS\r\n' );
fprintf(fid, '%d\r\n', size(IP_ics,1) );
fprintf(fid, '%d\t%.3f\t%.3f\t%.3f\t%.3f\r\n', IP_ics' );
fclose(fid);
