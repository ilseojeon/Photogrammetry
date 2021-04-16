close all;
clear all;
clc;

% Prob2 with a function of least squares estimation'LSE.m'
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% the observation equations for the horizontal network adjustment 
% to determine the locations of multiple cars simultaneouly. 
% The observations can be assumed as 
% (1) the absolute location of each car using its own GPS, 
% (2) the absolute locations of static objects (ex. traffic light poles) on the road from a high-precision road map, 
% and (3) the relative distances and (or) directions of the other cars 
% or the static objects measured from each car using its on-board sensors such as cameras, lidars, radars and so on.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Ilseo Jeon
% Dept. of Geoinformatics, The University of Seoul
% Adv. Photogrammetry Lab3 in Master Course
% 2020. 04. 20

%% real location of each car a, b, c
real_loc=[800,200;
          300,600;
          200,800];
no_data=3;

%% weight unit:m
G_Q = 20; % from GPS
O_Q = 10; % from relative distances 
T_Q  = 8; % from static objects(traffic light poles)

%% observation data
gps_loc=real_loc+normrnd(0, G_Q);  % 3x2

Oab=real_loc(2,:)-real_loc(1,:)+normrnd(0, O_Q); % 1x2
Obc=real_loc(3,:)-real_loc(2,:)+normrnd(0, O_Q); % 1x2
Oca=real_loc(1,:)-real_loc(3,:)+normrnd(0, O_Q); % 1x2

T(1,:)=real_loc(1,:)+normrnd(0, T_Q); % 1x2
T(2,:)=real_loc(2,:)+normrnd(0, T_Q); % 1x2

%% LSE
A = [1 0 0; 
     0 1 0;
     0 0 1;
     -1 1 0;
     0 -1 1;
     1 0 -1;
     1 0 0;
     0 1 0];  
y=[gps_loc;Oab;Obc;Oca;T(1,:);T(2,:)];
D=[G_Q, G_Q, G_Q, O_Q, O_Q, O_Q, T_Q, T_Q];
Dy=eye(size(A,1)).*D;

[kt_x, et_x, vct_x, Dkt_x] = LSE (A, y(:,1), Dy);
[kt_y, et_y, vct_y, Dkt_y] = LSE (A, y(:,2), Dy);

%% visualization
figure()
hold on
for i=1:no_data
    plot(real_loc(i,1),real_loc(i,2),'ro','MarkerSize', 10, 'DisplayName', sprintf('real location of car %d',i));
    plot(gps_loc(i,1),gps_loc(i,2),'g*', 'MarkerSize', 10,'DisplayName', sprintf('gps location of car %d',i));
    plot(kt_x(i), kt_y(i),'b.', 'MarkerSize', 10,'DisplayName', sprintf('adjusted location of car %d',i));
end
daspect([1, 1, 1])
grid on
legend
title('Horizontal Network Adjustment  with LSE') 