close all;
clear all;
clc;
format long g
% a straight line fitting with a function of least squares estimation'LSE.m'
% Ilseo Jeon
% Dept. of Geoinformatics, The University of Seoul
% Adv. Photogrammetry Lab3 in Master Course
% 2020. 04. 20

x=[10,12,16,18]';
y=[8,6,9,10]';

A=[x, ones(4,1)];
Dy=eye(size(y,1));

%LSE
[kt, et, vct, Dkt] = LSE (A, y, Dy);
xt=linspace(0,20);
yt=kt(1)*xt+kt(2);

figure();
plot(x,y,'o',xt,yt,'LineWidth',2);
grid on
title('A Straight Line Fitting with LSE') 
