%%%%% script will load and correct the RED CHANNEL - channel2 according to a set x,y offset
% 
% 1. load images filtered from Thunderstorm
% 2. identify gold fiducials and select ROI around
% 3. correct for additional sample drif by adding X and Y Offset
% 4. save corrected RED channel
%
%  Data: 20/04/15
%
%%%%%
%% Clear workspace

clear all, close all, clear all, clc

%% Open files and plot scatter

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

filenameC1='FOV2_20ms_gain500_FarRed_1_crop_TS_filtered_corr';          % -->  transformed far red channel, i.e. from Trans_2D_after_TS.m
filenameC1_2=[filenameC1 '.txt'];
filenameC2='FOV2_20ms_gain500_Red_1_crop_TS_filtered';                  % -->  red channel
filenameC2_2=[filenameC2 '.txt'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

channel1=importdata(filenameC1_2);
channel2=importdata(filenameC2_2);

scatter(channel1(:,1),channel1(:,2),1);  % far red channel
hold on;
scatter(channel2.data(:,2),channel2.data(:,3),1,'red'); % red channel 7,8 for data from PS


%% Select ROI and scatter plot it --> choose area with gold fiducial

upperx=10000; % max(all(:,1));
lowerx=9000;

% vx=find(channel1.data(:,1) < upperx & channel1.data(:,1) > lowerx);
% vxC2=find(channel2.data(:,7) < upperx & channel2.data(:,7) > lowerx);                   % 7,8 for data from PS
vx=find(channel1(:,1) < upperx & channel1(:,1) > lowerx);
vxC2=find(channel2.data(:,2) < upperx & channel2.data(:,2) > lowerx);  
subset=channel1(vx);
subsetC2=channel2.data(vxC2,2); % 1
subset(:,2)=channel1(vx,2);
subsetC2(:,2)=channel2.data(vxC2,3);  % 3

uppery=9500; % max(subset(:,2))
lowery=8500;

vy=find(subset(:,2) < uppery & subset(:,2) > lowery);
vyC2=find(subsetC2(:,2) < uppery & subsetC2(:,2) > lowery);
subset2=subset(vy);
subset2C2=subsetC2(vyC2);
subset2(:,2)=subset(vy,2);
subset2C2(:,2)=subsetC2(vyC2,2);


figure
scatter(subset2(:,1), subset2(:,2),1,'black'); hold on
scatter(subset2C2(:,1), subset2C2(:,2),1,'red')             % red channel

%% Plot 2D Histogramm to identify gold fiducials

% width=ceil(max(subset2(:,1))/10);
% height=ceil(max(subset2(:,2))/10);

c1=hist3([subset2(:,1), subset2(:,2)],[(upperx-lowerx)/100,(uppery-lowery)/100]);

c2=hist3([subset2C2(:,1), subset2C2(:,2)],[(upperx-lowerx)/100,(uppery-lowery)/100]);

figure('Position', [100 200 800 900])

subplot(3,2,1)
imagesc(c1);
colormap(jet);
caxis([0 5])
colorbar;

subplot(3,2,2)
imagesc(c2);
colormap(jet);
caxis([0 5])
colorbar;

subplot(3,2,3)
contour(c1);
colormap(jet);
caxis([0 3])
colorbar;

subplot(3,2,4)
contourf(c2);
colorbar;
colormap(jet);
caxis([0 3])
colorbar;

subplot(3,2,5)
scatter(subset2(:,1), subset2(:,2),1,'black')
subplot(3,2,6)
scatter(subset2C2(:,1), subset2C2(:,2),1,'red')


%% Shift channel 1 dataset according to gold bead

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
xOff= 0;
yOff= 80;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

subset2C2Cor(:,1)=subset2C2(:,1)+xOff;
subset2C2Cor(:,2)=subset2C2(:,2)+yOff;

figure('Position', [500 300 700 700])

subplot(2,2,1)
scatter(subset2C2Cor(:,1), subset2C2Cor(:,2),1,'black')
title('Before correction')

subplot(2,2,2)
scatter(subset2C2(:,1), subset2C2(:,2),1,'red')
title('After correction')

subplot(2,2,3)
scatter(subset2(:,1), subset2(:,2),1,'black'), hold on;
scatter(subset2C2(:,1), subset2C2(:,2),1,'red')
title('Overlay before correction')

subplot(2,2,4)
scatter(subset2(:,1), subset2(:,2),1,'black');hold on;
scatter(subset2C2Cor(:,1), subset2C2Cor(:,2),1,'red');
title('Overlay after correction')


%% Apply to the full dataset and save for colocalization analysis

channel2Cor(:,1)=channel2.data(:,2)+xOff;   % 7
channel2Cor(:,2)=channel2.data(:,3)+yOff;   % 8

figure
scatter(channel2Cor(:,1),channel2Cor(:,2),1,'red');  % red channel
hold on;
scatter(channel1(:,1),channel1(:,2),1,'black'); % far red channel

filename = ['Man_Corr_' filenameC2 '.txt'];

dlmwrite(filename, channel2Cor)




