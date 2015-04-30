%% Open files and plot scatter

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

filenameC1='FOV4_Gain_300_20ms_FarRed_1_crop_TS_full_corr';          % -->  transformed far red channel, i.e. from Trans_2D_after_TS.m
filenameC1_2=[filenameC1 '.txt'];
filenameC2='FOV4_Gain_300_20ms_Red_1_crop_TS_full';                  % -->  red channel
filenameC2_2=[filenameC2 '.txt'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

channel1=importdata(filenameC1_2);
channel2=importdata(filenameC2_2);

scatter(channel1(:,1),channel1(:,2),1);  % far red channel
hold on;
scatter(channel2.data(:,2),channel2.data(:,3),1,'red'); % red channel 7,8 for data from PS


%% Select ROI and scatter plot it --> choose area with gold fiducial

upperx=18000; % max(all(:,1));
lowerx=5000;

% vx=find(channel1.data(:,1) < upperx & channel1.data(:,1) > lowerx);
% vxC2=find(channel2.data(:,7) < upperx & channel2.data(:,7) > lowerx);                   % 7,8 for data from PS
vx=find(channel1(:,1) < upperx & channel1(:,1) > lowerx);
vxC2=find(channel2.data(:,2) < upperx & channel2.data(:,2) > lowerx);  
subset=channel1(vx);
subsetC2=channel2.data(vxC2,2); % 1
subset(:,2)=channel1(vx,2);
subsetC2(:,2)=channel2.data(vxC2,3);  % 3

uppery=9000; % max(subset(:,2))
lowery=3000;

vy=find(subset(:,2) < uppery & subset(:,2) > lowery);
vyC2=find(subsetC2(:,2) < uppery & subsetC2(:,2) > lowery);
subset2=subset(vy);
subset2C2=subsetC2(vyC2);
subset2(:,2)=subset(vy,2);
subset2C2(:,2)=subsetC2(vyC2,2);


figure
scatter(subset2(:,1), subset2(:,2),1,'black'); hold on
scatter(subset2C2(:,1), subset2C2(:,2),1,'red')     

%% 

xdim=(upperx-lowerx)/10%*100;
ydim=(uppery-lowery)/10%*100;

c=hist3([subset2(:,2), subset2(:,1)],[ydim xdim]);
d=hist3([subset2C2(:,2), subset2C2(:,1)],[ydim xdim]); % heigth x width



% 128x128 px --> 12.800 x 12.800 nm --> 10 nm/pixel grid
% 256x512 px --> 25.600 x 51.200 nm 

% save both channels

imwrite(c,[filenameC1 '_full_rendered_10nm_pxl.tiff']);
imwrite(d,[filenameC2 '_full_rendered_10nm_pxl.tiff']);
