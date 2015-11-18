%%%%% Loads the data after manual correct and outputs:
% 
%  1. rendered imaged of both channels with 10 nm/pxl
%  2. coordinate-based colocalization (CBC)
%     - saves CA histogram
%     - fraction of localizations having CA > than treshold
%     - mean CA
%
%  Data: 20/04/15
%%%%%
%% Clear workspace

clear all, close all, clear all, clc

%% if coming directly from analysis_1.m

pix=0.001;
all2(:,1)=pix.*channel1.data(:,1);
all2(:,2)=pix.*channel1.data(:,3);
all1(:,1)=pix.*channel2Cor(:,1);                  % corrected Red channel
all1(:,2)=pix.*channel2Cor(:,2);  

%% Open both channel localizations

filenameC1='Man_Corr_FOV1_Gain300_20ms_Red_1_crop_TS_filtered';             % -->  manually corrected red channel
filename_peaks1=[filenameC1 '.txt'];
filenameC2='FOV1_Gain300_20ms_FarRed_1_crop_TS_filtered_corr';                     % -->  transformed far red channel, i.e. from Trans_2D_after_TS.m
filename_peaks2=[filenameC2 '.txt'];

peaks1=dlmread(filename_peaks1);  % from manual correction
peaks2=dlmread(filename_peaks2); 
% peaks2=dlmread(filename_peaks2,'',2,0); % from auto correction RS output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% 256   x 128 pxl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% 25.6  x 12.8 mum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% 25600 x 128000 nm%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assign first channel

pix=0.001; % 0.1 for pxl to ?m -- 0.001 for nm to ?m
sdx=pix.*nonzeros(peaks1(:,1));% 3,20 --> rapidStorm 1  --> from analysis_1.m 1
sdy=pix.*nonzeros(peaks1(:,2));% 4,21 --> rapidStorm 3  --> from analysis_1.m 2
all1(:,1)=sdx;
all1(:,2)=sdy;

all1=unique(all1,'rows');

% asign second channel

sdx2=pix.*nonzeros(peaks2(:,1));% 3,20 --> rapidStorm 1
sdy2=pix.*nonzeros(peaks2(:,2));% 4,21 --> rapidStorm 3
all2(:,1)=sdx2;
all2(:,2)=sdy2;

all2=unique(all2,'rows');

% plot full dataset

figure
scatter(all1(:,1),all1(:,2),1)
hold on;
scatter(all2(:,1),all2(:,2),1)

%% Select region to analyze colocalization

%%%%%%%%%%%%%%%%%%%% select region %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

upperx=15; % max(all(:,1));
lowerx=5;

uppery=20; % max(subset(:,2))
lowery=5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% First channel

vx=find(all1(:,1) < upperx & all1(:,1) > lowerx);
subset=all1(vx);
subset(:,2)=all1(vx,2);

vy=find(subset(:,2) < uppery & subset(:,2) > lowery);
subset1st=subset(vy);
subset1st(:,2)=subset(vy,2);

clear subset vx vy

% Second channel

vx=find(all2(:,1) < upperx & all2(:,1) > lowerx);
subset=all2(vx);
subset(:,2)=all2(vx,2);

vy=find(subset(:,2) < uppery & subset(:,2) > lowery);
subset2nd=subset(vy);
subset2nd(:,2)=subset(vy,2);

% plot subset of dataset

figure
scatter(subset1st(:,1), subset1st(:,2),1,'red'); hold on;
scatter(subset2nd(:,1), subset2nd(:,2),1,'black');

%% save 2D Histogram of both channels

xdim=(upperx-lowerx)*100;
ydim=(uppery-lowery)*100;

c=hist3([subset1st(:,2), subset1st(:,1)],[ydim xdim]);
d=hist3([subset2nd(:,2), subset2nd(:,1)],[ydim xdim]); % heigth x width



% 128x128 px --> 12.800 x 12.800 nm --> 10 nm/pixel grid
% 256x512 px --> 25.600 x 51.200 nm 

% save both channels

imwrite(c,[filenameC1 '_rendered_10nm_pxl.tiff']);
imwrite(d,[filenameC2 '_rendered_10nm_pxl.tiff']);



%% Calculate NN for all localizations in a with respect to a and b within range minD - maxD

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

minD=0.1;%0.01
maxD=1;%0.1
b=subset1st; % Red channel --> subset1st        old:b
a=subset2nd; % Far Red channel --> subset2nd    old:a

list=minD:0.1:maxD;
NNab=cell(length(list),2);
NNaa=cell(length(list),2);

tic

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

parfor index=1:1:length(list)                       % NN search between a and b 
NNab{index,1}=vertcat(NNab{index,1},rangesearch(a, b, list(1,index)));
end

parfor index=1:1:length(list)                       
NNab{index,2}=vertcat(NNab{index,2},list(1,index));   
end

parfor index=1:1:length(list)                       % NN search between a and a 
NNaa{index,1}=vertcat(NNaa{index,1},rangesearch(a, a, list(1,index)));
end

parfor index=1:1:length(list)                       
NNaa{index,2}=vertcat(NNaa{index,2},list(1,index));   
end

toc

%% Quantify the total number of neighbors depending on the distance --> NN distribution

tic

NNab{1,3}=[]; % total number of NN for each distance N_A,B
NNaa{1,3}=[]; % total number of NN for each distance N_A,A

for index=1:length(NNab)
    
    for index2=1:length(NNab{index,1})
    
        NNab{index,3}=cat(1,NNab{index,3},length(NNab{index,1}{index2,1}));

        index2=index2+1;
    end

    NNab{index,4}=sum(NNab{index,3});
    index=index+1;
    
end

for index=1:length(NNaa)
    
    for index2=1:length(NNaa{index,1})
    
        NNaa{index,3}=cat(1,NNaa{index,3},length(NNaa{index,1}{index2,1}));

        index2=index2+1;
    end

    NNaa{index,4}=sum(NNaa{index,3});
    index=index+1;
    
end

NNab2(:,1)=cell2mat(NNab(:,2)); % radius
NNab2(:,2)=cell2mat(NNab(:,4)); % number of neighbors


NNaa2(:,1)=cell2mat(NNaa(:,2));
NNaa2(:,2)=cell2mat(NNaa(:,4));

toc
%% Calculate and normalize the localization ditribution for aa and ab

%   1. correct NN for the area of the search circle
%   2. normalize by the total number and the largest search area

rmaxab=max(NNab2(:,1));
target=find(NNab2(:,1)==rmaxab);
NRmax=NNab2(target);
NRmax=sum(NNab2(target,2));


for index=1:length(NNab2);
    

Dab(index,1)=(NNab2(index,2))/(pi*NNab2(index,1)^2)*((pi*((rmaxab)^2))/(NRmax));

%Dab(index,1)=((NNab2(index,2))*pi*((rmaxab)^2))/((pi*NNab2(index,1))*(NRmax));
%Dab(index,1)=(NNab2(index,2)*(rmaxab)^2)/((NRmax)*(NNab2(index,1)^2));
%Dab(index,1)=(NNab2(index,2)/NRmax)*(((rmaxab)^2)/(NNab2(index,1))^2);

end

clear NRmax 

rmaxaa=max(NNaa2(:,1));
target=find(NNaa2(:,1)==rmaxaa);
NRmax=NNaa2(target);
NRmax=sum(NNaa2(target,2));


for index=1:length(NNaa2);

        Daa(index,1)=(NNaa2(index,2))/(pi*NNaa2(index,1)^2)*((pi*((rmaxaa)^2))/(NRmax));
        %Daa(index,1)=((NNaa2(index,2))*pi*((rmaxab)^2))/((pi*NNaa2(index,1))*(NRmax));
        %Daa(index,1)=(NNaa2(index,2)*(rmaxab)^2)/((NRmax)*(NNaa2(index,1)^2));
        %Daa(index,1)=(NNaa2(index,2)/NRmax)*(((rmaxaa)^2)/(NNaa2(index,1))^2);

end

%% Calculate the correlation coefficient and colocalization value Ca for each localization

RHO = corr(Daa,Dab,'Type','Spearman');      % correlate Daa and Dab
[IDX,D] = knnsearch(b,a);                   % search nearest neighbor and output distance
Ca=RHO*exp((-D/rmaxab));                    % calculate Ca --> colocalization value for every single localization between -1 and 1

figure('Position',[100 500 900 250])

subplot(1,2,1)
scatter(b(:,1),b(:,2),1,'black');
hold on;
scatter(a(:,1),a(:,2),1,Ca);
% axis([lowerx upperx lowery uppery])
%axis([0 12.8 0 12.8])
colormap ('jet')
colorbar;
title(['Mean Ca = ',num2str(mean(Ca))])


subplot(1,2,2)
hist(Ca,50)
axis([-1 1 0 10000])
xlabel(['C_A colocalization'])
ylabel(['counts'])
title(['Spearman coefficient = ',num2str(RHO)])

% select and quantify the amount of localizations considered as colocalized
% i.e. above a certain treshold

tresh=0.6;
vx=find(Ca(:,1) > tresh);
colocalized=length(vx)/length(Ca)

%% Save results


filename = ['CBC_CA_10-100nm_' filenameC2 '.txt'];

dlmwrite(filename, Ca)

