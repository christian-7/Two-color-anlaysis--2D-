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

filenameC1='Man_Corr_FOV4_Gain300_20ms_Red_1_crop_TS_filtered2';                     % -->  manually corrected red channel
filename_peaks1=[filenameC1 '.txt'];
filenameC2='FOV4_Gain300_20ms_FarRed_1_crop_TS_filtered_corr';                      % -->  transformed far red channel, i.e. from Trans_2D_after_TS.m
filename_peaks2=[filenameC2 '.txt'];

peaks1=dlmread(filename_peaks1);  % from manual correction
peaks2=dlmread(filename_peaks2); 
% peaks2=dlmread(filename_peaks2,'',2,0); % from auto correction RS output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% 256   x 128 pxl    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% 25.6  x 12.8 mum   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% 25600 x 128000 nm  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

upperx=18; % max(all(:,1));
lowerx=7;

uppery=9; % max(subset(:,2))
lowery=4;

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

% 
% xdim=(upperx-lowerx)*100;
% ydim=(uppery-lowery)*100;
% 
% c=hist3([subset1st(:,2), subset1st(:,1)],[ydim xdim]);
% d=hist3([subset2nd(:,2), subset2nd(:,1)],[ydim xdim]); % heigth x width
% 
% 
% 
% % 128x128 px --> 12.800 x 12.800 nm --> 10 nm/pixel grid
% % 256x512 px --> 25.600 x 51.200 nm 
% 
% % save both channels
% 
% imwrite(c,[filenameC1 '_rendered_10nm_pxl.tiff']);
% imwrite(d,[filenameC2 '_rendered_10nm_pxl.tiff']);
% 


%% Calculate NN for all localizations in a with respect to a and b within range minD - maxD

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

minD=0.01;      % min radius
maxD=0.1;       % max radius
a=subset1st;    % Red channel --> subset1st        old:b
b=subset2nd;    % Far Red channel --> subset2nd    old:a

list=minD:0.01:maxD;
NNab=cell(length(list),2);
NNaa=cell(length(list),2);

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parfor index=1:length(list)                                               % NN search for each radius between a and b 
NNab{index,1}=vertcat(NNab{index,1},rangesearch(b, a, list(1,index)));
end

parfor index=1:length(list)                       
NNab{index,2}=vertcat(NNab{index,2},list(1,index));   
end

parfor index=1:length(list)                                               % NN searc for each radius between a and a 
NNaa{index,1}=vertcat(NNaa{index,1},rangesearch(a, a, list(1,index)));
end

parfor index=1:length(list)                       
NNaa{index,2}=vertcat(NNaa{index,2},list(1,index));   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc

%% Quantify the total number of neighbors depending on the distance for each point --> NN distribution

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preallocate empty cells

for index=1:length(NNab);
NNab{index,3}=[];
end

for index=1:length(NNaa);
NNaa{index,3}=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sum the total number of neighbors depending on the distance for each point 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for index=1:length(NNab)                % for each distance
        
    for index2=1:length(NNab{index,1})  % for each point
    
        NNab{index,3}=cat(1,NNab{index,3},length(NNab{index,1}{index2,1})); % total number of NN for each distance -> N_A,B
    end    
end

for index=1:length(NNaa)                % for each distance
      
    for index2=1:length(NNaa{index,1})  % for each point
        
        NNaa{index,3}=cat(1,NNaa{index,3},length(NNaa{index,1}{index2,1})); % total number of NN for each distance -> N_A,A      
    end   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc
%% Calculate and normalize the localization ditribution for aa and ab

%   1. correct NN for the area of the search circle
%   2. normalize by the total number and the largest search area
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preallocate empty cells

for index=1:length(NNab);
NNab{index,4}=[];
end

for index=1:length(NNaa);
NNaa{index,4}=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for index=1:length(NNab);

for index2=1:length(NNab{index,3});    

    NNab{index,4}(index2,1)=((NNab{index,3}(index2,1))/(NNab{10,3}(index2,1)))*(maxD^2/(NNab{index,2})^2); % -> D ai,b for each point (i)
    
end

end

for index=1:length(NNaa);

for index2=1:length(NNaa{index,3});    

    NNaa{index,4}(index2,1)=((NNaa{index,3}(index2,1))/(NNaa{10,3}(index2,1)))*(0.01/(NNaa{index,2})^2); % -> D ai,a for each point (i)

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc

%% Calculate the Spearman correlation coefficient and colocalization value Ca for each localization

RHO=zeros(length(NNaa{1,4}),1);
Ca=zeros(length(NNaa{1,4}),1);
[IDX,D] = knnsearch(b,a);                           % search nearest neighbor in b for each point in A and output distance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
for index=1:length(NNaa{1,4})
    
    Daa=zeros(10,1);
    Dab=zeros(10,1);

    for index2=1:length(NNaa)
        
        Daa(index2,1)=(NNaa{index2,4}(index,1));
        Dab(index2,1)=(NNab{index2,4}(index,1));
        
    end
    

    RHO(index,1)=corr(Daa,Dab,'Type','Spearman');   % correlate Daa and Dab                        
    Ca(index,1)=RHO(index)*exp((-D(index)/0.1));    % calculate Ca --> colocalization value for every single localization between -1 and 1
    
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position',[100 500 900 500])

subplot(1,2,1)
scatter(b(:,1),b(:,2),2,'black');
hold on;
scatter(a(:,1),a(:,2),2,Ca);
axis([lowerx upperx lowery uppery])
%axis([0 12.8 0 12.8])
colormap ('jet')
colorbar;
title(['Mean Ca = ',num2str(mean(Ca))])

subplot(1,2,2)
hist(Ca,10)
% axis([-1 1 0 100])
xlabel(['C_A colocalization'])
ylabel(['counts'])
title(['Spearman coefficient = ',num2str(mean(RHO))])

% select and quantify the amount of localizations considered as colocalized
% i.e. above a certain treshold

tresh=0;
vx=find(Ca(:,1) > tresh);
colocalized=length(vx)/length(Ca)

%% Save results


filename = ['CBC_CA_10-100nm_' filenameC2 '.txt'];

dlmwrite(filename, Ca)

