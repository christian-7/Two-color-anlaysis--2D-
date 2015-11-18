
clear 
clc
close all
%% Generate random X,Y data 

maxSim=20;
minSim=0;
locs=10000;

sim(:,1) = (maxSim-minSim).*rand(locs,1) + (0.2).*rand(locs,1);
sim(:,2) = 1*sim(:,1)+(2-minSim).*rand(locs,1);

sim2(:,1) = (maxSim-0).*rand(locs,1) + (0.2).*rand(locs,1);
sim2(:,2) = -1*sim2(:,1)+(2-minSim).*rand(locs,1)+20;

figure
scatter(sim(:,1),sim(:,2),1,'red'); hold on;
scatter(sim2(:,1),sim2(:,2),1,'black'); hold on;

all1=sim;
all2=sim2;

subset1st=sim;
subset2nd=sim2;
%% Generate random X,Y data 


maxSim=20;
minSim=0;
locs=10000;

result_x=(maxSim-minSim).*rand(locs,1) + (0.2).*rand(locs,1);
result_y=1*result_x+(2-minSim).*rand(locs,1);

result_x2=(maxSim-0).*rand(locs,1) + (0.2).*rand(locs,1);
result_y2=-1*result_x2+(2-minSim).*rand(locs,1)+20;

frames=transpose(1:10000);
int=zeros(10000,1)+1;

figure
scatter(result_x,result_y,1,'red'); hold on;
scatter(result_x2,result_y2,1,'black'); hold on;


%% Generate rapidStorm file for input in LAMA


pixelsz=1000;
localization_folder_name='K:\Christian\matlab\CBC';

        locResults(:,1) = double(result_x*pixelsz);
        locResults(:,2) = double(result_y*pixelsz);
        locResults(:,3) = double(frames);
        locResults(:,4) = double(int);
        
        locResults2(:,1) = double(result_x2*pixelsz);
        locResults2(:,2) = double(result_y2*pixelsz);
        locResults2(:,3) = double(frames);
        locResults2(:,4) = double(int);


% tabs
        hdrNames = {' # x[nm]', 'y[nm]',' t[frame]',' I[A.D. counts]'}; % header variables names
        outfileDirName = strcat( localization_folder_name , '\Channel1_', datestr(now,'yyyymmdd HHMM'), '.txt' ); % directory of outfile and name
        fmt = repmat('%s\t ', 1, length(hdrNames));
        fmt(end:end+1) = '\n';
 
%         dlmwrite(outfile, hdrNames, '');
        fid = fopen(outfileDirName, 'w'); % open output file to write headers
        fprintf(fid, fmt, hdrNames{:}); % write headers into the file
        fclose(fid);
%         fmt = [repmat('%g\t', size(locResults, 2)), '\n'];
%         fprintf(fid, fmt, transpose(locResults));
%         fclose(fid);
%        save( fullfile(localization_folder_name, strcat('ResultsTabSep', datestr(now,'yyyymmdd HHMM'))), 'locResults', '-ascii', '-double', '-tabs', '-append');
        dlmwrite(outfileDirName, locResults,  'delimiter', '\t', '-append');
        
        
% tabs
        hdrNames = {' # x[nm]', 'y[nm]',' t[frame]',' I[A.D. counts]'}; % header variables names
        outfileDirName = strcat( localization_folder_name , '\Channel2_', datestr(now,'yyyymmdd HHMM'), '.txt' ); % directory of outfile and name
        fmt = repmat('%s\t ', 1, length(hdrNames));
        fmt(end:end+1) = '\n';
 
%         dlmwrite(outfile, hdrNames, '');
        fid = fopen(outfileDirName, 'w'); % open output file to write headers
        fprintf(fid, fmt, hdrNames{:}); % write headers into the file
        fclose(fid);
%         fmt = [repmat('%g\t', size(locResults, 2)), '\n'];
%         fprintf(fid, fmt, transpose(locResults));
%         fclose(fid);
%        save( fullfile(localization_folder_name, strcat('ResultsTabSep', datestr(now,'yyyymmdd HHMM'))), 'locResults', '-ascii', '-double', '-tabs', '-append');
        dlmwrite(outfileDirName, locResults2,  'delimiter', '\t', '-append');