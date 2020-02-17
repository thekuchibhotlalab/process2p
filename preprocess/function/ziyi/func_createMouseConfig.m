function T = func_createMouseConfig(mouse, h5path)

cd(h5path);
files = dir('*.h5');
nFiles = length(files);
names = cell(nFiles,1);

for i=1:nFiles
    names{i} = files(i).name(1:end-3);
end
% get days and unique days 
for i = 1:length(names)
    tempSplit = strsplit(names{i},'_');
    days(i) = str2num(tempSplit{2});
end
uniqueDays = unique(days);
% get sessions 
sessionCount = 0;
lastDay = 0;
sessions = cell(1,length(uniqueDays));
for i = 1:length(names)
    tempSplit = strsplit(names{i},'_');
    currDay = str2num(tempSplit{2});
    if lastDay == currDay
        sessionCount = sessionCount + 1;
    else sessionCount = 1;
    end
    lastDay = currDay;
    sessions{find(currDay==uniqueDays)}= [sessions{find(currDay==uniqueDays)} sessionCount];
end

% get frames per plane for each file
nPlanes = 2;
nFrames_oneplane = zeros(length(names),2);
nFrames = zeros(length(names),1);
global info
for i=1:length(names)
    sbxread(names{i},1,1);
    nFrames(i) = info.max_idx;

    if mod(nFrames(i),2)
        nFrames_oneplane(i,:) = [round(nFrames(i)/nPlanes) round(nFrames(i)/nPlanes)-1];
    else
        nFrames_oneplane(i,:) = [nFrames(i)/nPlanes nFrames(i)/nPlanes];
    end
end


% then create the txt file
%fileID = fopen([mouse '_config.txt'], 'wt');
T = cell(nFiles,12);
count = 0;
for i = 1:length(uniqueDays)
    daySessCount = 0;
    for j = 1:length(sessions{i})
        count = count + 1;
        daySessCount = daySessCount +1;
        tempSplit = strsplit(names{count},'_');
        
        
        if uniqueDays(i) == 0 % for day 0, tuning
            if sessions{i}(j) == 1
                T(count,1:11) = {str2double(tempSplit{2}), sessions{i}(j),  'Tuning PuretoneBef',names{count}, '', 0, 0, 0, 0, 0, 0};
                %tempLine = ['Day ' str2double(str2double(tempSplit{2})) ' Session ' sessions{i}(j) ' Imaging ' names{count} ' Tuning PuretoneBef'];
            elseif sessions{i}(j) >= 2
                T(count,1:11) = {str2double(tempSplit{2}), sessions{i}(j),  'Tuning WhiteNoiseBef',names{count}, '', 0, 0, 0, 0, 0, 0};
                %tempLine = ['Day ' str2double(str2double(tempSplit{2})) ' Session ' sessions{i}(j) ' Imaging ' names{count} ' Tuning WhiteNoiseBef'];
            end
        elseif i == length(uniqueDays)-1 % post-learning D1 tuning curve sessions
            if sessions{i}(j) == 1
                T(count,1:11) = {str2double(tempSplit{2}), sessions{i}(j),  'Tuning PuretoneAftD1',names{count}, '', 0, 0, 0, 0, 0, 0};
                %tempLine = ['Day ' str2double(str2double(tempSplit{2})) ' Session ' sessions{i}(j) ' Imaging ' names{count} ' Tuning PuretoneAftD1'];
            elseif sessions{i}(j) >= 2
                T(count,1:11) = {str2double(tempSplit{2}), sessions{i}(j),  'Tuning WhiteNoiseAftD1',names{count}, '', 0, 0, 0, 0, 0, 0};
                %tempLine = ['Day ' str2double(str2double(tempSplit{2})) ' Session ' sessions{i}(j) ' Imaging ' names{count} ' Tuning WhiteNoiseAftD1'];
            end
        elseif i == length(uniqueDays) % post-learning D7 tuning curve sessions
            if sessions{i}(j) == 1
                T(count,1:11) = {str2double(tempSplit{2}), sessions{i}(j),  'Tuning PuretoneAftD7',names{count}, '', 0, 0, 0, 0, 0, 0};
                %tempLine = ['Day ' str2double(str2double(tempSplit{2})) ' Session ' sessions{i}(j) ' Imaging ' names{count} ' Tuning PuretoneAftD7'];
            elseif sessions{i}(j) >= 2
                T(count,1:11) = {str2double(tempSplit{2}), sessions{i}(j),  'Tuning WhiteNoiseAftD7',names{count}, '', 0, 0, 0, 0, 0, 0};
                %tempLine = ['Day ' str2double(str2double(tempSplit{2})) ' Session ' sessions{i}(j) ' Imaging ' names{count} ' Tuning WhiteNoiseAftD7'];
            end
        else % all behavioral days
            if j == length(sessions{i}) % baseline sessions
                T(count,1:11) = {str2double(tempSplit{2}), sessions{i}(j),  'Baseline',names{count}, '', 0, 0, 0, 0, 0, 0};
                %tempLine = ['Day ' str2double(str2double(tempSplit{2})) ' Session ' sessions{i}(j) ' Imaging ' names{count} ' Baseline ' temBehFileName];
            else % behavioral sessions    
                tempBehFileName = [tempSplit{1} '_' int2str(str2double(tempSplit{2})) 'v' int2str(daySessCount) '.txt']; % infer behavioral file name         
                T(count,1:11) = {str2double(tempSplit{2}), sessions{i}(j),  'Behavior',names{count}, tempBehFileName, 0, 0, 0, 0, 0, 0};
                %tempLine = ['Day ' str2double(tempSplit{2}) ' Session ' sessions{i}(j) ' Imaging ' names{count} ' Behavior ' temBehFileName];
            end
        end
        T(count,12) = {strjoin(strsplit(int2str(nFrames_oneplane(count,:))))};
        %tempLine = [tempLine ' Alignment Nan ROI Nan TC Nan Neuropil Nan Deconvolution Nan'];

        %fprintf(fileID,[tempLine '\n']);
    end


end
T = cell2table(T);
T.Properties.VariableNames = {'Day','Session','BehavType','ImagingFile',...
    'BehavFile','Alignment', 'ROI', 'TC', 'Neuropil', 'dff', 'Deconvolution',...
    'nFrames_oneplane'};
%fclose(fileID);
writetable(T,[mouse '_config.csv']);

    
    

end