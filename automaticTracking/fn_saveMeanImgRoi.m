function fn_saveMeanImgRoi(mouse)
dataPath = 'E:\excitData\';
disp(['Processing ' mouse]); tic;
T = func_loadMouseConfig(mouse,'root',[dataPath filesep 'config']);
allDay = unique(T.Day);
nDays = length(allDay);
nPlanes = 2;

%%
avg_day_temp = cell(nDays,nPlanes);
imgPath = [dataPath filesep mouse filesep 'meanImg'];
temp_refImg = cell(1,nPlanes);
for i = 1:nPlanes
    avg_session = load([imgPath filesep mouse '_MeanImgPerSessions_Plane' num2str(i-1) '.mat']);
    try
        avg_session = avg_session.sessionMeanImg;
    catch
        avg_session = avg_session.tosave;
        disp('MeanImg naming follows old convention');
    end
    
    for j=1:length(allDay)
        tempImg = [];
        this_day = allDay(j);
        sessionIdx_thisday = find(T.Day==this_day);
        for k=1:length(sessionIdx_thisday)
            tempImg(:,:,k) = avg_session{sessionIdx_thisday(k)};
        end
        avg_day_temp{j,i} = mean(tempImg,3)';
    end
    load([imgPath filesep mouse '_ops_plane' num2str(i-1) '.mat'],'ops');
    temp_refImg{i} = ops.refImg';
end
data_plane1.meanImg_day = fn_cell2mat(avg_day_temp(:,1),3);
data_plane2.meanImg_day = fn_cell2mat(avg_day_temp(:,2),3);
data_plane1.refImg = temp_refImg{1};
data_plane2.refImg = temp_refImg{2};
%%
roiPath = [dataPath filesep mouse filesep 'roi'];
xlen = 697; ylen = 403;

temp_roi = cell(1,nPlanes);
temp_roi_redrawn = cell(1,nPlanes);
temp_redrawnFlag = {1,1};
temp_ishere = cell(1,nPlanes);
for i = 1:nPlanes
    roisMask = loadRoi([roiPath filesep mouse '_roi' int2str(i-1) '.zip']);
    temp_roi{i} = roisMask;
    if exist([roiPath filesep 'roi_redrawn_plane' int2str(i-1) '_final.mat'])
        load([roiPath filesep 'roi_redrawn_plane' int2str(i-1) '_final.mat'], 'roi_redrawn');
        load([roiPath filesep 'ishere_plane' int2str(i-1) '_final.mat'], 'ishere');
 
        temp_roi_redrawn{i} = squeeze(roi_redrawn(:,:,2));
         
    elseif exist([roiPath filesep 'roi_redrawn_plane' int2str(i-1) '.mat'])
        load([roiPath filesep 'roi_redrawn_plane' int2str(i-1) '.mat'], 'roi_redrawn');
        load([roiPath filesep 'ishere_plane' int2str(i-1) '.mat'], 'ishere');
    else % nothing is redrawn
        load([roiPath filesep 'ishere_plane' int2str(i-1) '.mat'], 'ishere');
        temp_redrawnFlag{i} = 0;
    end
    temp_ishere{i} = ishere;
end
data_plane1.ishere = temp_ishere{1};
data_plane2.ishere = temp_ishere{2};
data_plane1.roi = temp_roi{1};
data_plane2.roi = temp_roi{2};
data_plane1.roi_redrawn = temp_redrawnFlag{1};
data_plane2.roi_redrawn = temp_redrawnFlag{2};
data_plane1.redrawFlag = temp_redrawnFlag{1};
data_plane2.redrawFlag = temp_redrawnFlag{2};

croi = zeros(length(data_plane1.roi),2);
for i = 1:length(data_plane1.roi)
    stat = regionprops(true(size(data_plane1.roi{i})),data_plane1.roi{i},'WeightedCentroid');
    try
        croi(i,:) = round(stat.WeightedCentroid);
    catch
        %croi(i,:) = (stat.Centroid);
        disp('ahhh')
    end
end 
data_plane1.croi = croi;

croi = zeros(length(data_plane2.roi),2);
for i = 1:length(data_plane2.roi)
    stat = regionprops(true(size(data_plane2.roi{i})),data_plane2.roi{i},'WeightedCentroid');
    try
        croi(i,:) = round(stat.WeightedCentroid);
    catch
        %croi(i,:) = (stat.Centroid);
        disp('ahhh')
    end
end 
data_plane2.croi = croi;

roiData = data_plane1;
save(['J:\automatic_tracking\' mouse '_plane1.mat'],'roiData');
roiData = data_plane2;
save(['J:\automatic_tracking\' mouse '_plane2.mat'],'roiData');
t = toc; disp(['Time Elapsed = ' num2str(t,'%.1f') ' sec.'] )
end
%%
%{
for k=1:nDays
    avg_day{k} = enhancedImage(avg_day{k},[1 size(avg_day{k},1)],[1 size(avg_day{k},2)]);
end
%}

function roisMask = loadRoi(filename)


xlen = 697; ylen = 403;
tempRoi = ReadImageJROI(filename);
roisCoord = cell(1,length(tempRoi));

for k = 1:length(tempRoi) % number of neurons in this plane
    % flip coord 1 and 2 here since roi mask flips x and y
    %xlen = 800; ylen = 600;
    %tempMask = roipoly(zeros(xlen,ylen),tempRoi{k}.mnCoordinates(:,2),...
    %    tempRoi{k}.mnCoordinates(:,1)); 
    tempMask = roipoly(zeros(xlen,ylen),tempRoi{k}.mnCoordinates(:,1),...
        tempRoi{k}.mnCoordinates(:,2)); 
    roisMask{k} = logical(tempMask);
    %roisCoord{k} = tempRoi{k}.mnCoordinates;
end

end