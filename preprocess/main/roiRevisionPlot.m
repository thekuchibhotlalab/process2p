function roiRevisionPlot(varargin)
%%roiTrackingBin - Description
%
% Syntax: output = roiTrackingBin(input)
%
% Long description

p = func_createInputParser();
p.parse(varargin{:});
sep = '\';
%---------CHECK NUMBER OF FRAMES IN SBX FILE-----------
global info
[nFuncChannel, functionalChannel, roiType] = func_getFuncChanRoiType(varargin{:});
filenames = strsplit(p.Results.filename);
nFiles = length(filenames);
nPlanes = str2double(p.Results.nPlanes);
mouse = p.Results.mouse;
sbxpath = p.Results.sbxpath;
suite2ppath = p.Results.suite2ppath;
h5path = p.Results.h5path;
datapath = p.Results.datapath;
behavpath = p.Results.behavpath;
tcFile = p.Results.tcFile; tcFile = reshape(strsplit(tcFile),nFuncChannel,nPlanes);
roiFile = p.Results.roiFile; roiFile = reshape(strsplit(roiFile),nFuncChannel,nPlanes);
nFrames_oneplane = p.Results.nFrames_oneplane;
nFrames_oneplane = cumsum(nFrames_oneplane,1);
nFrames_oneplane = [zeros(1,nPlanes);nFrames_oneplane];
rootpath = p.Results.root;


warning('off');

if nFuncChannel == 1
    prompt = {'Enter Plane (1/2):'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'1'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    planeToDo = str2double(answer{1});
    chanToDo = 1;
    avg_session = load([datapath sep 'meanImg' sep mouse '_MeanImgPerSessions_Plane' num2str(planeToDo-1) '.mat']);
else
    prompt = {'Enter Plane (1/2):','Enter Channel (green/red):'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'1','green'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    planeToDo = str2double(answer{1});
    chanToDo = find(strcmp(answer{2},functionalChannel));
    avg_session = load([datapath sep 'meanImg' sep mouse...
        '_MeanImgPerSessions_Plane' num2str(planeToDo-1) '_' functionalChannel{chanToDo} '.mat']);
end
try
    avg_session = avg_session.sessionMeanImg;
catch
    avg_session = avg_session.tosave;
    disp('MeanImg naming follows old convention');
end

revisePath = uigetdir();

mkdir([revisePath sep 'checkCell' sep]);
mkdir([revisePath sep 'checkFill' sep]);

plotRevise = true;
plotFilled = true;

%---------DFF CALCULATION-----------
imgData = func_loadMouseConfig(mouse,'root',rootpath);
allDay = unique(imgData.Day);
nDays = length(allDay);
avg_day = cell(nDays,1);

for i=1:length(allDay)
    this_day = allDay(i);
    sessionIdx_thisday = find(imgData.Day==this_day);
    for k=1:length(sessionIdx_thisday)
        tempImg(:,:,k) = avg_session{sessionIdx_thisday(k)};
    end
    avg_day{i,1} = mean(tempImg,3);
end
%---------ROI TRACKING CODE-----------
warning('off')
margins = [0.02 0.002];

i = planeToDo;

if ~exist([revisePath sep 'ishere_plane' num2str(i-1) '.mat'],'file')
    plotRevise = false;
    disp('Not plotting revised')
else
    ishere = load([revisePath sep  'ishere_plane' num2str(i-1) '.mat']);
    ishere = ishere.ishere;
    disp('Loading New ishere')
end

if ~exist([revisePath sep 'filled_plane' num2str(i-1) '.mat'],'file')
    plotFilled = false;
    disp('Not plotting filled')

else
    filled = load([revisePath sep  'filled_plane' num2str(i-1) '.mat']);
    filled = filled.filled; 
    disp('Loading New filled')
end

if ~exist([revisePath sep 'roi_redrawn_plane' num2str(i-1) '.mat'],'file')
    plotRevise = false;
    disp('Not plotting revised')
else
    roi_redrawn = load([revisePath sep 'roi_redrawn_plane' num2str(i-1) '.mat']);
    roi_redrawn = roi_redrawn.roi_redrawn;
    disp('Loading New roi_redrawn')
end
roiName = [datapath sep roiFile{chanToDo,planeToDo}]; %[mouse '_roi'  int2str(i-1) '.zip'];
rois = ReadImageJROI(roiName); %read imagej rois
%data = load([suite2ppath sep 'plane' num2str(i-1) sep 'Fall.mat']);
data = load([datapath sep 'meanImg' sep mouse '_ops_plane' num2str(i-1) '.mat']);

nCells = size(ishere,1);

%%    
for j=1:nCells 
    
    localImgList = {};
    
    
    yroi = rois{1,j}.mnCoordinates(:,1); %freehand rois have the outlines in x-y coordinates
    xroi = rois{1,j}.mnCoordinates(:,2); % y from 1-697, x from 1-403
    croi = [mean(minmax(xroi')) mean(minmax(yroi'))]; % croi of size 403*697
    ylimm = [croi(1)-20 croi(1)+20]; %ylimm should be from 1-403
    xlimm = [croi(2)-20 croi(2)+20]; %xlim should be 1-697
    
    for k=1:nDays     
        extraMargin = 20;
        if (ylimm(1)-extraMargin) < 1  || (xlimm(1)-extraMargin) < 1 ||...
                (ylimm(2)+extraMargin)> size(avg_day{k},1) || (xlimm(2)+extraMargin) > size(avg_day{k},2)
            extraMargin = min([ylimm(1), xlimm(1) size(avg_day{k},1)-ylimm(2)   size(avg_day{k},2)-xlimm(1)])-1;
        end
        localImg = uint16(avg_day{k}(ylimm(1)-extraMargin:ylimm(2)+extraMargin, xlimm(1)-extraMargin:xlimm(2)+extraMargin));
        localImgList{k} = localImg; 
    end 

    saveFig = figure('visible','off');
    set(saveFig, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.04, 0.85, 0.96]);
    tempCol = ceil(nDays/4);
    tight_margins = [0.001 0.001];
    for d = 1:nDays
        subplot_tight(4,tempCol,d,tight_margins)
        hold on;
        imagesc(imadjust(adapthisteq(localImgList{d}, 'NBins', 256),[0 1],[0 1],0.2)); 
        xlim([extraMargin+1 xlimm(2)-xlimm(1)+extraMargin+1]);ylim([extraMargin+1 ylimm(2)-ylimm(1)+extraMargin+1]);
        title(['D' num2str(d)]);
        colormap gray;
        axis off
        redrawFlag = roi_redrawn{j,d,1};
        roixy = roi_redrawn{j,d,3};
        pat1 = patch(roixy(:,1)-round(xlimm(1)-extraMargin)+1,roixy(:,2)-round(ylimm(1)-extraMargin)+1,...
            'b','FaceColor','None');
        if redrawFlag==1
            pat1.EdgeColor = 'b';
        elseif  ishere(j,d)
            pat1.EdgeColor = 'g';
        else
            pat1.EdgeColor = 'r';
        end
    end
    %pause(2);
    saveas(saveFig,[revisePath sep 'checkCell' sep 'check_cell_plane' num2str(i-1) '_cell' int2str(j) '.png']);
    close(saveFig)
    
    fillFig = figure('visible','off');
    set(fillFig, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.04, 0.85, 0.96]);
    tempCol = ceil(nDays/4);
    tight_margins = [0.001 0.001];
    for d = 1:nDays
        subplot_tight(4,tempCol,d,tight_margins)
        hold on;
        imagesc(imadjust(adapthisteq(localImgList{d}, 'NBins', 256),[0 1],[0 1],0.2)); 
        xlim([extraMargin+1 xlimm(2)-xlimm(1)+extraMargin+1]);ylim([extraMargin+1 ylimm(2)-ylimm(1)+extraMargin+1]);
        title(['D' num2str(d)]);
        colormap gray;
        axis off
        roixy = roi_redrawn{j,d,3};
        pat1 = patch(roixy(:,1)-round(xlimm(1)-extraMargin)+1,roixy(:,2)-round(ylimm(1)-extraMargin)+1,...
            'b','FaceColor','None');
        if filled(j,d)==1
            pat1.EdgeColor = 'r';
        else
            pat1.EdgeColor = 'g';
        end
    end
    %pause(2);
    saveas(fillFig,[revisePath sep 'checkFill' sep 'check_fill_plane' num2str(i-1) '_cell' int2str(j) '.png']);
    close(fillFig)

end
end
