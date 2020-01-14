function roiTracking(varargin)
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
% Default values
acq = 30.98;
pretone = 1; % for PSTH tone, in s
posttone = 3; % for PSTH tone, in s
% Column names
SESSION = 1;
TRIAL = 2;
TONE = 3;
RESPONSE = 4;
H = 1; M = 2; FA = 3; CR = 4;
NOLICK_PERIOD = 5;
RESP_TIME = 6;
DELAY_AFTER_RESP = 7;
TOTAL_TRIAL_DUR_MINUS_RESP_TIME = 8;
LICKF = 9;
REWARDF = 10;
TOTAL_TRIAL_DUR = 11; % toc 
TONEF = 12;
CONTEXT=13;
sep = '\';
%load(['E:\KishoreLab\Shared\Matlab\preprocessing-newpipeline-draft\preprocessing\config\mouse\' mouse '_ImagingMetaData.mat']);
if nFuncChannel == 1
    prompt = {'Enter Plane:'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'1'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    planeToDo = str2double(answer{1});
    chanToDo = 1;
    avg_session = load([datapath sep 'meanImg' sep mouse '_MeanImgPerSessions_Plane' num2str(planeToDo-1) '.mat']);
else
    prompt = {'Enter Plane:','Enter Channel (green/red):'};
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
    avg_session = avg_session.sessionMeanImg;s
catch
    avg_session = avg_session.tosave;
    disp('MeanImg namign follows old convention');
end
%configPath = 'E:\KishoreLab\Shared\Matlab\preprocessing-newpipeline-draft\preprocessing\config';
%mouseDataPath = [configPath sep 'mouse' sep mouse];
roipath = [datapath sep 'roi' sep];
mkdir(roipath)
mkdir([roipath sep 'checkCell' sep]);
%mkdir([mouseDataPath sep 'checkCell' sep])

%cd(suite2ppath);

imgData = func_loadMouseConfig(mouse,'root',rootpath);
allDay = unique(imgData.Day);
behavDay = unique(imgData.Day(strcmp(imgData.BehavType,'Behavior')));

nDays = length(allDay);
avg_day = cell(nDays,1);
foil_fr = [];target_fr = [];

for i=1:length(allDay)
    this_day = allDay(i);
    sessionIdx_thisday = find(imgData.Day==this_day);
    for k=1:length(sessionIdx_thisday)
        tempImg(:,:,k) = avg_session{sessionIdx_thisday(k)};
    end
    avg_day{i,1} = mean(tempImg,3);
end

b_count = 0;
for i=1:length(allDay)
    this_day = allDay(i);
    sessionIdx_thisday_behav = find(imgData.Day == this_day & strcmp(imgData.BehavType,'Behavior'));
    sessionIdx_thisday = find(imgData.Day==this_day);
    list_b{i} = b_count+(1:length(sessionIdx_thisday_behav)); % update the number of behavioral files today
    list_f{i} = sessionIdx_thisday;
    b_count = b_count + length(sessionIdx_thisday_behav);
end

%---------TONE FRAME FOR EACH SESSION-----------
for i = 1:length(behavDay)
    this_day = behavDay(i);
    sessionIdx_thisday_behav = find(imgData.Day == this_day & strcmp(imgData.BehavType,'Behavior'));
    sessionIdx_thisday = find(imgData.Day==this_day);
    for j = 1:length(sessionIdx_thisday_behav)
        % REVISE THIS LINE
        behavMatrix = load([behavpath sep imgData.BehavFile{sessionIdx_thisday_behav(j)}]);
        nTrials = size(behavMatrix,1);

        foil = unique([behavMatrix(behavMatrix(:,RESPONSE)==FA,TONE);behavMatrix(behavMatrix(:,RESPONSE)==CR,TONE)]);
        target = unique([behavMatrix(behavMatrix(:,RESPONSE)==H,TONE);behavMatrix(behavMatrix(:,RESPONSE)==M,TONE)]);
                
        foil_frames = behavMatrix(ismember(behavMatrix(:,TONE),foil),TONEF);
        target_frames = behavMatrix(ismember(behavMatrix(:,TONE),target),TONEF);
        nFoil = size(foil_frames,1);
        nTarget = size(target_frames,1);
        
        foilFramePlane = [foil_frames foil_frames];
        foilFramePlane(logical(mod(foil_frames,2)),:) = [round(foil_frames(logical(mod(foil_frames,2)))/nPlanes) ...
            round(foil_frames(logical(mod(foil_frames,2)))/nPlanes)-1];
        foilFramePlane(~mod(foil_frames,2),:) = [foil_frames(~mod(foil_frames,2))/nPlanes ...
            foil_frames(~mod(foil_frames,2))/nPlanes];

        targetFramePlane = [target_frames target_frames];
        targetFramePlane(logical(mod(target_frames,2)),:) = [round(target_frames(logical(mod(target_frames,2)))/nPlanes) ...
            round(target_frames(logical(mod(target_frames,2)))/nPlanes)-1];
        targetFramePlane(~mod(target_frames,2),:) = [target_frames(~mod(target_frames,2))/nPlanes ...
            target_frames(~mod(target_frames,2))/nPlanes];

        if ~(j==1 && i==1)
            foilFramePlane = foilFramePlane+repelem(nFrames_oneplane(sessionIdx_thisday(j),:),nFoil,1);
            targetFramePlane = targetFramePlane+repelem(nFrames_oneplane(sessionIdx_thisday(j),:),nTarget,1);
        end
        foil_fr = [foil_fr;foilFramePlane];
        target_fr = [target_fr;targetFramePlane];
    end
    %list_b{i} = sessionIdx_thisday_behav;
    %list_f{i} = sessionIdx_thisday;
end
start_foil = foil_fr - pretone*round(acq); % - pretone sec before tone onset
start_target = target_fr - pretone*round(acq); 
nframes_psth = pretone*round(acq) + posttone*round(acq); 
start_foil = foil_fr - pretone*round(acq); % - pretone sec before tone onset
start_target = target_fr - pretone*round(acq); 
nframes_psth = pretone*round(acq) + posttone*round(acq); 

%---------DFF CALCULATION-----------
try
    mat = load([datapath sep tcFile{chanToDo,planeToDo}],'TC'); mat = mat.TC;
catch
    mat = load([datapath sep tcFile{chanToDo,planeToDo}],'tempTC'); mat = mat.tempTC;
    disp('TC naming follows old convention')
end

for j=1:nFiles
    submat = mat(:,nFrames_oneplane(j,planeToDo)+1:nFrames_oneplane(j+1,planeToDo));
    mat(:,nFrames_oneplane(j,planeToDo)+1:nFrames_oneplane(j+1,planeToDo),2) = ...
        (submat-median(submat,2))./median(submat,2);%./repmat(max(submat,[],2)-min(submat,[],2),1,size(submat,2));
end

nFrames_add = sum(nFrames_oneplane(2:end,:),2);
nFrames = diff([0;nFrames_add]);
baseline = ismember(diff(nFrames_add),9999);
baseline = logical([0;0;baseline]);
nBaseline = sum(baseline);
tuningsessions = ismember(nFrames,[4999;16999]);
tuningsessions = logical([0;tuningsessions]);
nTuning = sum(tuningsessions);
intervals_baseline = cell(nPlanes,1);
intervals_tc = cell(nPlanes,1);
for i=1:nPlanes
    intervals_baseline{i} = [nFrames_oneplane(circshift(baseline,-1),i) nFrames_oneplane(baseline,i)];
    intervals_tc{i} = [nFrames_oneplane(circshift(tuningsessions,-1),i) nFrames_oneplane(tuningsessions,i)];
end

%---------ROI TRACKING CODE-----------
ishere = [];
warning('off')
margins = [0.02 0.002];
%for i=1:nPlanes
i = planeToDo;
global tempIsCell
global roi_redrawn
global tempRoiCoord
global origRoiCoord
global rois
global currentCopy
global pat_day
global copyCellNum

zfluo = squeeze(mat(:,:,2)); % 1D=fluo; 2D=df/f
nCells = size(zfluo,1);

if ~exist([roipath sep 'ishere_plane' num2str(i-1) '.mat'],'file')
    ishere = nan(nCells,nDays);
    c = 1;
else
    
    ishere = load([roipath sep  'ishere_plane' num2str(i-1) '.mat']);
    ishere = ishere.ishere;       
    c = size(ishere(~any(isnan(ishere),2),:),1)+1;
    if c==0
        c = 1;
    end
end

if ~exist([roipath sep 'roi_redrawn_plane' num2str(i-1) '.mat'],'file')
    roi_redrawn =  cell(nCells,nDays,3);
else
    roi_redrawn = load([roipath sep 'roi_redrawn_plane' num2str(i-1) '.mat']);
    roi_redrawn = roi_redrawn.roi_redrawn;
end

%if ~exist([roipath sep 'roi_redrawnFlag_plane' num2str(i-1) '.mat'],'file')
%    roi_redrawnFlag =  zeros(nCells,nDays);
%else
%    roi_redrawnFlag = load([roipath sep 'roi_redrawnFlag_plane' num2str(i-1) '.mat']);
%    roi_redrawnFlag = roi_redrawnFlag.roi_redrawnFlag;
%end

roiName = [datapath sep roiFile{chanToDo,planeToDo}]; %[mouse '_roi'  int2str(i-1) '.zip'];
rois = ReadImageJROI(roiName); %read imagej rois
data = load([suite2ppath sep 'plane' num2str(i-1) sep 'Fall.mat']);
%%    
disp(['cell' int2str(c)])
for j=c:nCells  
    tic;
    tempCoord = [rois{1,j}.mnCoordinates(:,1),rois{1,j}.mnCoordinates(:,2)];
    bw = roipoly(zeros(data.ops.Lx,data.ops.Ly),rois{1,j}.mnCoordinates(:,2),rois{1,j}.mnCoordinates(:,1));
    [xt,yt] = find(bw==1);
    tempRoiCoord = cell(nDays,3);
    tempRoiCoord (1:nDays,1) = {0};
    tempRoiCoord (1:nDays,2) = {[xt,yt]};
    tempRoiCoord (1:nDays,3) = {tempCoord};
    origRoiCoord = tempRoiCoord;
    currentCopy = [];
    
    figdays = figure;         
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.04, 0.85, 0.96]);% Enlarge figure to full screen.

    % Plot norm traces
    subplot_tight(4,1,1,margins);hold on;title(['Cell ' num2str(j) ', raw traces']);
    plot(zfluo(j,:));        
    for k=1:nBaseline
        plot(intervals_baseline{i}(k,1)+1:intervals_baseline{i}(k,2),zfluo(j,intervals_baseline{i}(k,1)+1:intervals_baseline{i}(k,2)),'g');% color the baseline traces
    end
    temp = nFrames_oneplane;
    nFrames_oneplane = temp;
    PlotHVLines(nFrames_oneplane(2:end,i),'v','color',[0.8 0.8 0.8]);

    xlim([0 size(zfluo(j,:),2)]);
    for k=1:nTuning
        plot(intervals_tc{i}(k,1)+1:intervals_tc{i}(k,2),zfluo(j,intervals_tc{i}(k,1)+1:intervals_tc{i}(k,2)),'color',[0.8 0.8 0.8]);
    end        

    % Plot tone-evoked resp        
    %[~,whichsessions_foil] = InIntervals(start_foil(:,i),[nFrames_oneplane(2:end-1,i) nFrames_oneplane(3:end,i)]);
    %[~,whichsessions_target] = InIntervals(start_target(:,i),[nFrames_oneplane(2:end-1,i) nFrames_oneplane(3:end,i)]);
    [~,whichsessions_foil] = InIntervals(start_foil(:,i),[nFrames_oneplane(1:end-1,i) nFrames_oneplane(2:end,i)]);
    [~,whichsessions_target] = InIntervals(start_target(:,i),[nFrames_oneplane(1:end-1,i) nFrames_oneplane(2:end,i)]);
    for k=1:nDays
        if ~ismember(allDay(k),behavDay), continue, end
        subplot_tight(4,nDays*2,nDays*2+allDay(k)*2+1,margins);hold on;
        ok_foil = ismember(whichsessions_foil,list_b{k});
        l = sum(ok_foil);
        matidx = repelem(0:nframes_psth-1,l,1);
        idx = (matidx+repelem(start_foil(ok_foil,i),1,nframes_psth))';
        m_foil = zfluo(j,idx(:));
        traces_foil = reshape(m_foil,nframes_psth,l);
        plot(traces_foil,'color',[0.8 0.8 0.8]);plot(mean(traces_foil,2),'k','linewidth',2);
        title(['D' num2str(k)]);
        if allDay(k)==behavDay(1), ylimm1 = ylim; else, ylim(ylimm1); end
%             PlotIntervals([pretone*round(acq) pretone*round(acq)+0.1*round(acq)],'color','r');
        PlotHVLines(pretone*round(acq),'v','color','r','linewidth',1);
        if k>1, axis off; end

        subplot_tight(4,nDays*2,nDays*2+allDay(k)*2+2,margins);hold on;
        ok_target = ismember(whichsessions_target,list_b{k});
        l = sum(ok_target);
        matidx = repelem(0:nframes_psth-1,l,1);
        idx = (matidx+repelem(start_target(ok_target,i),1,nframes_psth))';
        m_target = zfluo(j,idx(:));
        traces_target = reshape(m_target,nframes_psth,l);
        plot(traces_target,'color',[0.8 0.8 0.8]);plot(mean(traces_target,2),'k','linewidth',2);
        ylim(ylimm1);
%             PlotIntervals([pretone*round(acq) pretone*round(acq)+0.1*round(acq)],'color','g');
        PlotHVLines(pretone*round(acq),'v','color','g','linewidth',1);
        title(['D' num2str(k)]);
        axis off;
    end

    % Plot this ROI in the field across days    
    yroi = rois{1,j}.mnCoordinates(:,1); %freehand rois have the outlines in x-y coordinates
    xroi = rois{1,j}.mnCoordinates(:,2);
    croi = [mean(minmax(xroi')) mean(minmax(yroi'))];
    ylimm = [croi(1)-20 croi(1)+20];
    xlimm = [croi(2)-20 croi(2)+20];
    nCol = round(nDays/2);
    % plot it first on the mean img, then across days

    subplot_tight(4,nCol,nCol*4-nCol*2,margins); hold on; 
    imagesc(data.ops.meanImgE);
%         scatter(yroi,xroi,'g.');
    patch(yroi,xroi,'g','FaceColor','none');
    plot(round(croi(2)),round(croi(1)),'r+');
    xlim(xlimm);ylim(ylimm);
    axis off;
    localImgList = {};
    marginFlag = 0;
    subplot_day = cell(1,nDays);
    for k=1:nDays
        if isempty(avg_day{k,i}), continue;end
        subplot_day{k} = subplot_tight(4,nCol,nCol*4-nCol*2+k,margins);      

%             subimg = avg_day{k,i}(ylimm(1):ylimm(2),xlimm(1):xlimm(2));
%             subimgadj = imadjust(adapthisteq(uint16(subimg),'NBins',512));
%             hold on;imagesc(subimgadj);
%             plot(croi(2)-xlimm(1),croi(1)-ylimm(1),'r+');
%             xlim([0 croi(2)-xlimm(1)+20]);ylim([0 croi(1)-ylimm(1)+20]);
%             colormap gray;

%             imgg = imadjust(adapthisteq(uint16(avg_day{k,i})));
%             title(['D' num2str(k)]);
%             hold on;imagesc(imgg);
%             colormap gray;
%             axis off;            
%             plot(croi(2),croi(1),'r+');
%             xlim(xlimm);ylim(ylimm);

        % Try new contrast adjustment method
        hold on;
        extraMargin = 20;
        if (ylimm(1)-extraMargin) < 1  || (xlimm(1)-extraMargin) < 1 ||...
                (ylimm(2)+extraMargin)> size(avg_day{k,i},1) || (xlimm(2)+extraMargin) > size(avg_day{k,i},2)
            extraMargin = min([ylimm(1), xlimm(1) size(avg_day{k,i},1)-ylimm(2)   size(avg_day{k,i},2)-xlimm(1)])-1;
            marginFlag = 1;
        end
        localImg = uint16(avg_day{k,i}(ylimm(1)-extraMargin:ylimm(2)+extraMargin, xlimm(1)-extraMargin:xlimm(2)+extraMargin));
        localImgList{k} = localImg;
        imagesc(imadjust(adapthisteq(localImg, 'NBins', 256),[0 1],[0 1],0.2)); 
        xlim([extraMargin+1 xlimm(2)-xlimm(1)+extraMargin+1]);ylim([extraMargin+1 ylimm(2)-ylimm(1)+extraMargin+1]);
        axis off;
        title(['D' num2str(k)]);
        colormap gray;
        plot(extraMargin+1+croi(2)-xlimm(1),extraMargin+1+croi(1)-ylimm(1),'r+');
        pat_day{k} = patch(yroi-round(xlimm(1)-extraMargin)+1, xroi-round(ylimm(1)-extraMargin)+1,'g', 'FaceColor','None');
        %pat.FaceAlpha = 0.3;
    end    
    if marginFlag == 1
        disp('cell on edge, adjust margin')
    end
    % Open GUI for selecting days 
%         listDays = {};
%         for k=1:nDays, listDays{k} = ['Day ' num2str(k)]; end
%         titleWindow = '';
%         dims = [1 20];  
%         definput = {};
%         for k=1:nDays, definput{k} = [num2str(1)]; end
%         answer = inputdlg(listDays,titleWindow,dims,definput);
%         answr = str2double(answer)';

    % open the NEW GUI for selecting days and redrawing ROIs

    %screenSize = get(0,'screensize'); 
    f_selection = figure; 
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.15, 0.96]);

    %set(gcf, 'Position',  [100, 100, 200, screenSize(4)- 200]);

    figureH = 0.96; %screenSize(4)-200;
    items = nDays;
    itemH = figureH/(2*items+3);
    itemHPos = linspace(0,figureH, 2*items+4);

    buttonHPos = itemHPos(2);
    barHPos = itemHPos(4:2:end);

    startLoc = 0.02;%10;
    leng = 0.16;%40;

    % get all the txts
    for k = 1:length(barHPos)-1
        cpanel = uicontrol(f_selection,'Style','text'); %cpanel.Position = [startLoc barHPos(end-k) leng itemH];
        cpanel.String = ['day' int2str(k)];
        set(cpanel,'Units', 'normalized','Position',[startLoc barHPos(end-k) leng itemH])
    end   

    % get all the editable 
    startLoc = 0.2; %60;
    leng = 0.18; %60;
    edit_box = {};
    for k = 1:length(barHPos)-1
        cpanel = uicontrol(f_selection,'Style','edit'); 
        %cpanel.Position = [startLoc barHPos(end-k) leng itemH];
        set(cpanel,'Units', 'normalized','Position',[startLoc barHPos(end-k) leng itemH])
        cpanel.String = '1';
        edit_box{k} = cpanel;
    end
    cpanel = uicontrol(f_selection,'Style','pushbutton','CallBack',{@continue_callback,f_selection, edit_box});
    %cpanel.Position = [startLoc buttonHPos leng itemH];
    set(cpanel,'Units', 'normalized','Position',[startLoc buttonHPos leng itemH])
    cpanel.String = 'next cell';

    % get all the redraw
    startLoc = 0.40; %130;
    leng = 0.18;%60;
    redraw_button = {};
    for k = 1:length(barHPos)-1
        cpanel = uicontrol(f_selection,'Style','togglebutton','CallBack',{@redraw_callback, k,j,{localImgList{k}, extraMargin, xlimm, ylimm, xroi, yroi, croi}});
        %cpanel.Position = [startLoc barHPos(end-k) leng itemH];
        set(cpanel,'Units', 'normalized','Position',[startLoc barHPos(end-k) leng itemH])
        cpanel.String = 'redraw';
        redraw_button{k} = cpanel;
    end
    cpanel = uicontrol(f_selection,'Style','pushbutton','CallBack',{@reject_callback,f_selection,items});
    %cpanel.Position = [startLoc buttonHPos leng itemH];
    set(cpanel,'Units', 'normalized','Position',[startLoc buttonHPos leng itemH])
    cpanel.String = 'reject all';
    
    % get all the copies
    startLoc = 0.60;
    leng = 0.18;
    copy_button = {};
    for k = 1:length(barHPos)-1
        cpanel = uicontrol(f_selection,'Style','pushbutton','CallBack',{@copy_callback, k,j,...
            {localImgList{k}, extraMargin, xlimm, ylimm, xroi, yroi, croi}});
        %cpanel.Position = [startLoc barHPos(end-k) leng itemH];
        set(cpanel,'Units', 'normalized','Position',[startLoc barHPos(end-k) leng itemH])
        cpanel.String = 'copy';
        copy_button{k} = cpanel;
    end
    cpanel = uicontrol(f_selection,'Style','pushbutton','CallBack',{@preview_callback,figdays,...
        subplot_day,{localImgList{k}, extraMargin, xlimm, ylimm, xroi, yroi, croi}});
    %cpanel.Position = [startLoc buttonHPos leng itemH];
    set(cpanel,'Units', 'normalized','Position',[startLoc buttonHPos leng itemH])
    cpanel.String = 'preview';
    
    % get all the pastes
    startLoc = 0.80;
    leng = 0.18;
    paste_button = {};
    for k = 1:length(barHPos)-1
        cpanel = uicontrol(f_selection,'Style','togglebutton','CallBack',{@paste_callback, k,j,...
            {localImgList{k}, extraMargin, xlimm, ylimm, xroi, yroi, croi}});
        %cpanel.Position = [startLoc barHPos(end-k) leng itemH];
        set(cpanel,'Units', 'normalized','Position',[startLoc barHPos(end-k) leng itemH])
        cpanel.String = 'paste';
        paste_button{k} = cpanel;
    end
    cpanel = uicontrol(f_selection,'Style','pushbutton','CallBack',{@previewB_callback,...
        {localImgList, extraMargin, xlimm, ylimm, xroi, yroi, croi}});
    set(cpanel,'Units', 'normalized','Position',[startLoc buttonHPos leng itemH])
    cpanel.String = 'previewB';
    
    
    uiwait(f_selection);
    %tempIsCell = str2double(tempIsCell);

    % draw red stuff on the rejected cells
    nop = find(tempIsCell==0);
    %figure(figdays);
    for this=1:length(nop)
        k = nop(this);
        %subplot(4,nCol,nCol*4-nCol*2+k);
        %patch(yroi,xroi,'r','EdgeColor','none');
        pat = patch(subplot_day{k},yroi-round(xlimm(1)-extraMargin)+1, xroi-round(ylimm(1)-extraMargin)+1,'r','EdgeColor','none');
        pat.FaceAlpha = 0.3;
    end
    pause(2);
    %ishere(j,:) = answr;
    ishere(j,:) = tempIsCell;
    
    for d = 1:nDays
        %tempCellDay = [j d];
        %existCellDay = cellfun(@isequal,roi_redrawn(:,1),tempCellDay,'UniformOutput',true);
        roi_redrawn(j,d,:) = tempRoiCoord(d,:);
        %roi_redrawnFlag(j,d) = tempRoiCoord{d,1};
    end
    
    %disp(['should be saved: ' num2str(tempIsCell)]);
    close(figdays);

    screensize = get( groot, 'Screensize' );
    saveFig = figure('Position',[screensize(3)*0.05,screensize(4)*0.4,...
        screensize(3)*0.9,screensize(4)*0.5]);

    tempCol = ceil(nDays/2);
    %if ~isempty(roi_redrawn)
    %    roiRedrawnCellDay = cell2mat(roi_redrawn(:,1));
    %    roiRedrawnCellIndex = find(roiRedrawnCellDay(:,1)==j);% changed by Celine, 12/13/2019
    %    roiRedrawnDay = roiRedrawnCellDay(roiRedrawnCellDay(:,1)==j,2);% changed by Celine, 12/13/2019
    %else 
%             roiRedrawnCellDay = [];% changed by Celine, 12/13/2019
    %    roiRedrawnDay = []; % changed by Celine, 12/13/2019
    %end
%         roiRedrawnCellIndex = find(roiRedrawnCellDay(:,1)==j);% changed by Celine, 12/13/2019
%         roiRedrawnDay = roiRedrawnCellDay(roiRedrawnCellDay(:,1)==j,2);% changed by Celine, 12/13/2019
    for d = 1:nDays
        subplot_tight(2,tempCol,d,margins)
        hold on;
        imagesc(imadjust(adapthisteq(localImgList{d}, 'NBins', 256),[0 1],[0 1],0.2)); 
        xlim([extraMargin+1 xlimm(2)-xlimm(1)+extraMargin+1]);ylim([extraMargin+1 ylimm(2)-ylimm(1)+extraMargin+1]);
        title(['D' num2str(d)]);
        colormap gray;
        axis off
        redrawFlag = tempRoiCoord{d,1};
        roixy = tempRoiCoord{d,3};
        pat1 = patch(roixy(:,1)-round(xlimm(1)-extraMargin)+1,roixy(:,2)-round(ylimm(1)-extraMargin)+1,...
            'b','FaceColor','None');
        if redrawFlag==1
            pat1.EdgeColor = 'b';
        elseif  tempIsCell(d)
            pat1.EdgeColor = 'g';
        else
            pat1.EdgeColor = 'r';
        end
        
        
        %if any(roiRedrawnDay==d)
        %    roiRedrawnDayIndex = (roiRedrawnDay==d);
        %    roiRedrawnIndex = roiRedrawnCellIndex(roiRedrawnDayIndex);
        %    roixy = roi_redrawn{roiRedrawnIndex,3};
            %roimask = zeros(data.ops.Lx,data.ops.Ly);
            %for r = 1:size(roixy,1); roimask(roixy(r,2)-round(xlimm(1)-extraMargin)+1,...
            %        roixy(r,1)-round(ylimm(1)-extraMargin)+1) = 1; end
            %[roimaskx,roimasky] = find(roimask==0);
            %bound = bwboundaries(roimask);
            %pat = patch(roixy(:,2)-round(xlimm(1)-extraMargin)+1, roixy(:,1)-round(ylimm(1)-extraMargin)+1,'b');
            %pat.FaceAlpha = 0.2;
            %for b = 1:size(bound{1},1)
            %    dist = (roimaskx - bound{1}(b,1)).^2 + (roimasky - bound{1}(b,2)).^2;
            %    [~,closeIndex] = min(dist);
            %    bound{1}(b,:) = [roimaskx(closeIndex) roimasky(closeIndex)];
            %end              
            %if length(bound)>1
            %    pat1 = patch(bound{1}(:,1),bound{1}(:,2),'b');
            %    pat1.FaceAlpha = 0.2;
            %    for b = 1:size(bound{2},1)
            %        dist = (roimaskx - bound{2}(b,1)).^2 + (roimasky - bound{2}(b,2)).^2;
            %        [~,closeIndex] = min(dist);
            %        bound{2}(b,:) = [roimaskx(closeIndex) roimasky(closeIndex)];
            %    end
            %    pat2 = patch(bound{2}(:,1),bound{2}(:,2),'b');
            %    pat2.FaceAlpha = 0.0;
            %else
            %pat1 = patch(roixy(:,1)-round(xlimm(1)-extraMargin)+1,roixy(:,2)-round(ylimm(1)-extraMargin)+1,'b','EdgeColor','None');
            %pat1.FaceAlpha = 0.2;
            %end
            %set(pat,'AlphaData',0.2)
            %pat.MarkerFaceAlpha = 0.2;
            %pat.MarkerEdgeAlpha = 0.2;
        %else
        %    if tempIsCell(d)
        %        pat = patch(yroi-round(xlimm(1)-extraMargin)+1, xroi-round(ylimm(1)-extraMargin)+1,'g','EdgeColor','None');
        %        pat.FaceAlpha = 0.2;
        %    else 
        %        pat = patch(yroi-round(xlimm(1)-extraMargin)+1, xroi-round(ylimm(1)-extraMargin)+1,'r','EdgeColor','None');
        %        pat.FaceAlpha = 0.2;
        %    end
        %end
    end
    pause(2);
    saveas(saveFig,[roipath sep 'checkCell' sep 'check_cell_plane' num2str(i-1) 'cell' int2str(j) '.png']);
    close(saveFig)

    % Save indiv plane
    save([roipath sep 'ishere_plane' num2str(i-1) '.mat'],'ishere');
    save([roipath sep 'roi_redrawn_plane' num2str(i-1) '.mat'],'roi_redrawn');
    %save([roipath sep 'roi_redrawnFlag_plane' num2str(i-1) '.mat'],'roi_redrawnFlag');

    %save([suite2ppath sep 'plane' num2str(i-1) '\roi_redrawn_plane' num2str(i-1) '.mat'],'roi_redrawn');
    %save([roipath sep 'roi_redrawn_plane' num2str(i-1) '.mat'],'roi_redrawn');
end


%try
%    %cd(suite2ppath);
%    save([roipath sep 'ishere_combined.mat'],'ishere');
%catch
%    keyboard
%end


end

% GUI functions

function  reject_callback(hObject,eventdata,f_selection,items)
    global tempIsCell 
    tempIsCell = zeros(1,items);
    %disp(tempIsCell)
    %disp('reject')
    close(f_selection)
end


function  continue_callback(hObject,eventdata,f_selection, edit_box)
    global tempIsCell
    tempIsCell = [];
    for i = 1:length(edit_box)
        tempStr = get(edit_box{i}, 'String');
        tempIsCell(i) = str2double(tempStr);
    end
    %disp('continue')
    %disp(tempIsCell)
    close(f_selection)
end

function redraw_callback(hObject,eventdata,nDay,nCell,imgData)
    global tempRoiCoord
    %global roi_redrawn
    %global rois
    global origRoiCoord
    button_state = get(hObject,'Value');
    if button_state % if triggered from off --> on, redraw roi
        h_redraw = figure;
        set(h_redraw, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.04, 0.54+0.15, 0.96]);
        localImg = imgData{1};
        extraMargin = imgData{2};
        xlimm = imgData{3};
        ylimm = imgData{4};
        xroi = imgData{5};
        yroi = imgData{6};
        croi = imgData{7};
        
        hold on;
        imagesc(imadjust(adapthisteq(localImg, 'NBins', 256),[0 1],[0 1],0.2)); 
        xlim([extraMargin+1 xlimm(2)-xlimm(1)+extraMargin+1]);ylim([extraMargin+1 ylimm(2)-ylimm(1)+extraMargin+1]);
        pat = patch(yroi-round(xlimm(1)-extraMargin)+1, xroi-round(ylimm(1)-extraMargin)+1, 'g', 'FaceColor','None');
        plot(extraMargin+1+croi(2)-xlimm(1),extraMargin+1+croi(1)-ylimm(1),'r+');
        %pat.FaceAlpha = 0.3;
        title(['D' num2str(nDay)]);
        colormap gray;
        axis off
        
        h_roi = imfreehand;
        try
            roiMask = h_roi.createMask;
            roiCoord = round(h_roi.getPosition);
            roiCoord_revised =  zeros(size(roiCoord));
            roiCoord_revised(:,1) =  roiCoord(:,1)+round(xlimm(1)-extraMargin) - 1;
            roiCoord_revised(:,2) =  roiCoord(:,2)+round(ylimm(1)-extraMargin) - 1;
            
            [roix, roiy] = find(roiMask);
        
            roix = roix + round(ylimm(1)-extraMargin) - 1;
            roiy = roiy + round(xlimm(1)-extraMargin) - 1;

            roixy = [roiy,roix];
            %disp(nDay)
            uiwait(h_redraw);

            %---------------- everything here can be in the main program 
            %tempCellDay = [nCell nDay];tempC = cell(size(roi_redrawn,1),1); tempC(:) = {tempCellDay};
            %if ~isempty(roi_redrawn)
            %    existCellDay = cellfun(@isequal,roi_redrawn(:,1),tempC,'UniformOutput',true);
            %else 
            %    existCellDay = [];
            %end
            %if sum(existCellDay) >= 1 % if the roi has been redrawn, repalce it
                %idx = find(existCellDay==1);
            %    disp('repeated roi in:')
            %    disp(roi_redrawn{existCellDay,1})
            %    roi_redrawn{existCellDay,2} = roixy;
            %    roi_redrawn{existCellDay,3} = roiCoord_revised;
            %else % if roi not redrawn, create it
            %    roi_redrawn = [roi_redrawn;{tempCellDay,roixy,roiCoord_revised}];
            %end
            %----------------
            tempRoiCoord{nDay,1} = 1;
            tempRoiCoord{nDay,2} = roixy;
            tempRoiCoord{nDay,3} = roiCoord_revised;
        catch
            disp('no roi redrawn')
            set(hObject,'Value',0); 
        end
        
    else % if triggered from on --> off, delete redrawn roi
        %---------------- everything here can be in the main program 
        %tempCellDay = [nCell nDay];tempC = cell(size(roi_redrawn,1),1); tempC(:) = {tempCellDay};
        %existCellDay = cellfun(@isequal,roi_redrawn(:,1),tempC,'UniformOutput',true);
        %if sum(existCellDay) >= 1 % if the roi has been redrawn, delete it
        %    roi_redrawn(existCellDay,:) = [];
        %end
        %temp = [rois{1,nCell}.mnCoordinates(:,1);rois{1,nCell}.mnCoordinates(:,2)];
        %----------------
        tempRoiCoord{nDay,1} = 0;
        tempRoiCoord(nDay,2:3) = origRoiCoord(nDay,2:3);
        %tempRoiCoord{nDay,3} = temp;
        
    end
end

function copy_callback(hObject,eventdata,nDay,nCell,imgData)
    global currentCopy
    global tempRoiCoord
    global copyCellNum
    currentCopy = tempRoiCoord(nDay,2:3);
    copyCellNum = nDay;
end

function paste_callback(hObject,eventdata,nDay,nCell,imgData)
    global currentCopy
    %global rois
    global tempRoiCoord
    %global roi_redrawn
    global origRoiCoord
    global copyCellNum
    button_state = get(hObject,'Value');
    if button_state
        tempRoiCoord{nDay,1} = 1;
        tempRoiCoord(nDay,2:3) = currentCopy;
        %---------------- everything here can be in the main program 
        %tempCellDay = [nCell nDay];tempC = cell(size(roi_redrawn,1),1); tempC(:) = {tempCellDay};
        %if ~isempty(roi_redrawn)
        %    existCellDay = cellfun(@isequal,roi_redrawn(:,1),tempC,'UniformOutput',true);
        %else 
        %    existCellDay = [];
        %end
        %if sum(existCellDay) >= 1 % if the roi has been redrawn, repalce it
            %idx = find(existCellDay==1);
        %    disp('repeated roi in:')
        %    disp(roi_redrawn{existCellDay,1})
        %    roi_redrawn(existCellDay,2:3) = currentCopy;
            %roi_redrawn{existCellDay,3} = roiCoord;
        %else % if roi not redrawn, create it
        %    roi_redrawn = [roi_redrawn;{tempCellDay,currentCopy{:} }];
        %end
        %----------------
        hObject.String = ['day' int2str(copyCellNum)];
        
    else
        %---------------- everything here can be in the main program 
        %tempCellDay = [nCell nDay];tempC = cell(size(roi_redrawn,1),1); tempC(:) = {tempCellDay};
        %existCellDay = cellfun(@isequal,roi_redrawn(:,1),tempC,'UniformOutput',true);
        %if sum(existCellDay) >= 1 % if the roi has been redrawn, delete it
        %    roi_redrawn(existCellDay,:) = [];
        %end
        %----------------
        tempRoiCoord{nDay,1} = 0;
        tempRoiCoord(nDay,2:3) = origRoiCoord(nDay,2:3);
        hObject.String = 'paste';
    end

end

function preview_callback(hObject,eventdata,figdays,subplot_day,imgData)
    global tempRoiCoord
    global pat_day
    localImg = imgData{1};
    extraMargin = imgData{2};
    xlimm = imgData{3};
    ylimm = imgData{4};
    xroi = imgData{5};
    yroi = imgData{6};
    croi = imgData{7};
    %figure(figdays);
    for k = 1:size(tempRoiCoord,1)
        temp = tempRoiCoord{k,3};
        delete(pat_day{k})
        pat_day{k} = patch(subplot_day{k},temp(:,1)-round(xlimm(1)-extraMargin)+1,...
            temp(:,2)-round(ylimm(1)-extraMargin)+1,'g','FaceColor','None');
    end
end

function previewB_callback(hObject,eventdata,imgData)
    global tempRoiCoord
    localImgList = imgData{1};
    extraMargin = imgData{2};
    xlimm = imgData{3};
    ylimm = imgData{4};
    xroi = imgData{5};
    yroi = imgData{6};
    croi = imgData{7};
    figPrevB = figure;
    margins = 0.001;
    set(figPrevB, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.04, 0.85, 0.96]);
    
    tempCol = ceil(size(tempRoiCoord,1)/4);
    for k = 1:size(tempRoiCoord,1)
        subplot_tight(4,tempCol,k,margins)
        hold on;
        imagesc(imadjust(adapthisteq(localImgList{k}, 'NBins', 256),[0 1],[0 1],0.2)); 
        xlim([extraMargin+1 xlimm(2)-xlimm(1)+extraMargin+1]);ylim([extraMargin+1 ylimm(2)-ylimm(1)+extraMargin+1]);
        title(['D' num2str(k)]);
        colormap gray;
        axis off
        temp = tempRoiCoord{k,3};
        patch(temp(:,1)-round(xlimm(1)-extraMargin)+1,...
            temp(:,2)-round(ylimm(1)-extraMargin)+1,'g','FaceColor','None');
    end
end
    