function output = roiTrackingBinNew(mouse)
% Code is last revised on 11/18/2019

%%roiTrackingBin - Description
%
% Syntax: output = roiTrackingBin(input)
%
% Long description

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
load(['config\mouse\' mouse '_ImagingMetaData.mat']);

cd(suite2ppath);

imgData = func_loadMouseConfig(mouse);
allDay = unique(imgData.Day);
behavDay = unique(imgData.Day(strcmp(imgData.BehavType,'Behavior')));

nDays = length(allDay);
avg_day = cell(nDays,2);
foil_fr = [];target_fr = [];
% get the mean image for each session
for j=1:nPlanes
    cd([suite2ppath sep 'plane' num2str(j-1)])
    avg_session = load([mouse '_MeanImgPerSessions_Plane' num2str(j) '.mat']);
    avg_session = avg_session.tosave;
    for i=1:length(allDay)
        %this_day = dayss_behav(i);
        this_day = allDay(i);
        %if length(num2str(this_day))==2, toCompare = {[mouse '_0' num2str(this_day)]};
        %else, toCompare = {[mouse '_00' num2str(this_day)]}; end
        %toCmpr = {};
        %for k=1:nFiles, toCmpr = [toCmpr;toCompare]; end       
        %ok = cellfun(@strfind,names,toCmpr,'UniformOutput',false); 
        %ok2 = nan(nFiles,1);
        %for k=1:nFiles, if isempty(ok{k}), continue, end; ok2(k) = ok{k}; end
        %n_rec = nansum(ok2);
        %list_rec = find(ok2==1);
        %disp(['D' num2str(i) ', files ' num2str(list_rec')]);

        sessionIdx_thisday = find(imgData.Day==this_day);
        for k=1:length(sessionIdx_thisday)
            tempImg(:,:,k) = avg_session{sessionIdx_thisday(k)};
        end
        avg_day{i,j} = mean(tempImg,3);
        
        sessionIdx_thisday_behav = find(imgData.Day == this_day & strcmp(imgData.BehavType,'Behavior'));
        sessionIdx_thisday = find(imgData.Day==this_day);
        list_b{i} = sessionIdx_thisday_behav;
        list_f{i} = sessionIdx_thisday;
    end
end

% get tone frames for each behavioral session
for i = 1:length(behavDay)
    this_day = behavDay(i);
    sessionIdx_thisday_behav = find(imgData.Day == this_day & strcmp(imgData.BehavType,'Behavior'));
    
    sessionIdx_thisday = find(imgData.Day==this_day);

    for j = 1:length(sessionIdx_thisday_behav)
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

% Compute everything possible before plotting for speed purpose
mat = cell(nPlanes,1);
    
for i=1:nPlanes
    cd([suite2ppath sep 'plane' num2str(i-1)]);
    tc = load([mouse '_TC_plane' num2str(i-1) '.mat']);
    mat{i} = tc.tempTC; 
    for j=1:nFiles
        submat = mat{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i));
        mat{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i),2) = ...
            (submat-median(submat,2))./median(submat,2);%./repmat(max(submat,[],2)-min(submat,[],2),1,size(submat,2));
    end 
end
    
baseline = ismember(diff(nFrames_add),9999);
baseline = logical([0;0;baseline]);
nBaseline = sum(baseline);
tuningsessions = ismember(nFrames,[4999;16999]);
tuningsessions = logical([0;tuningsessions]);
nTuning = sum(tuningsessions);
% xBaseline = cell(nPlanes,1);
% xTuning = cell(nPlanes,1);
intervals_baseline = cell(nPlanes,1);
intervals_tc = cell(nPlanes,1);

for i=1:nPlanes
    intervals_baseline{i} = [nFrames_oneplane(circshift(baseline,-1),i) nFrames_oneplane(baseline,i)];
    intervals_tc{i} = [nFrames_oneplane(circshift(tuningsessions,-1),i) nFrames_oneplane(tuningsessions,i)];
%     intervals_baseline = [nFrames_oneplane(circshift(baseline,-1),i) nFrames_oneplane(baseline,i)];
%     for k=1:sum(baseline)
%         xBaseline{i} = [xBaseline{i};(intervals_baseline(k,1)+1:intervals_baseline(k,2))'];
%     end 
%     intervals_tc = [nFrames_oneplane(circshift(tuningsessions,-1),i) nFrames_oneplane(tuningsessions,i)];
%     for k=1:sum(tuningsessions)
%         xTuning{i} = [xTuning{i};(intervals_tc(k,1)+1:intervals_tc(k,2))'];
%     end  
end


%% 
ishere = cell(nPlanes,1);
checkmode = false;
warning('off')

for i=1:nPlanes

    global tempIsCell
    global roi_redrawn
    roi_redrawn = {[],[]}; % CHANGE ROI REDRAWN TO LOAD WHEN THE FILE EXISTS!!!!!!!
    
    cd([suite2ppath sep 'plane' num2str(i-1)]);
    m = mat{i}; % 1D=fluo; 2D=df/f
    zfluo = squeeze(m(:,:,2));
    nCells = size(zfluo,1);
    if ~exist([suite2ppath sep 'plane' num2str(i-1) '\ishere_test_plane' num2str(i-1) '.mat'],'file')
        ishere{i} = nan(nCells,nDays);
        c = 1;
    else
        ishere = load([suite2ppath sep 'plane' num2str(i-1) '\ishere_test_plane' num2str(i-1) '.mat']);
        ishere = ishere.ishere;       
        c = size(ishere{i}(~any(isnan(ishere{i}),2),:),1);
        if c==0
            c = 1;
        end
    end
    if ~exist([suite2ppath sep 'plane' num2str(i-1) '\roi_redrawn_plane' num2str(i-1) '.mat'],'file')
        roi_redrawn =  {[],[]};
    else
        roi_redrawn = load([suite2ppath sep 'plane' num2str(i-1) '\roi_redrawn_plane' num2str(i-1) '.mat']);
        roi_redrawn = roi_redrawn.roi_redrawn;
    end
    
    roiName = [mouse '_roi' int2str(i-1) '.zip'];
    rois = ReadImageJROI(roiName); %read imagej rois
    data = load('Fall.mat');
%%    
    disp(c)
    for j=c:nCells  
        tic;
        figdays = figure;         
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);% Enlarge figure to full screen.
        
        % Plot norm traces
        subplot(4,1,1);hold on;title(['Cell ' num2str(j) ', raw traces']);
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
        [~,whichsessions_foil] = InIntervals(start_foil(:,i),[nFrames_oneplane(2:end-1,i) nFrames_oneplane(3:end,i)]);
        [~,whichsessions_target] = InIntervals(start_target(:,i),[nFrames_oneplane(2:end-1,i) nFrames_oneplane(3:end,i)]);
        for k=1:nDays
            if ~ismember(k-1,behavDay), continue, end
            subplot(4,nDays*2,nDays*2+(k-1)*2+1);hold on;
            ok_foil = ismember(whichsessions_foil,list_b{k});
            l = sum(ok_foil);
            matidx = repelem(0:nframes_psth-1,l,1);
            idx = (matidx+repelem(start_foil(ok_foil,i),1,nframes_psth))';
            m_foil = zfluo(j,idx(:));
            traces_foil = reshape(m_foil,nframes_psth,l);
            plot(traces_foil,'color',[0.8 0.8 0.8]);plot(mean(traces_foil,2),'k','linewidth',2);
            title(['D' num2str(k)]);
            if k-1==behavDay(1), ylimm1 = ylim; else, ylim(ylimm1); end
%             PlotIntervals([pretone*round(acq) pretone*round(acq)+0.1*round(acq)],'color','r');
            PlotHVLines(pretone*round(acq),'v','color','r','linewidth',1);
            if k>1, axis off; end
            
            subplot(4,nDays*2,nDays*2+(k-1)*2+2);hold on;
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
        toc;
        
        tic;
        % Plot this ROI in the field across days    
        yroi = rois{1,j}.mnCoordinates(:,1); %freehand rois have the outlines in x-y coordinates
        xroi = rois{1,j}.mnCoordinates(:,2);
        croi = [mean(minmax(xroi)) mean(minmax(yroi))];
        ylimm = [croi(1)-20 croi(1)+20];
        xlimm = [croi(2)-20 croi(2)+20];
        nCol = round(nDays/2);
        % plot it first on the mean img, then across days
        
        subplot(4,nCol,nCol*4-nCol*2); hold on; 
        imagesc(data.ops.meanImgE);
%         scatter(yroi,xroi,'g.');
        patch(yroi,xroi,'g','EdgeColor','none');
        plot(round(croi(2)),round(croi(1)),'r+');
        xlim(xlimm);ylim(ylimm);
        axis off;
        localImgList = {};
        for k=1:nDays
            if isempty(avg_day{k,i}), continue;end
            subplot(4,nCol,nCol*4-nCol*2+k);      
            
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
                disp('cell on edge, adjust margin')
            end
            localImg = uint16(avg_day{k,i}(ylimm(1)-extraMargin:ylimm(2)+extraMargin, xlimm(1)-extraMargin:xlimm(2)+extraMargin));
            localImgList{k} = localImg;
            imagesc(imadjust(adapthisteq(localImg, 'NBins', 256),[0 1],[0 1],0.2)); 
            xlim([extraMargin+1 xlimm(2)-xlimm(1)+extraMargin+1]);ylim([extraMargin+1 ylimm(2)-ylimm(1)+extraMargin+1]);
%             axis off;
            title(['D' num2str(k)]);
            colormap gray;
            plot(extraMargin+1+croi(2)-xlimm(1),extraMargin+1+croi(1)-ylimm(1),'r+');
        end    
        toc;
        
        tic;
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
        
        screenSize = get(0,'screensize'); 
        f_selection = figure; set(gcf, 'Position',  [100, 100, 200, screenSize(4)- 200]);

        figureH = screenSize(4)-200;
        items = nDays;
        itemH = figureH/(2*items+3);
        itemHPos = linspace(0,figureH, 2*items+4);

        buttonHPos = itemHPos(2);
        barHPos = itemHPos(4:2:end);
        
        startLoc = 10;
        leng = 40;

        % get all the txts
        for k = 1:length(barHPos)-1
            cpanel = uicontrol(f_selection,'Style','text'); cpanel.Position = [startLoc barHPos(end-k) leng itemH];cpanel.String = ['day' int2str(k)];
        end   

        % get all the editable 
        startLoc = 60;
        leng = 60;
        edit_box = {};
        for k = 1:length(barHPos)-1
            cpanel = uicontrol(f_selection,'Style','edit'); cpanel.Position = [startLoc barHPos(end-k) leng itemH];cpanel.String = '1';
            edit_box{k} = cpanel;
        end
        cpanel = uicontrol(f_selection,'Style','pushbutton','CallBack',{@continue_callback,f_selection, edit_box});cpanel.Position = [startLoc buttonHPos leng itemH];cpanel.String = 'next cell';

        % get all the redraw
        startLoc = 130;
        leng = 60;
        redraw_button = {};
        for k = 1:length(barHPos)-1
            cpanel = uicontrol(f_selection,'Style','togglebutton','CallBack',{@redraw_callback, k,j,{localImgList{k}, extraMargin, xlimm, ylimm, xroi, yroi, croi}}); cpanel.Position = [startLoc barHPos(end-k) leng itemH];cpanel.String = 'redraw';
            redraw_button{k} = cpanel;
        end
        cpanel = uicontrol(f_selection,'Style','pushbutton','CallBack',{@reject_callback,f_selection,items});cpanel.Position = [startLoc buttonHPos leng itemH];cpanel.String = 'reject all';
        uiwait(f_selection);
        %tempIsCell = str2double(tempIsCell);
        
        % draw red stuff on the rejected cells
        nop = find(tempIsCell==0);
        figure(figdays);
        for this=1:length(nop)
            k = nop(this);            
            subplot(4,nCol,nCol*4-nCol*2+k);   
            patch(yroi,xroi,'r','EdgeColor','none');
        end
        pause(2);
        %ishere{i}(j,:) = answr;
        ishere{i}(j,:) = tempIsCell;
        close(figdays);
        toc;
        
        tic;
        if checkmode
            pause();
            %close(figdays);
            continue
        end
        
        % Save indiv plane
        toc;
        tic;
        try
            save(['ishere_test_plane' num2str(i-1)'],'ishere');
        catch
            keyboard
            disp('check this!')
        end
        save([suite2ppath sep 'plane' num2str(i-1) '\roi_redrawn_plane' num2str(i-1) '.mat'],'roi_redrawn');
        toc;
    end
end

try
    cd(suite2ppath);
    save('ishere_combined.mat','ishere');
catch
    keyboard
end


end

% GUI functions

function  reject_callback(hObject,eventdata,f_selection,items)
    global tempIsCell 
    tempIsCell = zeros(1,items);
    disp(tempIsCell)
    disp('reject')
    close(f_selection)
end


function  continue_callback(hObject,eventdata,f_selection, edit_box)
    global tempIsCell
    tempIsCell = [];
    for i = 1:length(edit_box)
        tempStr = get(edit_box{i}, 'String');
        tempIsCell(i) = str2double(tempStr);
    end
    disp('continue')
    disp(tempIsCell)
    close(f_selection)
end

function redraw_callback(hObject,eventdata,nDay,nCell,imgData)
    global roi_redrawn
    button_state = get(hObject,'Value');
    if button_state % if triggered from off --> on, redraw roi
        h_redraw = figure;
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
        
        h_roi = imfreehand;
        try
            roiMask = h_roi.createMask;
            [roix, roiy] = find(roiMask);
        
            roix = roix + round(ylimm(1)-extraMargin) - 1;
            roiy = roiy + round(xlimm(1)-extraMargin) - 1;

            roixy = [roix,roiy];
            disp(nDay)
            uiwait(h_redraw);

            tempCellDay = [nCell nDay];tempC = cell(size(roi_redrawn,1),1); tempC(:) = {tempCellDay};
            existCellDay = cellfun(@isequal,roi_redrawn(:,1),tempC,'UniformOutput',true);
            if sum(existCellDay) >= 1 % if the roi has been redrawn, repalce it
                %idx = find(existCellDay==1);
                disp('repeated roi in:')
                disp(roi_redrawn{existCellDay,1})
                roi_redrawn{existCellDay,2} = roixy;
            else % if roi not redrawn, create it
                roi_redrawn = [roi_redrawn;{tempCellDay,roixy}];
            end
        catch
            disp('no roi redrawn')
            set(hObject,'Value',0); 
        end
        
        
    else % if triggered from on --> off, delete redrawn roi
        tempCellDay = [nCell nDay];tempC = cell(size(roi_redrawn,1),1); tempC(:) = {tempCellDay};
        existCellDay = cellfun(@isequal,roi_redrawn(:,1),tempC,'UniformOutput',true);
        if sum(existCellDay) >= 1 % if the roi has been redrawn, delete it
            roi_redrawn(existCellDay,:) = [];
        end
    end
end




