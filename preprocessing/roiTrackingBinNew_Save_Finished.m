function roiTrackingBinNew_Save_Finished(mouse)
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
%load(['config\mouse\' mouse '_ImagingMetaData.mat']);
suite2ppath = 'H:\celine\cd036\suite2p';

configPath = 'E:\KishoreLab\Shared\Matlab\preprocessing-newpipeline-draft\preprocessing\config';
mouseDataPath = [configPath sep 'mouse' sep mouse];
mkdir(mouseDataPath)
mkdir([mouseDataPath sep 'checkCell' sep])

cd(suite2ppath);

imgData = func_loadMouseConfig(mouse,'root','E:\KishoreLab\Shared\Matlab\preprocessing-newpipeline-draft\preprocessing\');
allDay = 0:18;
behavDay = 1:16;

nDays = length(allDay);
avg_day = cell(nDays,2);
nPlanes = 2;
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
        
        
    end
end

warning('off')

for i=1:nPlanes
    
    cd([suite2ppath sep 'plane' num2str(i-1)]);

    ishere = load([suite2ppath sep 'plane' num2str(i-1) sep  'ishere_plane' num2str(i-1) '.mat']);
    ishere = ishere.ishere{i};
    extraMargin = 20;
    nCells = size(ishere,1);
    
    roiName = [mouse '_roi' int2str(i-1) '.zip'];
    rois = ReadImageJROI(roiName); %read imagej rois
%%    
    
    for j=1:nCells  
        disp(j)
        tic;
        % Plot this ROI in the field across days    
        yroi = rois{1,j}.mnCoordinates(:,1); %freehand rois have the outlines in x-y coordinates
        xroi = rois{1,j}.mnCoordinates(:,2);
        croi = [mean(minmax(xroi)) mean(minmax(yroi))];
        ylimm = [croi(1)-20 croi(1)+20];
        xlimm = [croi(2)-20 croi(2)+20];
        nCol = round(nDays/2);
       
        screensize = get( groot, 'Screensize' );
        saveFig = figure('Position',[screensize(3)*0.05,screensize(4)*0.4,...
            screensize(3)*0.9,screensize(4)*0.5],'visible','on');
        
        tempCol = ceil(nDays/3);

        for d = 1:nDays
            subplot_tight(3,tempCol,d,0.02)
            hold on;
            localImg = uint16(avg_day{d,i}(ylimm(1)-extraMargin:ylimm(2)+extraMargin, xlimm(1)-extraMargin:xlimm(2)+extraMargin));
            localImgList{d} = localImg;
            imagesc(imadjust(adapthisteq(localImgList{d}, 'NBins', 256),[0 1],[0 1],0.2)); 
            xlim([extraMargin+1 xlimm(2)-xlimm(1)+extraMargin+1]);ylim([extraMargin+1 ylimm(2)-ylimm(1)+extraMargin+1]);
            title(['D' num2str(d)]);
            colormap gray;
            
            if ishere(j,d)==1
                pat = patch(yroi-round(xlimm(1)-extraMargin)+1, xroi-round(ylimm(1)-extraMargin)+1,'g','EdgeColor','None');
                pat.FaceAlpha = 0.2;
            elseif ishere(j,d)==0
                pat = patch(yroi-round(xlimm(1)-extraMargin)+1, xroi-round(ylimm(1)-extraMargin)+1,'r','EdgeColor','None');
                pat.FaceAlpha = 0.2;
            elseif ishere(j,d)==2
                pat = patch(yroi-round(xlimm(1)-extraMargin)+1, xroi-round(ylimm(1)-extraMargin)+1,'b','EdgeColor','None');
                pat.FaceAlpha = 0.2;
            end
            set(gca,'YTickLabel',[]);
            set(gca,'XTickLabel',[]);
            
        end
        %pause(2);
        saveas(saveFig,[mouseDataPath sep 'checkCell' sep 'check_cell_plane' num2str(i-1) 'cell' int2str(j) '.png']);
        close(saveFig)   
    end
end
warning('on')


