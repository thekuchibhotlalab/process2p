function RoiTrackingBinZZ(mouse)

global info;

if strcmp(mouse,'cd017')
    suite2ppath = 'H:\celine\cd017\suite2p\';
    h5path = 'W:\LabData4\celine\cd017\h5';
    sbxpath = 'W:\LabData4\celine\cd017\';
    behavpath = 'W:\LabData4\celine\cd017\behavior\';
end

if strcmp(mouse,'cd036')
    suite2ppath = 'H:\celine\cd036\suite2p\';
    h5path = 'H:\celine\cd036\';
    sbxpath = 'V:\LabData5\cd036\';
    behavpath = 'V:\LabData5\cd036\behavior\';
end

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

% Get the number of frames of all the files process in suite2p
cd(h5path);
files = dir('*.h5');
nFiles = length(files);
names = cell(nFiles,1);

for i=1:nFiles, names{i} = files(i).name(1:end-3); end

cd(sbxpath);
nFrames = nan(nFiles,1);nFrames_add = nan(nFiles,1);
sbxread(names{1},1,1);
nPlanes = info.otparam(3);
for i=1:nFiles
    sbxread(names{i},1,1);
    nFrames(i) = info.max_idx;    
    if i>1
        nFrames_add(i) = nFrames(i)+nFrames_add(i-1); 
    else
        nFrames_add(i) = nFrames(i); 
    end
    
    if mod(nFrames(i),2)
        if i>1
            nFrames_oneplane(i,:) = [round(nFrames(i)/nPlanes)+nFrames_oneplane(i-1,1) round(nFrames(i)/nPlanes)-1+nFrames_oneplane(i-1,2)];
        else
            nFrames_oneplane(i,:) = [round(nFrames(i)/nPlanes) round(nFrames(i)/nPlanes)-1];
        end
    else
        nFrames_oneplane(i,:) = [nFrames(i)/nPlanes+nFrames_oneplane(i-1,1) nFrames(i)/nPlanes+nFrames_oneplane(i-1,2)];
    end        
end
nFrames_oneplane = [zeros(1,nPlanes);nFrames_oneplane];% Here, the nb of frames / plane IS cumulative

% EXCEPTIONS 'NAMES' FOR CD017
if strcmp(mouse,'cd017'), names{end} = 'cd017_018_000'; end

%% 

% Load behavioral files and place them in the right order
cd(behavpath);
behav_files = dir('*_*v*.txt');
f=1;
while ~strcmp(behav_files(f).name,[mouse '_1v1.txt'])
    f=f+1;
end
behav_files = circshift(behav_files,-(f-1));
nSessions = size(behav_files,1);
names_behav = cell(nSessions,1);
for i=1:nSessions, names_behav{i} = behav_files(i).name; end

% Load the mean image saved for each session (from binary file)
% define day 1 and nb od days of recording, and find target and foil frames

% find which recording session correspond to the first day of behavior
% idea: check file date
cd(sbxpath);
sbxfiles = dir('*.sbx');
day1_behavior = behav_files(1).date(1:end-8);
good1day='';a=1;
while ~strcmp(good1day,day1_behavior)
    good1day = sbxfiles(a).date(1:end-8);
    a=a+1;
end
day1 = str2double(names{a-1}(7:9));
tcdays = day1-1;
% day1 = str2double(names{1}(7:9)); % of behavior
lastday_behavior = behav_files(end).date(1:end-8);
goodlastday='';a=0;
while ~strcmp(goodlastday,lastday_behavior)
    goodlastday = sbxfiles(end-a).date(1:end-8);
    a=a+1;
end
lastday = str2double(names{end-a-1}(7:9));
lastdaytc = str2double(names{end}(7:9));
if lastdaytc>lastday
    tcdays = [tcdays;(lastday+1:lastdaytc)'];
end
% lastday = str2double(names{end}(7:9));
dayss_behav = day1:lastday; 
dayss = sort([tcdays;dayss_behav']);
nDays_behav = length(dayss_behav);
nDays = nDays_behav + length(tcdays);
avg_day = cell(nDays,2);
foil_fr = [];target_fr = [];
list_b = cell(nDays,1);list_f = cell(nDays,1);
cd(suite2ppath);
for j=1:nPlanes
    cd([suite2ppath 'plane' num2str(j-1)])
    avg_session = load([mouse '_MeanImgPerSessions_Plane' num2str(j) '.mat']);
    avg_session = avg_session.tosave;
    for i=1:nDays
%         this_day = dayss_behav(i);
        this_day = dayss(i);
        if length(num2str(this_day))==2, toCompare = {[mouse '_0' num2str(this_day)]};
        else, toCompare = {[mouse '_00' num2str(this_day)]}; end
        toCmpr = {};
        for k=1:nFiles, toCmpr = [toCmpr;toCompare]; end       
        ok = cellfun(@strfind,names,toCmpr,'UniformOutput',false); 
        ok2 = nan(nFiles,1);
        for k=1:nFiles, if isempty(ok{k}), continue, end; ok2(k) = ok{k}; end
        n_rec = nansum(ok2);
        list_rec = find(ok2==1);
        disp(['D' num2str(i) ', files ' num2str(list_rec')]);

        m = [];
        for k=1:n_rec
            m(:,:,k) = avg_session{list_rec(k)};
        end
        avg_day{i,j} = mean(m,3);
%         figure;imagesc(mean(m,3)); colormap gray;
%     end
% end
        if ~ismember(this_day,dayss_behav), continue, end
        
        if j==1 
            toCompare = {[mouse '_' num2str(this_day) 'v']};toCmpr = {};
            for k=1:nSessions, toCmpr = [toCmpr;toCompare]; end      
            ok = cellfun(@strfind,names_behav,toCmpr,'UniformOutput',false); 
            ok2 = nan(nFiles,1);
            for k=1:nSessions, if isempty(ok{k}), continue, end; ok2(k) = ok{k}; end
            n_behav = nansum(ok2);
            list_behav = find(ok2==1); 

            for k=1:n_behav
                d = importdata([behav_files(list_behav(k)).folder '\' behav_files(list_behav(k)).name]);                
                nTrials = size(d,1);
                
                foil = unique([d(d(:,RESPONSE)==FA,TONE);d(d(:,RESPONSE)==CR,TONE)]);
                target = unique([d(d(:,RESPONSE)==H,TONE);d(d(:,RESPONSE)==M,TONE)]);
                        
                foil_frames = d(ismember(d(:,TONE),foil),TONEF);
                target_frames = d(ismember(d(:,TONE),target),TONEF);
                nFoil = size(foil_frames,1);
                nTarget = size(target_frames,1);
                
                foilFramePlane = [foil_frames foil_frames];
                foilFramePlane(logical(mod(foil_frames,2)),:) = [round(foil_frames(logical(mod(foil_frames,2)))/nPlanes) round(foil_frames(logical(mod(foil_frames,2)))/nPlanes)-1];
                foilFramePlane(~mod(foil_frames,2),:) = [foil_frames(~mod(foil_frames,2))/nPlanes foil_frames(~mod(foil_frames,2))/nPlanes];
        
                targetFramePlane = [target_frames target_frames];
                targetFramePlane(logical(mod(target_frames,2)),:) = [round(target_frames(logical(mod(target_frames,2)))/nPlanes) round(target_frames(logical(mod(target_frames,2)))/nPlanes)-1];
                targetFramePlane(~mod(target_frames,2),:) = [target_frames(~mod(target_frames,2))/nPlanes target_frames(~mod(target_frames,2))/nPlanes];
        
                if (i==1 && k==4)
                    foilFramePlane = foilFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nFoil,1);
                    targetFramePlane = targetFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nTarget,1);
                else
                    if ~(k==1 && i==1)
                        foilFramePlane = foilFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nFoil,1);
                        targetFramePlane = targetFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nTarget,1);
                    end
                end
                
                foil_fr = [foil_fr;foilFramePlane];
                target_fr = [target_fr;targetFramePlane];
                
%                 if (i==1 && k==4)
%                     disp(['Day ' num2str(i) ', behav ' num2str(list_behav(k)) ', file ' num2str(num2str(list_rec(k)))]);
%                     foil_frames = d(ismember(d(:,TONE),foil),TONEF);
%                     foil_frames = foil_frames+nFrames_add(list_rec(k));
%                     foil_fr = [foil_fr;foil_frames];
%                     target_frames = d(ismember(d(:,TONE),target),TONEF);
%                     target_frames = target_frames+nFrames_add(list_rec(k));
%                     target_fr = [target_fr;target_frames];
%                     disp(['1fr = ' num2str(target_frames(1)) ' add = ' num2str(nFrames_add(list_rec(k))) ', first-last target_frames = ' num2str(target_frames(1)) '-' num2str(target_frames(end))]);
%                     disp('');            
%                 else
%                     disp(['Day ' num2str(i) ', behav ' num2str(list_behav(k)) ', file ' num2str(num2str(list_rec(k)-1))]);
%                     foil_frames = d(ismember(d(:,TONE),foil),TONEF);
%                     if ~(k==1 && i==1), foil_frames = foil_frames+nFrames_add(list_rec(k)-1); end
%                     foil_fr = [foil_fr;foil_frames];
%                     target_frames = d(ismember(d(:,TONE),target),TONEF);
%                     if ~(k==1 && i==1), target_frames = target_frames+nFrames_add(list_rec(k)-1); end
%                     target_fr = [target_fr;target_frames];
%                     if ~(k==1 && i==1)
%                         disp(['1fr = ' num2str(target_frames(1)) ' add = ' num2str(nFrames_add(list_rec(k)-1)) ', first-last target_frames = ' num2str(target_frames(1)) '-' num2str(target_frames(end))]);
%                     end
%                     disp('');
%                 end
                list_b{i} = list_behav;
                list_f{i} = list_rec;
            end  
        end    
    end
end
% checktargetandfoil = false;
% if checktargetandfoil
%     figure;subplot(2,2,1);plot(foil_fr,'.');
%     subplot(2,2,2);plot(target_fr,'.');
% end

%%

% Check that tone onset frame is always even, otherwise correct frame
% number according to each plane

% nFoils = length(foil_fr);nTargets = length(target_fr);
% if (sum(~mod(foil_fr,2))~=nFoils) || (sum(~mod(target_fr,2))~=nTargets)
%     start_foil = nan(nFoils,nPlanes);
%     start_foil(logical(mod(foil_fr,2)),:) = [round(foil_fr(logical(mod(foil_fr,2)))/nPlanes) round(foil_fr(logical(mod(foil_fr,2)))/nPlanes)-1];
%     start_foil(~mod(foil_fr,2),:) = [foil_fr(~mod(foil_fr,2))/nPlanes foil_fr(~mod(foil_fr,2))/nPlanes];
%     start_foil = start_foil - pretone*round(acq); % - pretone sec before tone onset
%     
%     start_target = nan(nTargets,nPlanes);
%     start_target(logical(mod(target_fr,2)),:) = [round(target_fr(logical(mod(target_fr,2)))/nPlanes) round(target_fr(logical(mod(target_fr,2)))/nPlanes)-1];
%     start_target(~mod(target_fr,2),:) = [target_fr(~mod(target_fr,2))/nPlanes target_fr(~mod(target_fr,2))/nPlanes];
%     start_target = start_target - pretone*round(acq); % - pretone sec before tone onset
%     
%     nframes_psth = pretone*round(acq) + posttone*round(acq);
% else
%     start_foil = round(foil_fr/nPlanes) - pretone*round(acq); % - pretone sec before tone onset
%     start_target = round(target_fr/nPlanes) - pretone*round(acq); 
%     nframes_psth = pretone*round(acq) + posttone*round(acq);        
% end

% start_foil = round(foil_fr/nPlanes) - pretone*round(acq); % - pretone sec before tone onset
% start_target = round(target_fr/nPlanes) - pretone*round(acq); 
% nframes_psth = pretone*round(acq) + posttone*round(acq); 

start_foil = foil_fr - pretone*round(acq); % - pretone sec before tone onset
start_target = target_fr - pretone*round(acq); 
nframes_psth = pretone*round(acq) + posttone*round(acq); 

% Compute everything possible before plotting for speed purpose
mat = cell(nPlanes,1);
for i=1:nPlanes
    cd([suite2ppath 'plane' num2str(i-1)]);
    tc = load([mouse '_TC_plane' num2str(i-1) '.mat']);
    mat{i} = tc.tempTC; 
    for j=1:nFiles
        submat = mat{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i));
        mat{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i),2) = ...
            (submat-median(submat,2))./median(submat,2);%./repmat(max(submat,[],2)-min(submat,[],2),1,size(submat,2));
    end 
end

% Identify baseline and tuning curve sessions
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
% Plot normalized fluo, ROI in the field and PSTH across days
ishere = cell(nPlanes,1);
checkmode = false;
warning('off')

for i=2:nPlanes
    cd([suite2ppath 'plane' num2str(i-1)]);
    m = mat{i}; % 1D=fluo; 2D=df/f
    zfluo = squeeze(m(:,:,2));
    nCells = size(zfluo,1);
    if ~exist([suite2ppath 'plane' num2str(i-1) '\ishere_plane' num2str(i) '.mat'],'file')
        ishere{i} = nan(nCells,nDays);
        c = 1;
    else
        ishere = load([suite2ppath 'plane' num2str(i-1) '\ishere_plane' num2str(i) '.mat']);
        ishere2 = ishere.answr;
        ishere = cell(nPlanes,1);
        ishere{i} = ishere2;
        c = size(ishere{i},2);
    end
    roiName = [mouse '_roi' int2str(i-1) '.zip'];
    rois = ReadImageJROI(roiName); %read imagej rois
    data = load('Fall.mat');
%%    
    for j=902:nCells  
        figdays = figure;         
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);% Enlarge figure to full screen.
        
        % Plot norm traces
        subplot(4,1,1);hold on;title(['Cell ' num2str(j) ', raw traces']);
        plot(zfluo(j,:));        
        for k=1:nBaseline
            plot(intervals_baseline{i}(k,1)+1:intervals_baseline{i}(k,2),zfluo(j,intervals_baseline{i}(k,1)+1:intervals_baseline{i}(k,2)),'g');% color the baseline traces
        end
        PlotHVLines(nFrames_oneplane(2:end,i),'v','color',[0.8 0.8 0.8]);
        xlim([0 size(zfluo(j,:),2)]);
        for k=1:nTuning
            plot(intervals_tc{i}(k,1)+1:intervals_tc{i}(k,2),zfluo(j,intervals_tc{i}(k,1)+1:intervals_tc{i}(k,2)),'color',[0.8 0.8 0.8]);
        end        
        
        % Plot tone-evoked resp        
        [~,whichsessions_foil] = InIntervals(start_foil(:,i),[nFrames_oneplane(2:end-1,i) nFrames_oneplane(3:end,i)]);
        [~,whichsessions_target] = InIntervals(start_target(:,i),[nFrames_oneplane(2:end-1,i) nFrames_oneplane(3:end,i)]);
        for k=1:nDays
            if ~ismember(k-1,dayss_behav), continue, end
            subplot(4,nDays*2,nDays*2+(k-1)*2+1);hold on;
            ok_foil = ismember(whichsessions_foil,list_b{k});
            l = sum(ok_foil);
            matidx = repelem(0:nframes_psth-1,l,1);
            idx = (matidx+repelem(start_foil(ok_foil,i),1,nframes_psth))';
            m_foil = zfluo(j,idx(:));
            traces_foil = reshape(m_foil,nframes_psth,l);
            plot(traces_foil,'color',[0.8 0.8 0.8]);plot(mean(traces_foil,2),'k','linewidth',2);
            title(['D' num2str(k)]);
            if k-1==dayss_behav(1), ylimm1 = ylim; else, ylim(ylimm1); end
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
            imagesc(imadjust(adapthisteq(localImg, 'NBins', 256),[0 1],[0 1],0.2)); 
            xlim([extraMargin+1 xlimm(2)-xlimm(1)+extraMargin+1]);ylim([extraMargin+1 ylimm(2)-ylimm(1)+extraMargin+1]);
%             axis off;
            title(['D' num2str(k)]);
            colormap gray;
            plot(extraMargin+1+croi(2)-xlimm(1),extraMargin+1+croi(1)-ylimm(1),'r+');
        end    
        
        if checkmode
            pause();
            close(figdays);
            continue
        end
        
        % Open GUI for selecting days 
        listDays = {};
        for k=1:nDays, listDays{k} = ['Day ' num2str(k)]; end
        titleWindow = '';
        dims = [1 20];  
        definput = {};
        for k=1:nDays, definput{k} = [num2str(1)]; end
        answer = inputdlg(listDays,titleWindow,dims,definput);
        answr = str2double(answer)';
        
        nop = find(answr==0);
        figure(figdays);
        for this=1:length(nop)
            k = nop(this);            
            subplot(4,nCol,nCol*4-nCol*2+k);   
            patch(yroi,xroi,'r','EdgeColor','none');
        end
        pause(2);
        ishere{i}(j,:) = answr;
        close(figdays);
        
        % Save indiv plane
        try
            if j==1
                if exist(['ishere_plane' num2str(i) '.mat'],'file')
                    continue
                end
                save(['ishere_plane' num2str(i) '.mat'],'answr');   
            else
                save(['ishere_plane' num2str(i) '.mat'],'answr','-append');   
            end
        catch
            keyboard
        end
    end
end


% Save the matrix of both planes
try
    cd(suite2ppath);
    save('ishere_combined.mat','ishere');
catch
    keyboard
end