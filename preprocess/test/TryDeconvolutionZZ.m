% mouse = 'cd017';
% suite2ppath = 'H:\celine\cd017\suite2p\';
% 
% signals = cell(nPlanes,1);
% for i=1:nPlanes
%    cd([suite2ppath 'plane' num2str(i-1)]);
%    tc = load([mouse '_TC_plane' num2str(i-1) '.mat']);
%    signals{i} = tc.tempTC;
%    for j=1:nFiles
%        submat = signals{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i));
%        signals{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i),2) = ...
%            (submat-median(submat,2))./median(submat,2);%./repmat(max(submat,[],2)-min(submat,[],2),1,size(submat,2));
%    end
% end

mouse = 'cd036';
global info;

if strcmp(mouse,'cd017')
    suite2ppath = 'I:\celine\cd017\suite2p\';
    h5path = 'V:\LabData4\celine\cd017\h5\'; % h5path = 'H:\celine\cd017\h5\';
    sbxpath = 'V:\LabData4\celine\cd017\';
    behavpath = 'V:\LabData4\celine\cd017\behavior\';
end

if strcmp(mouse,'cd036')
    suite2ppath = 'I:\celine\cd036\suite2p\';
    h5path = 'U:\LabData5\cd036\'; 
    sbxpath = 'U:\LabData5\cd036\';
    behavpath = 'U:\LabData5\cd036\behavior\';
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
nFrames_oneplane = nan(nFiles,nPlanes);
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

% nFrames_oneplane = nan(nFiles,nPlanes);
% nFrames_oneplane(logical(mod(nFrames_add,2)),:) = [round(nFrames_add(logical(mod(nFrames_add,2)))/nPlanes) round(nFrames_add(logical(mod(nFrames_add,2)))/nPlanes)-1];
% nFrames_oneplane(~mod(nFrames_add,2),:) = [nFrames_add(~mod(nFrames_add,2))/nPlanes nFrames_add(~mod(nFrames_add,2))/nPlanes];
nFrames_oneplane = [zeros(1,nPlanes);nFrames_oneplane];

% EXCEPTIONS 'NAMES' FOR CD017
if strcmp(mouse,'cd017'), names{end} = 'cd017_018_000'; end

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
%%

dayss_behav = day1:lastday; 
dayss = sort([tcdays;dayss_behav']);
nDays_behav = length(dayss_behav);
nDays = nDays_behav + length(tcdays);
matrix= [];
list_b = cell(nDays,1);list_f = cell(nDays,1);
previousTrial = 0;
cd(suite2ppath);
for i=1:nDays
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

    if ~ismember(this_day,dayss_behav), continue, end

    toCompare = {[mouse '_' num2str(this_day) 'v']};toCmpr = {};
    for k=1:nSessions, toCmpr = [toCmpr;toCompare]; end      
    ok = cellfun(@strfind,names_behav,toCmpr,'UniformOutput',false); 
    ok2 = nan(nFiles,1);
    for k=1:nSessions, if isempty(ok{k}), continue, end; ok2(k) = ok{k}; end
    n_behav = nansum(ok2);
    list_behav = find(ok2==1); 

    
    for k=1:n_behav
        d = importdata([behav_files(list_behav(k)).folder '\' behav_files(list_behav(k)).name]);                

        if (i==1 && k==1)
            foil = unique([d(d(:,RESPONSE)==FA,TONE);d(d(:,RESPONSE)==CR,TONE)]);
            target = unique([d(d(:,RESPONSE)==H,TONE);d(d(:,RESPONSE)==M,TONE)]);
        end
        
        trialnb = d(:,TRIAL) + previousTrial;
        nTrials = length(trialnb);
        
        contxt = d(:,CONTEXT);
        resp = d(:,RESPONSE);
        toneFrame = d(:,TONEF);
        rewardFrame = d(:,LICKF);
        
        toneFramePlane = [toneFrame toneFrame];
        toneFramePlane(logical(mod(toneFrame,2)),:) = [round(toneFrame(logical(mod(toneFrame,2)))/nPlanes) round(toneFrame(logical(mod(toneFrame,2)))/nPlanes)-1];
        toneFramePlane(~mod(toneFrame,2),:) = [toneFrame(~mod(toneFrame,2))/nPlanes toneFrame(~mod(toneFrame,2))/nPlanes];

        rewardFramePlane = [rewardFrame rewardFrame];
        rewardFramePlane(logical(mod(rewardFrame,2)),:) = [round(rewardFrame(logical(mod(rewardFrame,2)))/nPlanes) round(rewardFrame(logical(mod(rewardFrame,2)))/nPlanes)-1];
        rewardFramePlane(~mod(rewardFrame,2),:) = [rewardFrame(~mod(rewardFrame,2))/nPlanes rewardFrame(~mod(rewardFrame,2))/nPlanes];
        
        rewardFramePlane(rewardFrame==1000000,:) = nan;
        
        if (i==1 && k==4)
%             toneFrame = toneFrame+nFrames_add(list_rec(k));            
            toneFramePlane = toneFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nTrials,1);
            rewardFramePlane = rewardFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nTrials,1);
        else
            if ~(k==1 && i==1)
                toneFramePlane = toneFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nTrials,1);
                rewardFramePlane = rewardFramePlane+repelem(nFrames_oneplane(list_rec(k),:),nTrials,1);
            end
        end            
        matrix = [matrix;i*ones(nTrials,1) list_behav(k)*ones(nTrials,1) trialnb contxt resp toneFramePlane rewardFramePlane];
        
        list_b{i} = list_behav;
        list_f{i} = list_rec;
        previousTrial = previousTrial+nTrials;
    end    
end


%% Column names of the new matrix

DAY = 1;
BLOC = 2;
TRIAL = 3;
CONTEXT = 4;
RESP = 5;
% TONEF = 6;

TONEF1 = 6; % tone frame plane 1
TONEF2 = 7; % tone frame plane 2
REW1 = 8;
REW2 = 9;

% nPlanes = 2;
% tone_fr = matrix(:,TONEF);
% nTones = size(matrix,1);
% if sum(~mod(tone_fr,2))~=nTones
%     start_tone = nan(nTones,nPlanes);
% %     start_tone = [round(tone_fr/nPlanes) round(tone_fr/nPlanes)];
%     start_tone(logical(mod(tone_fr,2)),:) = [round(tone_fr(logical(mod(tone_fr,2)))/nPlanes) round(tone_fr(logical(mod(tone_fr,2)))/nPlanes)-1];
%     start_tone(~mod(tone_fr,2),:) = [tone_fr(~mod(tone_fr,2))/nPlanes tone_fr(~mod(tone_fr,2))/nPlanes];
%     matrix(:,TONEF1) = start_tone(:,1);
%     matrix(:,TONEF2) = start_tone(:,2);
% end

% matrix(:,TONEF1) = matrix(:,TONEF);
% matrix(:,TONEF2) = matrix(:,TONEF);


%% Load cell activity + tracking matrix
cd(suite2ppath);
% ishere = load('ishere_combined.mat');
% ishere = ishere.ishere;
if strcmp(mouse,'cd017')
    cd([suite2ppath 'plane0'])
    ishere1 = load('ishere_plane0.mat');
    ishere1 = ishere1.ishere;
    cd([suite2ppath 'plane1'])
    ishere2 = load('ishere.mat');
    ishere2 = ishere2.ishere;
    ishere = ishere1;
    ishere{2} = ishere2{2};
elseif strcmp(mouse,'cd036')
    cd([suite2ppath 'plane0'])
    ishere1 = load('ishere_plane0.mat');
    ishere1 = ishere1.ishere;
    cd([suite2ppath 'plane1'])
    ishere2 = load('ishere_plane1.mat');
    ishere2 = ishere2.ishere;
    ishere = ishere1;
    ishere{2} = ishere2{2};
end

nbtrialperdays = Accumulate(matrix(:,DAY));
startfrom = dayss_behav(1)+1;
nbtrialperdays = nbtrialperdays(startfrom:end);

signals = cell(nPlanes,1);
signals_baselineremoved = cell(nPlanes,1);
for i=1:nPlanes
    cd([suite2ppath 'plane' num2str(i-1)]);
    tc = load([mouse '_TC_plane' num2str(i-1) '.mat']);
    signals{i} = tc.tempTC; 
    for j=1:nFiles
        submat = signals{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i));  
        %signals{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i),2) = ...
        %    (submat-median(submat,2))./median(submat,2);%./repmat(max(submat,[],2)-min(submat,[],2),1,size(submat,2));
        %temp = (submat-median(submat,2))./median(submat,2);
        %baseline = movmean(temp,200,2);
        %signals_baselineremoved{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i),2) = ...
        %   temp- baseline;
        
        tempDff = (submat-quantile(submat,0.5,2))./quantile(submat,0.5,2);
        %baseline = movmean(tempDff,500,2);
        %tempDff = tempDff - baseline;
        signals{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i),2) = ...
            tempDff;%./repmat(max(submat,[],2)-min(submat,[],2),1,size(submat,2));
        
        tempStd = std(tempDff,0,2);
        signals{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i),3) = ...
            (tempDff - mean(tempDff,2)) ./ tempStd;%./repmat(max(submat,[],2)-min(submat,[],2),1,size(submat,2));
        
    end 
end

%% get rid of baseline
allSessList = list_f;
allSessList{1} = [1;2;3;4];
allSessList{18} = [71;72];
allSessList{19} = [73;74];
signals_nothereremoved = {};
for i=1:nPlanes
    %cd([suite2ppath 'plane' num2str(i-1)]);
    %tc = load([mouse '_TC_plane' num2str(i-1) '.mat']);
    %signals_baselineremoved{i} = signals{i};
    signals_nothereremoved{i} = signals{i}(:,:,2);
    for k = 1:size(signals{i},1)
        try
            goodSessionList = cell2mat(allSessList(logical(ishere{i}(k,:))));
        catch
            disp('ahh')
        end
        for j=1:nFiles
            submat = signals{i}(k,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i),2);
            %baseline = movmean(submat,200,2);
            %figure; for k = 1:4; subplot(4,1,k);plot(submat(k,:));hold on; plot(baseline(k,:));end
            %signals_baselineremoved{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i),2) = ...
            %    submat- baseline;%./repmat(max(submat,[],2)-min(submat,[],2),1,size(submat,2));
            if any(goodSessionList==j)
                signals_nothereremoved{i}(k,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i)) = ...
                   submat;
            else
                signals_nothereremoved{i}(k,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i)) = nan;
            end
        end
    end
end


%%
c_deconvolve_by_trial = {};
s_deconvolve_by_trial = {};
params = {};
%baselineFile = [6 11 16 21 26 31 36 41 46 51 56 61 66 71 76 82];
%baselineFile = [8 13 18 22 26 30 34 38 42 46 50 54 58 62 66 70 71 73];
%cd036
for i=1:nPlanes
    %cd([suite2ppath 'plane' num2str(i-1)]);
    %tc = load([mouse '_TC_plane' num2str(i-1) '.mat']);
    c_deconvolve_by_trial{i} = zeros(size(signals{i},1),size(signals{i},2));
    s_deconvolve_by_trial{i} = zeros(size(signals{i},1),size(signals{i},2));
    for j=1:nFiles
       submat = signals{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i),2);
       submat(:,1) = 0; % correct of the first frame
       %baseline = movmean(submat,200,2);
       tic
       for k = 1:size(submat,1)
            
            %[c, s, options] = deconvolveCa(submat(k,:),'foopsi','ar1', 0.9381, 'optimize_b');
            try
                [c, s, options] = deconvolveCa(submat(k,:),'constrained','ar1', 0.9381, 'optimize_b','optimize_pars');% time x one cell
                c_deconvolve_by_trial{i}(k,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i)) = c;
                s_deconvolve_by_trial{i}(k,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i)) = s;
                params{i,j,k} = options;
            catch
                c_deconvolve_by_trial{i}(k,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i)) = nan;
                s_deconvolve_by_trial{i}(k,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i)) = nan;
                params{i,j,k} = nan;
                disp([k j])
            end
       end
       toc;
    end
end
save('cd017_TC_deconvolve_on_both_plane_by_trial_params_thresholded_fix_decay_fit_b.mat','params','-v7.3')
save('cd017_TC_deconvolve_on_both_plane_by_trial_s_thresholded_fix_decay_fit_b.mat','s_deconvolve_by_trial','-v7.3')
save('cd017_TC_deconvolve_on_both_plane_by_trial_c_thresholded_fix_decay_fit_b.mat','c_deconvolve_by_trial','-v7.3')
%% save restricted traces
c_deconvolve_all_session = {};
s_deconvolve_all_session = {};
restrict_traces = {};
for p=1:nPlanes
    tic;
    signals_thisplane = signals_baselineremoved{p};
    nCells = size(signals_thisplane,1); 
%     for i=86:nCells
    for i=1:size(signals_thisplane,1)
        tic;
        trace = signals_thisplane(i,:,2);
        %figure;plot(trace);title(['cell ' num2str(i)]);
        %load([suite2ppath '\plane' int2str(p-1) '\ishere_plane'  int2str(p-1) '.mat'])
        %here = ishere{p}(i,:);
      
        nFrames = size(trace,2);
        restrict = 180000:280000; 
        %figure;plot(trace(restrict));title(['cell ' num2str(i)]);
        [c, s, options] = deconvolveCa(trace(restrict),'foopsi','ar1', 0.9381, 'optimize_b'); % time x one cell
        restrict_traces = [restrict_traces trace(restrict)];
        
        % outputs:
        %   c: T x 1 vector, denoised trace
        %   s: T x 1 vector, deconvolved signal
        %   b: fluorescence baseline
        %   kernel: struct variable containing the parameters for the selected
        %       convolution model
        %   lambda: Optimal Lagrange multiplier for noise constraint under L1 penalty
        %     """olves the noise constrained sparse nonnegat

        xlimm = [1 3000];
        
        c_deconvolve_all_session = [c_deconvolve_all_session c];
        s_deconvolve_all_session = [s_deconvolve_all_session s];
        
%         figure;
%         
%         subplot(3,1,1); hold on;
%         plot(trace(restrict));
%         title('Raw trace (median centered)');
%         xlabel('frames'); ylabel('ddf');
%         xlim(xlimm)
%         
%         subplot(3,1,2); hold on;
%         plot(c);
%         title('Denoised trace');
%         xlabel('frames'); ylabel('ddf');
%         xlim(xlimm)
%         
%         subplot(3,1,3); hold on;
%         plot(s);
%         title('Deconvolved trace');
%         xlabel('frames'); ylabel('spike rate');
%         xlim(xlimm)
%         pause;
        toc;
    end
    toc;
end

%%
%save('cd017_TC_deconvolve_on_both_plane_all_sess_remove_neuropil_s_foopsi_fix_decay_fit_b.mat','s_deconvolve_all_session','-v7.3')
%save('cd017_TC_deconvolve_on_both_plane_all_sess_remove_neuropil_c_foopsi_fix_decay_fit_b.mat','c_deconvolve_all_session','-v7.3')
%save('cd017_TC_deconvolve_on_both_plane_remove_neuropil_trace.mat','c_deconvolve_all_session','-v7.3')
%save('cd017_TC_deconvolve_on_both_plane_by_trial_remove_neuropil_s_foopsi_fix_decay_fit_b.mat','s_deconvolve_by_trial','-v7.3')
%save('cd017_TC_deconvolve_on_both_plane_by_trial_remove_neuropil_c_foopsi_fix_decay_fit_b.mat','c_deconvolve_by_trial','-v7.3')
dff{1} = signals{1}(:,:,2);dff{2} = signals{2}(:,:,2);
save('cd017_TC_dff_remove_neuropil','dff','-v7.3')
%%

c_deconvolve_all_session = zeros(length(restrict_traces), length(restrict_traces{1}));
s_deconvolve_all_session = zeros(length(restrict_traces), length(restrict_traces{1}));


for i=1:length(restrict_traces)
    tic;
    
    restrict_traces = [restrict_traces trace(restrict)];
    [c, s, b, g, active_set] = foopsi_oasisAR1(restrict_traces{i}, 0.9381, 0.1, true,...
    false, [], 100);
    % outputs:
    %   c: T x 1 vector, denoised trace
    %   s: T x 1 vector, deconvolved signal
    %   b: fluorescence baseline
    %   kernel: struct variable containing the parameters for the selected
    %       convolution model
    %   lambda: Optimal Lagrange multiplier for noise constraint under L1 penalty
    %     """olves the noise constrained sparse nonnegat

    xlimm = [1 3000];

    %c_deconvolve_all_session = [c_deconvolve_all_session c];
    %s_deconvolve_all_session = [s_deconvolve_all_session s];

%         figure;
%         
%         subplot(3,1,1); hold on;
%         plot(trace(restrict));
%         title('Raw trace (median centered)');
%         xlabel('frames'); ylabel('ddf');
%         xlim(xlimm)
%         
%         subplot(3,1,2); hold on;
%         plot(c);
%         title('Denoised trace');
%         xlabel('frames'); ylabel('ddf');
%         xlim(xlimm)
%         
%         subplot(3,1,3); hold on;
%         plot(s);
%         title('Deconvolved trace');
%         xlabel('frames'); ylabel('spike rate');
%         xlim(xlimm)
%         pause;
    toc;
end
