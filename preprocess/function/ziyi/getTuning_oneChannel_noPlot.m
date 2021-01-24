function getTuning_oneChannel_noPlot(TC,savePath,funcParam)

%---------DECLARE SOME PARAMETERS-----------
global nPlanes
sep = '\';
sumTC = sum(TC(1:end-1,:),1); TC(:,sumTC==0) = 5000 + randi(100,size(TC,1),sum(sumTC==0));
nNeuron = size(TC,2);
eval(funcParam);
if ~exist('tonePeakLim','var'); tonePeakLim = 0.66;end
%---------SHIFT THE TC FOR PRETONE PERIOD-----------
TC_original = TC; % keep a copy of original TC
TC = circshift(TC,1-toneOnset,1); % shift 1 on 1st axix so that tone is on frame 1 instead of 0
%---------GAUSSIAN FILTER TO SMOOTH THE TRACES-----------
if smoothWindow ~= 0; TC = smoothdata(TC,1,'gaussian',smoothWindow); end
if ~exist('testByPeak','var'); testByPeak = true; end
TCpretone = circshift(TC,pretoneFrames,1);
%---------COMPUTE DFF-----------
%baseline = prctile(TC,50,1);
%baseline = mean(TC,1);
%TC = TC ./ repmat(baseline,[size(TC,1) 1]);
%TCpretone = TCpretone ./ repmat(baseline,[size(TC,1) 1]);
%---------RESHAPE TC AND PRETONE-----------
TC=reshape(TC,nFramesPerTone,nTones,nTrials,nNeuron); 
TCpretone = reshape(TCpretone,nFramesPerTone,nTones,nTrials,nNeuron);
%---------REORDER TC TO ALIGN WITH TONE FREQUENCY ORDER-----------
toneorder = [45254.834 8000 13454.34264 4756.82846 5656.854249,...
    22627.417 64000 53817.37058 4000 9513.65692,...
    16000 6727.171322 19027.31384 26908.6852 32000,...
    11313.7085 38054.62768];
toneindex = [9;4;5;12;2;10;16;3;11;13;6;14;15;17;1;8;7];
TC_reorder=zeros(size(TC));
TCpretone_reorder = zeros(size(TCpretone));
for x=1:nTones
    index=toneindex(x);%+1; %the "+1" is because the first 20 frames are 38 kHz that wraps around
    TC_reorder(:,x,:,:)=TC(:,index,:,:);
    TCpretone_reorder(:,x,:,:)=TCpretone(:,index,:,:);
end
%---------PLOT ALL TONE EVOKED-----------
trialMean = squeeze(nanmean(TCpretone_reorder,3));
trialMeanTrialTC = reshape(trialMean,[nFramesPerTrial,nNeuron]);
trialSEM = squeeze(nanstd(TCpretone_reorder,0,3)) ./sqrt(nTrials);
trialMedian = squeeze(nanmedian(TCpretone_reorder,3));
trialMedianTrialTC = reshape(trialMedian,[nFramesPerTrial,nNeuron]);
toneMean = squeeze(nanmean(nanmean(TCpretone_reorder,2),3));
%---------PLOT ALL TONE EVOKED-----------
frameAxis = pretoneFrames:20:nFramesPerTone;
frameLabel = cellfun(@num2str,num2cell(frameAxis-pretoneFrames),'UniformOutput',false);
psthFig = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.04, 0.85, 0.96]);% Enlarge figure to full screen.
for i = 1:nTones
    subplot(3,6,i)
    imagesc(squeeze(trialMean(:,i,:))');
    caxis([prctile(trialMean(:),5) prctile(trialMean(:),95)])
    title([num2str(round(toneorder(toneindex(i))),'%d') 'HZ'])
    if i>12
        xlabel('frames')
        frameAxis = pretoneFrames:10:nFramesPerTone;
        frameLabel = cellfun(@num2str,num2cell(frameAxis-pretoneFrames),'UniformOutput',false);
        xticks(frameAxis)
        xticklabels(frameLabel)
    else
        xticklabels([])
    end
    if mod(i,6) == 1
        ylabel('neurons') 
    else
        yticklabels([])
    end
end
subplot(3,6,18)
imagesc(toneMean')
xticks(frameAxis)
xticklabels(frameLabel)
xlabel('frames')
yticklabels([])
caxis([prctile(trialMean(:),5) prctile(trialMean(:),95)])
title('all tone')
saveas(gcf,[ savePath ...
        '/populationTonePSTH.png']);

%---------PLOT PEAK FRAME-----------
trialMean = squeeze(nanmean(TCpretone_reorder,3));

for i = 1:nTones
    [maxValue,tempPeakIdx] = max(squeeze(trialMean(pretoneFrames+1:pretoneFrames+ceil(frameRate*tonePeakLim),i,:)),[],1); 
    peakActAvg(:,i) = maxValue; peakFrames(:,i) =  tempPeakIdx; 
end

[maxValue,peakIndex] = max(toneMean(pretoneFrames+1:pretoneFrames+ceil(frameRate*tonePeakLim),:),[],1);
%[maxValue,peakIndex] = max(toneMean(pretoneFrames+1:pretoneFrames+1+frameRate,:),[],1);
latFig = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.04, 0.45, 0.56]);% Enlarge figure to full screen.
subplot(2,2,1)
histogram(peakIndex)
xlabel('peak frame')
ylabel('frequency')
title('peak frame of all neurons')
subplot(2,2,2)
scatter(peakIndex,maxValue,15,'filled'); hold on; xlimm = xlim; 
plot(xlimm, [prctile(maxValue,50) prctile(maxValue,50)],'--','Color',[0.8 0.8 0.8]);
xlabel('peak frame')
ylabel('peak amplitude')
title('peak frame & amplitude')
subplot(2,2,3)
[~,sortIndex] = sort(peakIndex);
imagesc(toneMean(:,sortIndex)')
xticks(frameAxis)
xticklabels(frameLabel)
ylabel('roi')
title('average dff of all tones')
caxis([prctile(toneMean(:),5) prctile(toneMean(:),95)])
subplot(2,2,4)
plot(toneMean,'Color',[0.8 0.8 0.8]);hold on;
plot(mean(toneMean,2))
xticks(frameAxis)
xticklabels(frameLabel)
ylabel('dff')
title('population average dff')
xlim([0 nFramesPerTone])
saveas(gcf,[ savePath ...
        '/populationTonePeak.png']);

%---------GET ANOVA ACTIVITY-----------
peakANOVA = nan((nTrials)*nTones,nTones+1,nNeuron);
peakANOVACorr = nan((nTrials)*nTones,nTones+1,nNeuron);
pretoneMean = (nanmean(TCpretone_reorder(1:pretoneFrames,:,:,:)));
TCpretone_reorderCorr = TCpretone_reorder - repmat(pretoneMean, [nFramesPerTone, 1, 1, 1]);

if testByPeak
    for i = 1:nTones
        for j = 1:nNeuron
            % Use the peak activity
            % Although most times we are taking 1 frame, but use nanmean here to allow multiple frames
            peakANOVA(((i-1)*(nTrials)+1):(i*(nTrials)),i,j) = ...
                squeeze(nanmean(TCpretone_reorder(pretoneFrames+peakFrames(j,i),i,:,j),1));
            peakANOVA(((i-1)*(nTrials)+1):(i*(nTrials)),end,j) = ...
                squeeze(nanmean(TCpretone_reorder(pretoneFrames,i,:,j),1));

            peakANOVACorr(((i-1)*(nTrials)+1):(i*(nTrials)),i,j) = ...
                squeeze(nanmean(TCpretone_reorderCorr(pretoneFrames+peakFrames(j,i),i,:,j),1));
            peakANOVACorr(((i-1)*(nTrials)+1):(i*(nTrials)),end,j) = ...
                squeeze(nanmean(TCpretone_reorderCorr(pretoneFrames,i,:,j),1));
        end
    end
else
    for i = 1:nTones
        for j = 1:nNeuron
        % Use a window to average activity
        % Use nanmean here to allow multiple frames
        peakANOVA(((i-1)*(nTrials)+1):(i*(nTrials)),i,j) = ...
            squeeze(nanmean(TCpretone_reorder(pretoneFrames+(1:floor(tonePeakLim*frameRate)),i,:,j),1));
        peakANOVA(((i-1)*(nTrials)+1):(i*(nTrials)),end,j) = ...
            squeeze(nanmean(TCpretone_reorder(pretoneFrames+(-ceil(tonePeakLim*frameRate):-1),i,:,j),1));
            
        peakANOVACorr(((i-1)*(nTrials)+1):(i*(nTrials)),i,j) = ...
            squeeze(nanmean(TCpretone_reorderCorr(pretoneFrames+(1:floor(tonePeakLim*frameRate)),i,:,j),1));
        peakANOVACorr(((i-1)*(nTrials)+1):(i*(nTrials)),end,j) = ...
            squeeze(nanmean(TCpretone_reorderCorr(pretoneFrames+(-ceil(tonePeakLim*frameRate):-1),i,:,j),1));
        end
    end
    
end


for i = 1:nNeuron
    tempToneAct = zeros(1,nTones);
    for j = 1:nTones
        if testByPeak
            tempToneAct(j) = trialMedian(pretoneFrames+peakFrames(i,j),j,i);
        else
            tempToneAct(j) = nanmean(trialMedian(pretoneFrames+(1:floor(tonePeakLim*frameRate)),j,i),1);
        end
    end
    [~, tempTuningPeak] = max(tempToneAct);
    tuningPeak(i) = tempTuningPeak;
end
%---------DECLARE ALL MATRICES FOR SIGNIFICANT TESTS-----------
pairedTestTrials = startTrial:nTrials; % trial 1 does not have pair for baseline
% paired ttest for peak-base>0
ttestToneP = zeros(nTones,nNeuron);
ttestToneH = zeros(nTones,nNeuron);
ttestAlpha = 0.05;
% sign rank test for peak-base>0
signrankToneP = zeros(nTones,nNeuron);
signrankToneH = zeros(nTones,nNeuron);
signrankAlpha = 0.05;
% anova test for any group difference + at least one tone above baseline
anovaPeak = zeros(2,nNeuron);
anovaSignifTone = zeros(nTones,nNeuron);
anovaPeakCorr = zeros(2,nNeuron);
anovaSignifToneCorr = zeros(nTones,nNeuron);
% roc analysis for tone response
nShuffle = 100;
rocAuc = zeros(nTones,nNeuron);
rocTpr = zeros(nTones,nNeuron);
rocAucZscore = zeros(nTones,nNeuron);
%---------RUN SIGNIFICANCE TESTS-----------
peakAct = zeros(nTones,nTrials,nNeuron);
for i = 1:nNeuron
    tic;
    tempPeakAct = zeros(nTones,nTrials);
    for j = 1:nTones
        if testByPeak
            tempPeakAct(j,:) = squeeze(TCpretone_reorder(pretoneFrames+peakFrames(i,j),j,:,i));
        else
            tempPeakAct(j,:) = squeeze(nanmean(TCpretone_reorder(pretoneFrames+(1:floor(tonePeakLim*frameRate)),j,:,i),1));
        end
    end
    peakAct(:,:,i) = tempPeakAct;
    % Take first 5 frames before tone onset. (Since smoothwindow is generally 5)
    if testByPeak
        baseAct = squeeze(TCpretone_reorder(pretoneFrames,:,:,i)); 
    else
        baseAct = squeeze(nanmean(TCpretone_reorder(pretoneFrames+(-ceil(tonePeakLim*frameRate):-1),:,:,i),1)); 
    end
    % Do paired tests
    for j = 1:nTones
        try
            [h,p] = ttest(tempPeakAct(j,pairedTestTrials),baseAct(j,pairedTestTrials),'alpha',ttestAlpha,'tail','right');
            ttestToneP(j,i) = p; 
            ttestToneH(j,i) = h; 
        catch 
            disp(['Error in Cell ' int2str(i) 'Tone ' int2str(j) ' t test']);
            ttestToneP(j,i) = 1; % set p at 1
            ttestToneH(j,i) = 0; % set hypothesis wrong
        end
        
        try 
            [p,h,stats] = signrank(tempPeakAct(j,pairedTestTrials),baseAct(j,pairedTestTrials),'alpha',signrankAlpha,'tail','right'); 
            signrankToneP(j,i) = p;
            signrankToneH(j,i) = h;
        catch 
            disp(['Error in Cell ' int2str(i) 'Tone ' int2str(j) ' signrank test']);
            signrankToneP(j,i) = 1; % set p at 1
            signrankToneH(j,i) = 0; % set hypothesis wrong
        end
        
    end
    % Do ANOVA test
    
    groupNames = cell(1,nTones+1);
    groupNames{end} = 'Baseline';
    for j = 1:nTones
        groupNames{j} = int2str(toneorder(toneindex(j)));
    end
    [p,h,stats] = anova1(peakANOVA(:,:,i),groupNames,'off');
    anovaPeak (1,i) = p;
    [results,~,~,~] = multcompare(stats,'Display','off');
    results = results(results(:,2)==(nTones+1),6);
    anovaSignifTone(:,i) = (results<0.05);
    anovaPeak (2,i) = p<0.05 && sum(results<0.05)>0;

    [p,h,stats] = anova1(peakANOVACorr(:,:,i),groupNames,'off');
    anovaPeakCorr (1,i) = p;

    [results,~,~,~] = multcompare(stats,'Display','off');
    results = results(results(:,2)==(nTones+1),6);
    anovaSignifToneCorr(:,i) = (results<0.05);
    anovaPeakCorr (2,i) = p<0.05 && sum(results<0.05)>0;
    
    % Do ROC analysis
    % take only 5 frames before 
    baseAct = squeeze(TCpretone_reorder(pretoneFrames-4:pretoneFrames,:,:,i));
    for j = 1:nTones
        rocAct = [baseAct(:)' tempPeakAct(j,:)];
        rocLabel = [zeros(1,numel(baseAct)) ones(1,nTrials)];
        [tpr, fpr, threshold] = roc(rocLabel, rocAct);
        tempAuc = trapz([0 fpr 1],[0 tpr 1]);
        rocAuc(j,i) = tempAuc;
        tempFpr = find(fpr<0.05);
        rocTpr(j,i) = tpr(tempFpr(end));
        shuffAuc = zeros(1,nShuffle);
        for k = 1:nShuffle
            shuffLabel = rocLabel(randperm(length(rocLabel)));
            [tprShuff, fprShuff, thresholdShuff] = roc(shuffLabel, rocAct);
            shuffAuc(k) = trapz([0 fprShuff 1],[0 tprShuff 1]);
        end
        rocAucZscore(j,i) = (tempAuc - mean(shuffAuc)) / std(shuffAuc);
    end
    toc;
end

%---------SELECT CRITERIA FOR RESPONSIVE CELL-----------
responsiveCellFlag = anovaPeakCorr(2,:);

tuningPeakIndexMedian = zeros(1,nNeuron);
tuningPeakIndexMean = zeros(1,nNeuron);
popTuningPeakMedian = zeros(nTones,1);
popTuningPeakMean = zeros(nTones,1);
popTuningMedian = zeros(nTones,nNeuron);
popTuningMean = zeros(nTones,nNeuron);
for i = 1:nNeuron
    if responsiveCellFlag(i)

        [~,tempTuningPeakIndex] = max(squeeze(trialMedian(...
            pretoneFrames+peakFrames(i),:,i)));
        tuningPeakIndexMedian(i) = tempTuningPeakIndex;

        [~,tempTuningPeakIndex] = max(squeeze(trialMean(...
            pretoneFrames+peakFrames(i),:,i)));
        tuningPeakIndexMean(i) = tempTuningPeakIndex;

        popTuningMedian(:,i) = squeeze(trialMedian(pretoneFrames+peakFrames(i),:,i));
        popTuningMean(:,i) = squeeze(trialMean(pretoneFrames+peakFrames(i),:,i));
    else
        tuningPeakIndexMedian(i) = nan;
        tuningPeakIndexMean(i) = nan;
        popTuningMedian(:,i) = nan;
        popTuningMean(:,i) = nan;
    end
end

for i = 1:nTones
    popTuningPeakMedian(i) = sum(tuningPeakIndexMedian==i);
    popTuningPeakMean(i) = sum(tuningPeakIndexMean==i);
end

allDataName = {
    'peakFrames','peakActAvg','anovaPeak', 'anovaSignifTone',...
    'anovaPeakCorr','anovaSignifToneCorr',...
    'popTuningPeakMean', 'popTuningPeakMedian',...
    'rocAuc','rocTpr', 'rocAucZscore',...
    'signrankToneH', 'signrankToneP','signrankAlpha',...
    'ttestToneH', 'ttestToneP','ttestAlpha',...
    'tuningPeak', 'tuningPeakIndexMean', 'tuningPeakIndexMedian','peakAct'...
    'TC_original','TCpretone_reorder','TCpretone_reorderCorr',...
    'trialMean','trialMedian'};
save([savePath '/tuning.mat'],allDataName{:});

end