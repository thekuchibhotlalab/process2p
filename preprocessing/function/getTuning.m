function getTuning(varargin)
%myFun - Description
%
% Syntax: Func_getTuning(input)
%
% Long description

p = func_createInputParser();
p.parse(varargin{:});
sep = '\';
%---------CHECK NUMBER OF FRAMES IN SBX FILE-----------
global nTones nTrials nFramesPerTone nFramesPerTrial startTrial nFrames nPlanes

nPlanes = str2double(p.Results.nPlanes);
filename = p.Results.filename;

if iscell(filename) && length(filename) >1
    disp('ERROR - More than one session selected!')
end

%cd(p.Results.sbxpath);
data = load([p.Results.sbxpath sep filename '.mat']); 
infosbx = data.info;
if isempty(infosbx.otparam)
    check_nPlanes = 1;
else
    check_nPlanes = infosbx.otparam(3);
end
if check_nPlanes ~= nPlanes
    disp('ERROR - nPlanes in the parameter not consistent with sbx file.')
    pause;
end

%---------CHECK NUMBER OF CHANNELS-----------
if iscell(p.Results.functionalChannel)
    nFuncChannel = length(p.Results.functionalChannel);
    functionalChannel = p.Results.functionalChannel;
    roiType = p.Results.roiType{1};
else
    nFuncChannel = 1;
    functionalChannel = {p.Results.functionalChannel};
    roiType = p.Results.roiType;
end

%---------GET CALCIUM TRACES-----------
switch p.Results.roiMethod
    case 'suite2p'
        [suite2pTC, neuronEachPlane, roisCoord] = func_loadTCsuite2p(varargin{:},'file',1);
        roisCoord = {roisCoord};
    case 'manual'
        roiFileSplit = strsplit(p.Results.roiFile);
        roiFile = reshape(roiFileSplit,nFuncChannel,nPlanes);
        roisCoord = cell(1,nFuncChannel);
        for i = 1:nFuncChannel
            for j = 1:nPlanes
                roiName = [p.Results.datapath sep roiFile{i,j}];
                roisCoord{i} = [roisCoord{i} ReadImageJROI(roiName)];
            end
        end
end
switch p.Results.tcMethod
    case 'suite2p'
        TC = {suite2pTC};
        neuronEachPlane = {neuronEachPlane};
    case 'manual'
        [TC, neuronEachPlane] = func_loadTCmanual(varargin{:},'file',1);
        
end
%---------CHECK IF NUMBER OF CHANNELs IS CORRECT-----------
if length(TC) ~= nFuncChannel
    disp('ERROR - Number of functional channels not correct')
end
%---------GET TUNING DATA FOR EACH CHANNEL-----------
if length(TC) == 1
    %CREATE SAVE FOLDER
    savePath = [p.Results.savepath sep p.Results.filename '_Tuning'];
    if exist(savePath,'dir') ~= 7
        mkdir(savePath);
        mkdir([savePath '/singleNeuron']);
        mkdir([savePath '/population']);
    end
    %PROCESS TUNING DATA
    getTuning_oneChannel(TC{1},neuronEachPlane{1},roisCoord{1},savePath,p.Results.suite2ppath);

elseif length(TC) == 2
    for i = 1:nFuncChannel        
        %CREATE SAVE FOLDER
        savePath = [p.Results.savepath sep p.Results.filename '_Tuning' sep 'chan' int2str(i)];
        if exist(savePath,'dir') ~= 7
            mkdir(savePath);
            mkdir([savePath '/singleNeuron']);
            mkdir([savePath '/population']);
        end
        %PROCESS TUNING DATA
        getTuning_oneChannel(TC{i},neuronEachPlane{i},roisCoord{i},savePath,p.Results.suite2ppath);
    end
end
%---------END OF GET TUNING FUNCTION-----------
end

%---------FUNCTION TO PROCESS CALCIUM TRACES-----------
function getTuning_oneChannel(TC,neuronEachPlane,roisCoord,savePath,suite2ppath)

%---------DECLARE SOME PARAMETERS-----------
global nPlanes
nTones = 17;
nTrials = 10;
nFramesPerTone = 100/nPlanes; % 50
nFramesPerTrial = nFramesPerTone * nTones; % 850
startTrial = 2; % the first tone is on frame 0
nFrames = nFramesPerTrial*nTrials; % 4250
frameRate = 30.98/nPlanes; % in the future, do not hard code this.
pretoneFrames = 10;
baselineFrames = 5;
nNeuron = size(TC,2);
smoothWindow = 3;

%---------SHIFT THE TC FOR PRETONE PERIOD-----------
TC_original = TC; % keep a copy of original TC
TC = circshift(TC,1,1); % shift 1 on 1st axix so that tone is on frame 1 instead of 0
%---------GAUSSIAN FILTER TO SMOOTH THE TRACES-----------
TC = smoothdata(TC,1,'gaussian',smoothWindow);
TCpretone = circshift(TC,pretoneFrames,1);
%---------COMPUTE DFF-----------
baseline = prctile(TC,50,1);
TC = TC ./ repmat(baseline,[size(TC,1) 1]);
TCpretone = TCpretone ./ repmat(baseline,[size(TC,1) 1]);
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
%---------PLOT PEAK FRAME-----------
[maxValue,peakIndex] = max(toneMean(pretoneFrames+1:pretoneFrames+1+frameRate,:),[],1);
latFig = figure;
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
%
individualPeakFrameFlag = true;
if individualPeakFrameFlag
    peakFrames = peakIndex;
else
    peakFrames = ones(size(peakIndex)) * mode(peakIndex);
end

%startFrame = 4;%peakIndex-ceil(6/nPlanes);
%endFrame = 8;%peakIndex+ceil(6/nPlanes);
%peakFrames  = startFrame : endFrame;



%TCtone = TCreorder;
%TCpretone = TCpretoneReorder;
%TCreorder=permute(reshape(TCreorder,nFramesPerTrial,nTrials,nNeuron),[2 1 3]); % now TC reorder is the size of trials*frame a trial * neuron

% use the same basline for all trials, different 
% have sem in the tuning curve data
%Creorder = TCreorder(startTrial:end,:,:); % [nTrials-1,(nFramesPerTone*nTones),nNeuron]
%baseline = prctile(reshape(TCreorder, [(nTrials-startTrial+1)*nFramesPerTone*nTones,nNeuron]),25);
%TCreorder = TCreorder ./ repmat(reshape(baseline,[1 1 nNeuron]),nTrials-startTrial+1,nFramesPerTone*nTones,1);

%TCpretone = TCpretone(:,:,startTrial:end,:);
%TCpretone = TCpretone ./ repmat(reshape(baseline,[1 1 1 nNeuron]),nFramesPerTone,nTones,nTrials-startTrial+1,1);

%TCtrialMean = squeeze(nanmean(TC_reorder));
%TCtrialMedian = squeeze(nanmedian(TC_reorder));
%TCtrialSEM = squeeze(nanstd(TCreorder)./sqrt(nTrials));

%TCtoneMean = reshape(TCtrialMean,[nFramesPerTone,nTones,nNeuron]);
%TCtoneMedian = reshape(TCtrialMedian,[nFramesPerTone,nTones,nNeuron]);

%tuningData = reshape(TCreorder, [nTrials-startTrial+1, nFramesPerTone,nTones,nNeuron]);
%tuningData = reshape(tuningData(:,peakFrames,:,:),[(nTrials-startTrial+1)*length(peakFrames),nTones,nNeuron]);
%tuningMean = squeeze(nanmean(tuningData));
%tuningMedian = squeeze(nanmedian(tuningData));
%tuningSEM = squeeze(nanstd(tuningData)/sqrt(size(tuningData,1)));

% first do a Anova to test if the cell is responsive to any tones
% need 18 groups, each group with 8 data points for 8 trials
% the baseline group is averaged from 8 trials 

peakANOVA = nan((nTrials)*nTones,nTones+1,nNeuron);
peakANOVACorr = nan((nTrials)*nTones,nTones+1,nNeuron);
%peakANOVA(:,1,:) = squeeze(nanmean(nanmean(TCpretone(1:pretoneFrames,:,startTrial:nTrials,:))));
pretoneMean = (nanmean(TCpretone_reorder(1:pretoneFrames,:,:,:)));
TCpretone_reorderCorr = TCpretone_reorder - repmat(pretoneMean, [nFramesPerTone, 1, 1, 1]);
for i = 1:nTones
    for j = 1:nNeuron
        % Although most times we are taking 1 frame, but use nanmean here to allow multiple frames
        peakANOVA(((i-1)*(nTrials)+1):(i*(nTrials)),i,j) = ...
            squeeze(nanmean(TCpretone_reorder(pretoneFrames+peakFrames(j),i,:,j),1));
        peakANOVA(((i-1)*(nTrials)+1):(i*(nTrials)),end,j) = ...
            squeeze(nanmean(TCpretone_reorder(pretoneFrames,i,:,j),1));

        peakANOVACorr(((i-1)*(nTrials)+1):(i*(nTrials)),i,j) = ...
            squeeze(nanmean(TCpretone_reorderCorr(pretoneFrames+peakFrames(j),i,:,j),1));
        peakANOVACorr(((i-1)*(nTrials)+1):(i*(nTrials)),end,j) = ...
            squeeze(nanmean(TCpretone_reorderCorr(pretoneFrames,i,:,j),1));
    end
end

for i = 1:nNeuron
    [~,tempTuningPeak] = max(squeeze(trialMedian(pretoneFrames+peakFrames(i),:,i)));
    tuningPeak(i) = tempTuningPeak;
end
%---------DECLARE ALL MATRICES FOR SIGNIFICANT TESTS-----------
pairedTestTrials = 2:10; % trial 1 does not have pair for baseline
% paired ttest for peak-base>0
ttestToneP = zeros(nTones,nNeuron);
ttestToneH = zeros(nTones,nNeuron);
ttestAlpha = 0.005;
% sign rank test for peak-base>0
signrankToneP = zeros(nTones,nNeuron);
signrankToneH = zeros(nTones,nNeuron);
signrankAlpha = 0.005;
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
for i = 1:nNeuron
    tic;
    peakAct = squeeze(TCpretone_reorder(pretoneFrames+peakFrames(i),:,:,i));
    baseAct = squeeze(TCpretone_reorder(pretoneFrames,:,:,i));
    % Do paired tests
    for j = 1:nTones
        [h,p] = ttest(peakAct(j,pairedTestTrials),baseAct(j,pairedTestTrials),'alpha',ttestAlpha,'tail','right');
        ttestToneP(j,i) = p; 
        ttestToneH(j,i) = h; 
        [p,h,stats] = signrank(peakAct(j,pairedTestTrials),baseAct(j,pairedTestTrials),'alpha',signrankAlpha,'tail','right');
        signrankToneP(j,i) = p;
        signrankToneH(j,i) = h;
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
    for j = 1:nTones
        baseAct = squeeze(TCpretone_reorder(pretoneFrames-4:pretoneFrames,:,:,i));
        rocAct = [baseAct(:)' peakAct(j,:)];
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
allDataName = {'ttestToneP','ttestToneH','ttestAlpha','signrankToneP',...
'signrankToneH','signrankAlpha','anovaPeak','anovaSignifTone',...
'anovaPeakCorr','anovaSignifToneCorr','rocAuc','rocTpr','rocAucZscore'};
save([savePath '/population/tuning.mat'],allDataName{:});
%---------PLOT SIGNICANT TEST RESULTS ON POPULATION LEVEL-----------
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

%popTuningPeakMean= zeros(nTones,1);
%[~,peakIndex] = max(squeeze(trialMean(:,responsiveCellFlag)));
%for i = 1:nTones
%    popTuningPeakMean(i) = sum(peakIndex==i);
%end

%popTuningMedian = squeeze(nanmean(trialMedian(:,responsiveCellFlag),3));
%popTuningMean = squeeze(nanmean(trialMean(:,responsiveCellFlag),3));

%popTuningStd = std(tuningMean');

figure;

subplot(2,2,2)
toneRespCount = sum(anovaSignifToneCorr,2);
freqAxis = log2(sort(toneorder));
plot(freqAxis, toneRespCount,'LineWidth',2);
set(gca, 'XTick', log2([4000 8000 16000 32000 64000]));
set(gca, 'XTickLabel', [4 8 16 32 64]);
xlim([log2(4000) log2(64000)]);
xlabel('Frequency (kHz)');
ylabel('Cell Count')
title('Cell Count by Significance to Each Tone')


subplot(2,2,4)
freqAxis = log2(sort(toneorder));
plot(freqAxis, popTuningPeakMedian,'LineWidth',2); hold on;
plot(freqAxis, popTuningPeakMean,'LineWidth',2); 
set(gca, 'XTick', log2([4000 8000 16000 32000 64000]));
set(gca, 'XTickLabel', [4 8 16 32 64]);
xlim([log2(4000) log2(64000)]);
xlabel('Frequency (kHz)');
ylabel('Cell Count')
title('Cell Count by Peak Tone')
legend('Median','Mean')
% saveas(gcf,[filePath '/' savePath...
%         '/population/peakCount.png']);

subplot(2,2,3)
%errorbar(freqAxis, popTuning, popTuningStd,'LineWidth',2);
plot(freqAxis, nanmean(popTuningMedian,2),'LineWidth',2); hold on;
plot(freqAxis, nanmean(popTuningMean,2),'LineWidth',2);
set(gca, 'XTick', log2([4000 8000 16000 32000 64000]));
set(gca, 'XTickLabel', [4 8 16 32 64]);
xlim([log2(4000) log2(64000)]);
xlabel('Frequency (kHz)');
ylabel('mean F/F0');
legend('Median','Mean')
title('Mean Tuning Activity')

subplot(2,2,1)
percentResp = [sum(responsiveCellFlag), sum(~responsiveCellFlag)] ./ nNeuron;
pie([sum(responsiveCellFlag), sum(~responsiveCellFlag)],{['resp ' int2str(sum(responsiveCellFlag))],...
    ['not resp ' int2str(sum(~responsiveCellFlag))]});
title([int2str(percentResp(1)*100) '% cell Responsive'])
saveas(gcf,[ savePath ...
        '/population/populationTuning.png']);
%saveas(gcf,[ savePath ...
%        '/population/populationTuning.m']);
  
%%

%[~,peakIndex] = max(tuningMedian);
neuronPlane = cumsum(neuronEachPlane);
neuronPlane = [0 neuronPlane];
figure;
for i = 1:nPlanes
    cd([suite2ppath '\plane' num2str(i-1)]);
    data = load('Fall.mat'); 
    refImg = data.ops.meanImg;
    
    
    colormapIndex = round(linspace(1,64,17));
    C = colormap('jet');
    
    %if i==1
    %    rois = roisCoord(1:neuronEachPlane(i));
    %else
    %    rois = roisCoord(neuronEachPlane(i-1)+1:neuronPlane(end)); 
    %end
    
    subplot(2,2,i);           
    imagesc(refImg);colormap gray;hold on;
    ylim([0 size(refImg,1)]);xlim([0 size(refImg,2)]);
    
    subplot(2,2,2+i);             
    imagesc(refImg);colormap gray;hold on;
    ylim([0 size(refImg,1)]);xlim([0 size(refImg,2)]);

    for j=(neuronPlane(i)+1):neuronPlane(i+1)
        x = roisCoord{1,j}.mnCoordinates(:,1); %freehand rois have the outlines in x-y coordinates
        y = roisCoord{1,j}.mnCoordinates(:,2); %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
        if responsiveCellFlag(j)
             plot(x,y,'.','color',C(colormapIndex(peakIndex(j-neuronPlane(i))),:),'MarkerSize',1);
        else 
            plot(x,y,'.','color',[0.8 0.8 0.8],'MarkerSize',1);
        end
    end
    title(['Plane' int2str(i)])
end

saveas(gcf,[savePath...
        '/population/tuningMap.png']);
%saveas(gcf,[tuningFolderName...
%    '/population/tuningMap.m']);

toneRespTable = zeros(nTones, 4);
toneRespTable(:,1) = sort(toneorder);
toneRespTable(:,2) = sum(neuronRespTable(1:end-1,2,:));
toneRespTable(:,3) = popTuningMedian;
toneRespTable(:,4) = popTuningMean;
save([savePath '/population/toneRespTable.mat'],'toneRespTable');



%% single neuron analysis 
%plotNeuron = nNeuron;
%cellPerPlane = zeros(1,nPlanes);
if true
    currentNeuron = [0 neuronEachPlane(1)];
    for j = 1:nPlanes
    %     roiName = [ filePath '/' roiFolderName '/' fileName '_roi' int2str(j),suffix,'.zip'];
        %roiName = [ filePath '/' roiFolderName '/' fileName '_roi' int2str(j) '.zip'];
        %roiName = [filePath '/' roiFolderName '/' fileName '_ot_' num2str(whichPlane-1,'%03d') '_rois.zip'];
        %rois = ReadImageJROI(roiName);
        %cellThisPlane = length(rois);
        %startingIndex = sum(cellPerPlane);
        %cellPerPlane(j) = cellThisPlane;
        for i = 1:neuronEachPlane(j)
            tic;

            tuningFig = figure('visible','off');
            cellIndex = currentNeuron(j) + i;
            % REMEMBER TO CHANGE THIS LINE
            cd([suite2ppath '\plane' num2str(j-1)]);
            data = load('Fall.mat'); 
            refImg = data.ops.meanImg;
            subplot(2,2,1)
            imagesc(refImg);colormap gray;hold on;

            %x=roisCoord{i+currentNeuron(j)}.xpix; %freehand rois have the outlines in x-y coordinates
            %y=roisCoord{i+currentNeuron(j)}.ypix; %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
            %plot(x,y,'.','color',[0.8 0.8 0.8],'MarkerSize',3);

            x=roisCoord{1,cellIndex}.mnCoordinates(:,1); %freehand rois have the outlines in x-y coordinates
            y=roisCoord{1,cellIndex}.mnCoordinates(:,2); %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
            patch(x,y,'g','EdgeColor','none');
            title(['Cell #' int2str(cellIndex) ' Plane #' int2str(j) ' FoV'])

            subplot(2,2,2)
            plot(TC(:,cellIndex),'LineWidth',0.5);
            xlim([1 size(TC,1)]);
            xlabel('Frames')
            ylabel('Fluorescence')
            title('Raw Fluorescence')

            subplot(2,2,3)
            for k = 1:nTrials-startTrial+1
                tempSampleRate = 2;
                timeAxis = (0:tempSampleRate:nFramesPerTrial-1) * nPlanes / 30;
                plot(timeAxis,TC_reorder(k,1:tempSampleRate:nFramesPerTrial,cellIndex),'color',[0.9 0.9 0.9],'LineWidth',0.5); hold on
            end

            timeAxis = (0:nFramesPerTrial-1) * nPlanes / 30;
            p1 = plot(timeAxis, trialMedian(peakFrames(cellIndex),:,cellIndex),'LineWidth',1.2,'color',[0.0000 0.4470 0.7410]);
            p2 = plot(timeAxis, trialMean(peakFrames(cellIndex),:,cellIndex),'LineWidth',1.2,'color',[0.8500 0.3250 0.0980]);
            yMax = 0.7 * max(trialMean(peakFrames(cellIndex),:,cellIndex)) + 0.3 * max(max(TCreorder(:,1:tempSampleRate:nFramesPerTrial,cellIndex)));
            yMin = 0.7 * min(trialMean(peakFrames(cellIndex),:,cellIndex)) + 0.3 * min(min(TCreorder(:,1:tempSampleRate:nFramesPerTrial,cellIndex)));
            onsetTime = (1:nFramesPerTone:nFramesPerTrial) * nPlanes / 30;
            for k = 1:nTones
                plot([onsetTime(k) onsetTime(k)],[yMin yMax],'LineWidth',0.5,'color',[0.4 0.4 0.4]);
            end
            legend([p1, p2],'Median','Mean');
            xlim([0 timeAxis(end)])
            ylim([yMin yMax])
            ylabel('F/F0')
            xlabel('Time(s)')
            title('Mean Activity Across Trials')

            signifTonePoint = neuronRespTable(1:end-1,[1 cellIndex+1]);
            signifTuningMedian = tuningMedian(:,cellIndex);
            signifTuningMedian = signifTuningMedian(signifTonePoint(:,2)==1);
            signifTonePoint = signifTonePoint(signifTonePoint(:,2)==1,1);
            subplot(2,2,4)
            freqAxis = log2(sort(toneorder));
            magicNum = sqrt(pi/2);
            errorbar(freqAxis, tuningMedian(:,cellIndex),magicNum * tuningSEM(:,cellIndex),'LineWidth',2,'color',[0.0000 0.4470 0.7410]); hold on;
            errorbar(freqAxis, tuningMean(:,cellIndex),tuningSEM(:,cellIndex),'LineWidth',2,'color',[0.8500 0.3250 0.0980]);
            scatter(log2(signifTonePoint),signifTuningMedian,40,[0.2 0.2 0.2],'*');
            set(gca, 'XTick', log2([4000 8000 16000 32000 64000]));
            set(gca, 'XTickLabel', [4 8 16 32 64]);
            xlabel('Frequency (kHz)');
            ylabel('F/F0');
            xlim([log2(4000) log2(64000)]);
            legend('Median','Mean');
            title('Tuning Curve');

            if responsiveCellFlag(cellIndex)
                ylimit = get(gca,'YLim');
                scatter(freqAxis(2),ylimit(2)*0.9+ylimit(1)*0.1,40,[0 0 0],'*')
            end

            saveas(tuningFig,[savePath ...
                '/singleNeuron/Neuron' num2str(cellIndex,'%03d') '.png']);
                
            close(tuningFig);
            timeElapsed = toc;

            disp(['Cell #' int2str(cellIndex) ' Plane #' int2str(j) ' Time=' num2str(timeElapsed,'%03f') ' secs'])

        end
    end
end



end