function getTuningFromSuite2p(mouse,celltype)
warning off
if strcmp(mouse,'se051')
    suite2ppath = 'X:\LabData5\se051\h5_files\day1\suite2p\';
    h5path = 'X:\LabData5\se051\h5_files\day1\';
end

if strcmp(mouse,'cd041')
    suite2ppath = 'U:\LabData4\celine\cd041\suite2p\';
    h5path = 'U:\LabData4\celine\cd041\';
end

if strcmp(mouse,'cd042')
    suite2ppath = 'U:\LabData4\celine\cd042\suite2p\';
    h5path = 'U:\LabData4\celine\cd042\';
end

if strcmp(mouse,'cd043')
    suite2ppath = 'U:\LabData4\celine\cd043\suite2p\';
    h5path = 'U:\LabData4\celine\cd043\';
end

if strcmp(mouse,'cd044')
    suite2ppath = 'U:\LabData4\celine\cd044\suite2p\';
    h5path = 'U:\LabData4\celine\cd044\';
end

if strcmp(mouse,'se058')
    suite2ppath = 'E:\se058\suite2p\';
    h5path = 'E:\se058\';
end

if strcmp(mouse,'zz008')
    suite2ppath = 'C:\Users\zzhu34\Documents\tempdata\zz008\suite2p\';
    h5path = 'C:\Users\zzhu34\Documents\tempdata\zz008\';
end

if strcmp(mouse,'zz012FOV0')
    suite2ppath = 'C:\Users\zzhu34\Documents\tempdata\zz012\zz012_000_000_FOV0\suite2p\';
    h5path = 'C:\Users\zzhu34\Documents\tempdata\zz012\zz012_000_000_FOV0\';
end


if strcmp(mouse,'zz012FOV1')
    suite2ppath = 'C:\Users\zzhu34\Documents\tempdata\zz012\zz012_000_002_FOV1\suite2p\';
    h5path = 'C:\Users\zzhu34\Documents\tempdata\zz012\zz012_000_002_FOV1';
end

if strcmp(mouse,'zz012FOV2')
    suite2ppath = 'C:\Users\zzhu34\Documents\tempdata\zz012\zz012_000_004_FOV2\suite2p\';
    h5path = 'C:\Users\zzhu34\Documents\tempdata\zz012\zz012_000_004_FOV2';
end

if strcmp(mouse,'zz012FOV3')
    suite2ppath = 'C:\Users\zzhu34\Documents\tempdata\zz012\zz012_000_006_FOV3\suite2p\';
    h5path = 'C:\Users\zzhu34\Documents\tempdata\zz012\zz012_000_006_FOV3';
end

if strcmp(mouse,'zz012FOVFINAL14')
    suite2ppath = 'D:\suite2p_temp\zz012\zz012_000_009_FOVFINAL14\suite2p\';
    h5path = 'D:\suite2p_temp\zz012\zz012_000_009_FOVFINAL14\';
end

if strcmp(mouse,'zz012FOVFINAL17')
    suite2ppath = 'C:\Users\zzhu34\Documents\tempdata\zz012\zz012_000_010_FOVFINAL17\suite2p\';
    h5path = 'C:\Users\zzhu34\Documents\tempdata\zz012\zz012_000_010_FOVFINAL17\';
end

cd(h5path);
files = dir('*.h5');
name = files(1).name(1:end-3);
data = load([name '.mat']); 
infosbx = data.info;

cd(suite2ppath)

tuningFolderName = ['TuningCurve_' name];
if exist([suite2ppath '/' tuningFolderName],'dir') ~= 7
    mkdir(tuningFolderName);
    mkdir([tuningFolderName '/singleNeuron']);
    mkdir([tuningFolderName '/population']);
end

% Default values
nTones = 17;
if isempty(infosbx.otparam)
    nPlanes = 1;
else
    nPlanes = infosbx.otparam(3);
end
nFramesPerTone = 100/nPlanes; % 25
nFramesPerTrial = nFramesPerTone * nTones; % 425
nTrials = 10;
startTrial = 2; % get rid of first 2 trials where overall fluor is high
nFrames = nFramesPerTrial*nTrials; % 4250

% Get TC
TC = []; neuronEachPlane = nan(nPlanes,1);
roisCoord = [];
for i=1:nPlanes
    cd([suite2ppath 'plane' num2str(i-1)]);
    data = load('Fall.mat'); 
    tc = data.F;
    %tc = data.spks;
    
    if strcmp(celltype,'axon')
        clear iscelltmp
        iscelltmp = zeros(length(data.iscell(:,1)),1);
        iscelltmp((data.iscell(:,1)==0)) = 1;
        data.iscell(:,1) = iscelltmp;
    end
    
    ok = data.iscell(:,1);
    
    if size(TC,1)>0,
        TC = [TC [tc(logical(ok),:)';nan(1,sum(ok))]];
    else
        TC = [TC tc(logical(ok),:)']; % time by neuron
    end
    neuronEachPlane(i) = sum(ok);
    roisCoord = [roisCoord;data.stat(logical(ok))'];
end

nNeuron = size(TC,2);
    
completeTC = nan(nFrames+100,nNeuron); % just so that if TC size is a lot larger, it can still work
%completeTC(2:size(TC,1)+1,:) = TC; % shift one frame so that all the tone onset is on 101, 201, 301...
completeTC(1:size(TC,1),:) = TC; % shift one frame so that all the tone onset is on 101, 201, 301...
TC = completeTC(1:nFrames,:);
    
allNeuronMean = reshape(TC,nFramesPerTone,nTones,nTrials,nNeuron);
allNeuronMean = allNeuronMean(:,:,startTrial:end,:);
allNeuronMean = squeeze(nanmean(nanmean(nanmean(allNeuronMean,2),3),4));
    
[~,peakIndex] = max(smoothdata(allNeuronMean,'movmean',4));

startFrame = 4;%peakIndex-ceil(6/nPlanes);
endFrame = 8;%peakIndex+ceil(6/nPlanes);
peakFrames  = startFrame : endFrame;

% UIFlag = true;
% while UIFlag
% 
%     f = figure;
%     subplot('Position',[.1 .3 .8 .6])
%     plot(allNeuronMean);title('mean activity across neuron and tones'); hold on;
%     
%     plot(peakFrames,allNeuronMean(peakFrames),'LineWidth',2)
%     title(sprintf(['mean activity across neuron and tones.\nFrames Chosen: '...
%         int2str(peakFrames(1)) '-' int2str(peakFrames(end))]))
% 
%     p = uipanel(f,'Title','Control Panel',...
%                     'Position',[0.1 0.05 .8 .2],'Fontsize',12);
%     h1 = uicontrol(p,'Units','normalized','Position',...
%         [.2 .5 .6 .4],'String', sprintf('Confirm current frame selection'),'Fontsize',12,'Callback',@confirmCallBack);
%     h2 = uicontrol(p,'Units','normalized','Position',...
%         [.2 .1 .6 .4],'String',sprintf('Redo frame selection'),'Fontsize',12,'Callback',@reselectCallBack);
%     uiresume(gcbf);
%     UIFlag = false;
% 
%     uiwait(f);
%     close(f);
% end
 
%pretoneFrames = length(peakFrames);
pretoneFrames = 10;
baselineFrames = 5;

completeTCpretone = nan(nFrames+100,nNeuron); % just so that if TC size is a lot larger, it can still work
completeTCpretone(2+pretoneFrames:size(TC,1)+1+pretoneFrames,:) = TC; % shift one frame so that all the tone onset is on 101, 201, 301...
TCpretoneTemp = completeTCpretone(1:nFrames,:);

%TC(1:2,:)=NaN;
% tone onset frames are originally 0 100 200 300 400 etc., but shifted one to 101, 201, etc. in the index
toneorder = [45254.834 8000 13454.34264 4756.82846 5656.854249,...
    22627.417 64000 53817.37058 4000 9513.65692,...
    16000 6727.171322 19027.31384 26908.6852 32000,...
    11313.7085 38054.62768];
%toneorder = circshift(toneorder,1);
toneindex = [9;4;5;12;2;10;16;3;11;13;6;14;15;17;1;8;7];


% function confirmCallBack(hObject, eventdata, handles)
%     uiresume(gcbf);
%     UIFlag = false;
% end
% 
% function reselectCallBack(hObject, eventdata, handles)
%     
%     prompt = {'start frame','end frame'};
%     UItitle = 'Input';
%     dims = [1 35;1 35];
%     definput = {int2str(startFrame),int2str(endFrame)};
%     answer = inputdlg(prompt,UItitle,dims,definput);
%     startFrame = str2double(answer{1});
%     endFrame = str2double(answer{2});
%     uiresume(gcbf);
%     
% end

% cd(suite2ppath)
% roiFolderName = 'ROI';
% refFolderName = 'ROI/RefImg';
% neuronEachPlane = zeros(1,nPlanes); 
% for i = 1:nPlanes 
%     roiName = [filePath '/' roiFolderName '/' fileName '_roi' int2str(i) '.zip'];% CD 7/10/2019
%      
%     rois = ReadImageJROI(roiName);
%     neuronEachPlane(i) = length(rois);
% end
%%
TCreshape=reshape(TC,nFramesPerTone,nTones,nTrials,nNeuron); %reshape(tc,25,17,9) 25 = frames per tone, 17 = number of tones, 9 = trials
TCpretoneReshape = reshape(TCpretoneTemp,nFramesPerTone,nTones,nTrials,nNeuron);
TCreorder=zeros(size(TCreshape));
TCpretoneReorder = zeros(size(TCpretoneReshape));
for x=1:nTones
    index=toneindex(x);%+1; %the "+1" is because the first 20 frames are 38 kHz that wraps around
    TCreorder(:,x,:,:)=TCreshape(:,index,:,:);
    TCpretoneReorder(:,x,:,:)=TCpretoneReshape(:,index,:,:);
end
TCtone = TCreorder;
TCpretone = TCpretoneReorder;
TCreorder=permute(reshape(TCreorder,nFramesPerTrial,nTrials,nNeuron),[2 1 3]); % now TC reorder is the size of trials*frame a trial * neuron

% use the same basline for all trials, different 
% have sem in the tuning curve data
TCreorder = TCreorder(startTrial:end,:,:); % [nTrials-1,(nFramesPerTone*nTones),nNeuron]
baseline = prctile(reshape(TCreorder, [(nTrials-startTrial+1)*nFramesPerTone*nTones,nNeuron]),25);
TCreorder = TCreorder ./ repmat(reshape(baseline,[1 1 nNeuron]),nTrials-startTrial+1,nFramesPerTone*nTones,1);

TCpretone = TCpretone(:,:,startTrial:end,:);
TCpretone = TCpretone ./ repmat(reshape(baseline,[1 1 1 nNeuron]),nFramesPerTone,nTones,nTrials-startTrial+1,1);

TCtrialMean = squeeze(nanmean(TCreorder));
TCtrialMedian = squeeze(nanmedian(TCreorder));
TCtrialSEM = squeeze(nanstd(TCreorder)./sqrt(nTrials));

TCtoneMean = reshape(TCtrialMean,[nFramesPerTone,nTones,nNeuron]);
TCtoneMedian = reshape(TCtrialMedian,[nFramesPerTone,nTones,nNeuron]);

tuningData = reshape(TCreorder, [nTrials-startTrial+1, nFramesPerTone,nTones,nNeuron]);
tuningData = reshape(tuningData(:,peakFrames,:,:),[(nTrials-startTrial+1)*length(peakFrames),nTones,nNeuron]);
tuningMean = squeeze(nanmean(tuningData));
tuningMedian = squeeze(nanmedian(tuningData));
tuningSEM = squeeze(nanstd(tuningData)/sqrt(size(tuningData,1)));

% first do a Anova to test if the cell is responsive to any tones
% need 18 groups, each group with 8 data points for 8 trials
% the baseline group is averaged from 8 trials 

peakActANOVA = nan((nTrials-startTrial+1)*nTones,nTones+1,nNeuron);
%peakActANOVA(:,1,:) = squeeze(nanmean(nanmean(TCpretone(1:pretoneFrames,:,startTrial:nTrials,:))));
pretoneMean = (nanmean(TCpretone(1:pretoneFrames,:,:,:)));
TCpretoneCorr = TCpretone - repmat(pretoneMean, [nFramesPerTone, 1, 1, 1]);
for i = 1:nTones
    
    
    peakActANOVA(((i-1)*(nTrials-startTrial+1)+1):(i*(nTrials-startTrial+1)),i,:) = squeeze(nanmean(TCpretoneCorr(pretoneFrames+peakFrames,i,:,:))) ;
    peakActANOVA(((i-1)*(nTrials-startTrial+1)+1):(i*(nTrials-startTrial+1)),end,:) = squeeze(nanmean(TCpretoneCorr((pretoneFrames-baselineFrames+1):pretoneFrames,i,:,:)));
end

% then do a t test to test if the cell is any-responsive, peak tone responsive vs. baseline
[~,peakIndex] = max(tuningMedian);
ttestPeakTone = zeros(2,nNeuron);
ttestAllTone = zeros(2,nNeuron);
anovaAllTone = zeros(1,nNeuron);
Table = zeros(nTones+1, nNeuron +1);
neuronRespTable(:,1) = [sort(toneorder) 0];

signRankP = zeros(nNeuron,nTones);
% note fix the t-test here 
for i = 1:nNeuron
    peakAct = TCpretone(peakFrames+pretoneFrames,:,:,i);% i.e. 3,(8*17)
    peakActAvg = squeeze( nanmean(peakAct,1)); 
    baseAct = TCpretone(1:pretoneFrames,:,:,i);
    baseActAvg = squeeze( nanmean(baseAct,1));
    
    %[h,p] = ttest(peakAct(:),baseAct(:));
    %ttestAllTone(:,i) = [h;p]; 
        
    %peakActPrefTone = squeeze(peakAct(peakIndex(i),:));
    %baseActPrefTone = squeeze(baseAct(peakIndex(i),:));
    %[h,p] = ttest(peakActPrefTone,baseActPrefTone);
    %ttestPeakTone(:,i) = [h;p];

    
    for j = 1:nTones
        [p,h,stats] = signrank(peakActAvg(j,:),baseActAvg(j,:));
        %signRankH (i,j) = h;
        signRankP(i,j) = p;
    end
    
    groupNames = cell(1,nTones+1);
    groupNames{end} = 'Baseline';
    for j = 1:nTones
        groupNames{j} = int2str(toneorder(toneindex(j)));
    end
    [p,~,stats] = anova1(peakActANOVA(:,:,i),groupNames,'off');
    anovaAllTone (i) = p;
    
    for j = 1:nTones
        rocAct = [baseAct(:)' peakActAvg(j,:)];
        rocLabel = [zeros(1,numel(baseAct)) ones(1,nTrials)];
        [tp, fpr, threshold] = roc(rocAct, rocLabel);
        auc = trapz([0 fpr 1],[0 tpr 1]);
        auc_all(i,j) = auc;
        %nTruePos = ;
    end

    [results,~,~,~] = multcompare(stats);
    results = results(results(:,2)==18,6);
    signifTone = (results<0.05);
    
    neuronRespTable(1:nTones,i+1) = signifTone;
    neuronRespTable(end,i+1) = p < 0.05;

end

%responsiveCellFlag = ttestAllTone(2,:) <= 0.01;% | (ttestPeakTone(1,:) == 1) | (anovaAllTone < 0.05);
responsiveCellFlag = anovaAllTone < 0.05;

save([suite2ppath '/' tuningFolderName '/population/anovaAllTone.mat'],'anovaAllTone');
save([suite2ppath '/' tuningFolderName '/population/responsiveCellFlag.mat'],'responsiveCellFlag');
save([suite2ppath '/' tuningFolderName '/population/ttestPeakTone.mat'],'ttestPeakTone');
save([suite2ppath '/' tuningFolderName '/population/ttestAllTone.mat'],'ttestAllTone');
save([suite2ppath '/' tuningFolderName '/population/responsiveCellFlag.mat'],'responsiveCellFlag');
save([suite2ppath '/' tuningFolderName '/population/neuronRespTable.mat'],'neuronRespTable');
save([suite2ppath '/' tuningFolderName '/population/tuningMean.mat'],'tuningMean'); % JL 05/08/19 for tuning and bandwidth analysis
save([suite2ppath '/' tuningFolderName '/population/tuningMedian.mat'],'tuningMedian'); % JL 05/08/19 for tuning and bandwidth analysis
save([suite2ppath '/' tuningFolderName '/population/tuningData.mat'],'tuningData'); % JL 05/08/19 for tuning and bandwidth analysis
%% population analysis

popTuningPeakMedian = zeros(nTones,1);
[~,peakIndex] = max(tuningMedian(:,responsiveCellFlag));
for i = 1:nTones
    popTuningPeakMedian(i) = sum(peakIndex==i);
end

popTuningPeakMean= zeros(nTones,1);
[~,peakIndex] = max(tuningMean(:,responsiveCellFlag));
for i = 1:nTones
    popTuningPeakMean(i) = sum(peakIndex==i);
end

popTuningMedian = nanmean(tuningMedian(:,responsiveCellFlag)');
popTuningMean = nanmean(tuningMean(:,responsiveCellFlag)');

%popTuningStd = std(tuningMean');

figure;

subplot(2,2,2)
toneRespCount = sum(neuronRespTable(1:end-1,2:end),2);
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
% saveas(gcf,[filePath '/' tuningFolderName...
%         '/population/peakCount.png']);

subplot(2,2,3)
%errorbar(freqAxis, popTuning, popTuningStd,'LineWidth',2);
plot(freqAxis, popTuningMedian,'LineWidth',2); hold on;
plot(freqAxis, popTuningMean,'LineWidth',2);
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
saveas(gcf,[suite2ppath '/' tuningFolderName ...
        '/population/populationTuning.png']);
saveas(gcf,[suite2ppath '/' tuningFolderName ...
        '/population/populationTuning.m']);
  
%%

[~,peakIndex] = max(tuningMedian);
neuronPlane = cumsum(neuronEachPlane);
neuronPlane = [0 neuronPlane'];
figure;
for i = 1:nPlanes
    cd([suite2ppath 'plane' num2str(i-1)]);
    data = load('Fall.mat'); 
    refImg = data.ops.meanImg;
    
    
    colormapIndex = round(linspace(1,64,17));
    C = colormap('jet');
    
    if i==1
        rois = roisCoord(1:neuronEachPlane(i));
    else
        rois = roisCoord(neuronEachPlane(i-1)+1:neuronPlane(end)); 
    end
    
    subplot(2,2,i);           
    imagesc(refImg);colormap gray;hold on;
    ylim([0 size(refImg,1)]);xlim([0 size(refImg,2)]);
    
    subplot(2,2,2+i);             
    imagesc(refImg);colormap gray;hold on;
    ylim([0 size(refImg,1)]);xlim([0 size(refImg,2)]);

    for j=(neuronPlane(i)+1):neuronPlane(i+1)
        x=rois{j-neuronPlane(i)}.xpix; %freehand rois have the outlines in x-y coordinates
        y=rois{j-neuronPlane(i)}.ypix; %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
        if responsiveCellFlag(j)
             plot(x,y,'.','color',C(colormapIndex(peakIndex(j-neuronPlane(i))),:),'MarkerSize',1);
        else 
            plot(x,y,'.','color',[0.8 0.8 0.8],'MarkerSize',1);
        end
    end
    title(['Plane' int2str(i)])
end

saveas(gcf,[suite2ppath '/' tuningFolderName...
        '/population/tuningMap.png']);
saveas(gcf,[suite2ppath '/' tuningFolderName...
    '/population/tuningMap.m']);

toneRespTable = zeros(nTones, 4);
toneRespTable(:,1) = sort(toneorder);
toneRespTable(:,2) = sum(neuronRespTable(1:end-1,2,:));
toneRespTable(:,3) = popTuningMedian;
toneRespTable(:,4) = popTuningMean;
save([suite2ppath '/' tuningFolderName '/population/toneRespTable.mat'],'toneRespTable');



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
        for i = 1:neuronEachPlane
            tic;

            figure('visible','off');
            cellIndex = currentNeuron(j) + i;
            % REMEMBER TO CHANGE THIS LINE
            cd([suite2ppath 'plane' num2str(j-1)]);
            data = load('Fall.mat'); 
            refImg = data.ops.meanImg;
            subplot(2,2,1)
            imagesc(refImg);colormap gray;hold on;

            x=roisCoord{i+currentNeuron(j)}.xpix; %freehand rois have the outlines in x-y coordinates
            y=roisCoord{i+currentNeuron(j)}.ypix; %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
            plot(x,y,'.','color',[0.8 0.8 0.8],'MarkerSize',3);

            %x=rois{1,i}.mnCoordinates(:,1); %freehand rois have the outlines in x-y coordinates
            %y=rois{1,i}.mnCoordinates(:,2); %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
            %patch(x,y,'g','EdgeColor','none');
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
                plot(timeAxis,TCreorder(k,1:tempSampleRate:nFramesPerTrial,cellIndex),'color',[0.9 0.9 0.9],'LineWidth',0.5); hold on
            end

            timeAxis = (0:nFramesPerTrial-1) * nPlanes / 30;
            p1 = plot(timeAxis, TCtrialMedian(:,cellIndex),'LineWidth',1.2,'color',[0.0000 0.4470 0.7410]);
            p2 = plot(timeAxis, TCtrialMean(:,cellIndex),'LineWidth',1.2,'color',[0.8500 0.3250 0.0980]);
            yMax = 0.7 * max(TCtrialMean(:,cellIndex)) + 0.3 * max(max(TCreorder(:,1:tempSampleRate:nFramesPerTrial,cellIndex)));
            yMin = 0.7 * min(TCtrialMean(:,cellIndex)) + 0.3 * min(min(TCreorder(:,1:tempSampleRate:nFramesPerTrial,cellIndex)));
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

            saveas(gcf,[suite2ppath '/' tuningFolderName ...
                '/singleNeuron/Neuron' num2str(cellIndex,'%03d') '.png']);
            saveas(gcf,[suite2ppath '/' tuningFolderName ...
                '/singleNeuron/Neuron' num2str(cellIndex,'%03d') '.m']);
            close all;
            timeElapsed = toc;

            disp(['Cell #' int2str(cellIndex) ' Plane #' int2str(j) ' Time=' num2str(timeElapsed,'%03f') ' secs'])

        end
    end
end

warning on

end