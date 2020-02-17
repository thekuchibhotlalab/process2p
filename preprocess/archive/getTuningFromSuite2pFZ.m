function getTuningFromSuite2pFZ(mouse,site,celltype,file,channel)%,plane)

suite2ppath = ['Q:\Fangchen\' mouse '\' site '\suite2p\'];
h5path = ['Q:\Fangchen\' mouse '\' site '\'];

cd(h5path);
files = dir('*.h5');
name = files(1).name(1:end-3);
data = load([name '.mat']); 
infosbx = data.info;

cd(suite2ppath)
tuningFolderName = ['TuningCurve_' name '_' channel];
if exist([suite2ppath '/' tuningFolderName],'dir') ~= 7
    mkdir(tuningFolderName);
    mkdir([tuningFolderName '/singleNeuron']);
    mkdir([tuningFolderName '/population']);
end

% Set parameters
nTones = 17;
if isempty(infosbx.otparam)
    nPlanes = 1;
else
    nPlanes = infosbx.otparam(3);
end
%Note: we are presenting each 17 tones 10 times each to the animal
nFramesPerTone = 100/nPlanes; % 50 for 2 planes
nFramesPerTrial = nFramesPerTone * nTones; % 850 for 17 tones
nTrials = 10;
startTrial = 1; % don't get rid of first 2 trials where overall fluor is high
nFrames = nFramesPerTrial*nTrials; % 8500

% Get tuning curve (TC)
% Initialize matrices
TC = []; 
neuronEachPlane = zeros(nPlanes,1);
roisCoord = [];
for i=1:nPlanes 
    %selecting ROIs
    if strcmp(file,'suite2p')   % use suite2p ROIs
        cd([suite2ppath 'plane' num2str(i-1)]);
        data = load('Fall.mat');
        tc = data.F;
        roisCoord = [roisCoord;data.stat(logical(ok))'];
    elseif strcmp(file,'handroi')   % use handdrawn ROIs
        cd([suite2ppath 'plane' num2str(i-1)]);
        data = load([mouse '_TC_plane'  num2str(i-1) '_' channel '.mat']);
        tc = data.tempTC;
        roiName = [suite2ppath 'plane' num2str(i-1) '\' 'roi_' channel '.zip'];
        roisCoord = [roisCoord ReadImageJROI(roiName)];
    end
    
    %selecting cell type (only if using suite2p)
    if strcmp(celltype,'axon') %labelled 'not cell' in suite2p
        clear iscelltmp
        iscelltmp = zeros(length(data.iscell(:,1)),1);
        iscelltmp((data.iscell(:,1)==0)) = 1;
        data.iscell(:,1) = iscelltmp;
    elseif strcmp(celltype,'cell')
        clear iscelltmp
        data.iscell(1:size(tc,1),1) = 1;
    end
    
    ok = data.iscell(:,1); %Obtain iscell vector that index all the relevant ROIs
    
    TC = [TC tc(logical(ok),:)']; % Get tuning curve for each neuron and concatenate to TC vector
    
    neuronEachPlane(i) = sum(ok);

end

nNeuron = size(TC,2); %obtain total number of neurons
    
completeTC = nan(nFrames+100,nNeuron); % just so that if TC size is a lot larger, it can still work
completeTC(2:size(TC,1)+1,:) = TC; % shift one frame so that all the tone onset is on 101, 201, 301...
TC = completeTC(1:nFrames,:); % Take only the first 8500 frames which has tone evoked responses

%Find mean activity of all neurons
allNeuronMean = reshape(TC,nFramesPerTone,nTones,nTrials,nNeuron);
allNeuronMean = allNeuronMean(:,:,startTrial:end,:);
allNeuronMean = squeeze(nanmean(nanmean(nanmean(allNeuronMean,2),3),4));

%Find frames with peak average neuron response using smoothdata with a
%moving window (4 frames)
[~,peakIndex] = max(smoothdata(allNeuronMean,'movmean',4)); %argmax

if nPlanes==1
    startFrame = peakIndex-3;
    endFrame = peakIndex +3;
else
    startFrame = peakIndex-ceil(6/nPlanes);
    endFrame = peakIndex+ceil(6/nPlanes);
end
peakFrames  = startFrame : endFrame; %set frames to obtain peak activity from

%Plotting figure of mean activity and peak activities
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
%         [.2 .5 .6 .4],'String', sprintf('Confirm current frame selection'),'Fontsize',12);
%     h2 = uicontrol(p,'Units','normalized','Position',...
%         [.2 .1 .6 .4],'String',sprintf('Redo frame selection'),'Fontsize',12);
%     uiresume(gcbf);
%     UIFlag = false;
% 
%     uiwait(f);
%     close(f);
% end
 
pretoneFrames = length(peakFrames);
completeTCpretone = nan(nFrames+100,nNeuron); % just so that if TC size is a lot larger, it can still work
completeTCpretone(2+pretoneFrames:size(TC,1)+1+pretoneFrames,:) = TC; % shift one frame so that all the tone onset is on 101, 201, 301...
%Add 'pretone' frames to calculate local baseline
TCpretoneTemp = completeTCpretone(1:nFrames,:);

% Frequency of tones presented (always in the same order)
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

%%
TCreshape=reshape(TC,nFramesPerTone,nTones,nTrials,nNeuron);
%FC plot mean activity of neuron
% oneNeuron = squeeze(nanmean(TCreshape(:,:,:,1), 3)); 
% % mean over trials for one ROI
% figure; 
% for i = 1:17
%   subplot(3,6,i);
%   plot([1:50], oneNeuron(:,i)); axis([0 50 2200 3000]);
% end; 
% subplot(3,6,18); plot([1:50], squeeze(nanmean(oneNeuron, 2))); axis([0 50 2200 3000]);
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
baseline = prctile(reshape(TCreorder, [(nTrials-startTrial+1)*nFramesPerTone*nTones,nNeuron]),25); %Baseline is calculated using bottom 25%
TCreorder = TCreorder ./ repmat(reshape(baseline,[1 1 nNeuron]),nTrials-startTrial+1,nFramesPerTone*nTones,1);

%get mean response for each trial
TCtrialMean = squeeze(nanmean(TCreorder));
TCtrialMedian = squeeze(nanmedian(TCreorder));
TCtrialSEM = squeeze(nanstd(TCreorder)./sqrt(nTrials));

%get mean response for each tone
TCtoneMean = reshape(TCtrialMean,[nFramesPerTone,nTones,nNeuron]);
TCtoneMedian = reshape(TCtrialMedian,[nFramesPerTone,nTones,nNeuron]);

%get mean response for peak frames
tuningData = reshape(TCreorder, [nTrials-startTrial+1, nFramesPerTone,nTones,nNeuron]);
tuningData = reshape(tuningData(:,peakFrames,:,:),[(nTrials-startTrial+1)*length(peakFrames),nTones,nNeuron]);
tuningMean = squeeze(nanmean(tuningData));
tuningMedian = squeeze(nanmedian(tuningData));
tuningSEM = squeeze(nanstd(tuningData)/sqrt(size(tuningData,1))); %not n-1?

% first do a Anova to test if the cell is responsive to any tonesbn 
% need 18 groups, each group with 8 data points for 8 trials
% the baseline group is averaged from 8 trials 

peakActANOVA = zeros(nTrials-startTrial+1,nTones+1,nNeuron);
peakActANOVA(:,1,:) = squeeze(nanmean(nanmean(TCpretone(1:pretoneFrames,:,startTrial:nTrials,:))));
for i = 2:nTones+1
    peakActANOVA(:,i,:) = squeeze(mean(TCtone(peakFrames,i-1,startTrial:nTrials,:)));
end

% then do a t test to test if the cell is any-responsive, peak tone responsive vs. baseline
[~,peakIndex] = max(tuningMedian);
ttestPeakTone = zeros(2,nNeuron);
ttestAllTone = zeros(2,nNeuron);
anovaAllTone = zeros(1,nNeuron);
neuronRespTable = zeros(nTones+1, nNeuron +1);
neuronRespTable(:,1) = [sort(toneorder) 0];

% note fix the t-test here 
for i = 1:nNeuron
    peakAct = TCpretone(peakFrames+pretoneFrames,:,startTrial:nTrials,i);% i.e. 3,(8*17)
    peakAct = squeeze( nanmean(peakAct,1)); 
    baseAct = TCpretone(1:pretoneFrames,:,startTrial:nTrials,i);
    baseAct = squeeze( nanmean(baseAct,1));

    [h,p] = ttest(peakAct(:),baseAct(:));
    ttestAllTone(:,i) = [h;p]; 
        
    peakActPrefTone = squeeze(peakAct(peakIndex(i),:));
    baseActPrefTone = squeeze(baseAct(peakIndex(i),:));
    [h,p] = ttest(peakActPrefTone,baseActPrefTone);
    ttestPeakTone(:,i) = [h;p];

    groupNames = cell(1,nTones+1);
    groupNames{1} = 'Baseline';
    for j = 1:nTones
        groupNames{j+1} = int2str(toneorder(toneindex(j)));
    end
    [p,~,stats] = anova1(peakActANOVA(:,:,i),groupNames,'off');
    anovaAllTone (i) = p;
    
    [results,~,~,~] = multcompare(stats);
    results = results(results(:,1)==1,6);
    signifTone = (results<0.05);
    
    neuronRespTable(1:nTones,i+1) = signifTone;
    neuronRespTable(end,i+1) = p < 0.05;

end

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
save([suite2ppath '/' tuningFolderName '/population/TCtoneMean.mat'],'TCtoneMean'); % JL 05/08/19 for tuning and bandwidth analysis

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
% plot(freqAxis, popTuningPeakMedian,'LineWidth',2); hold on;
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
% plot(freqAxis, popTuningMean,'LineWidth',2);
set(gca, 'XTick', log2([4000 8000 16000 32000 64000]));
set(gca, 'XTickLabel', [4 8 16 32 64]);
xlim([log2(4000) log2(64000)]);
xlabel('Frequency (kHz)');
ylabel('mean F/F0');
legend('Median')
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
for i=1:nPlanes
    cd([suite2ppath 'plane' num2str(i-1)]);
    data = load('Fall.mat'); 
    refImg = data.ops.meanImg; % this is the image for the green channel not the red channel
    colormapIndex = round(linspace(1,64,17));
    C = colormap('jet');
    
    %select ROIs for different planes
    if i==1
        rois = roisCoord(1:neuronEachPlane(i));
    else
        rois = roisCoord(neuronEachPlane(i-1)+1:neuronPlane(end)); 
    end
    
    %create plots for population
    subplot(2,2,i);       
    imagesc(uint16(refImg));colormap gray;hold on;
    ylim([0 size(refImg,1)]);xlim([0 size(refImg,2)]);

    for j=(neuronPlane(i)+1):neuronPlane(i+1)
        if strcmp(file,'suite2p')
            x=rois{j-neuronPlane(i)}.xpix; %freehand rois have the outlines in x-y coordinates
            y=rois{j-neuronPlane(i)}.ypix; %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
        elseif strcmp(file,'handroi')
            x = rois{1,j-neuronPlane(i)}.mnCoordinates(:,1); %freehand rois have the outlines in x-y coordinates
            y = rois{1,j-neuronPlane(i)}.mnCoordinates(:,2);
        end
        
        if responsiveCellFlag(j)
            plot(x,y,'.','color',[0.9290, 0.6940, 0.1250],'MarkerSize',1);
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
currentNeuron = [0 neuronEachPlane(1)];
for j = 1:nPlanes
    cd([suite2ppath 'plane' num2str(j-1)]);
    data = load('Fall.mat'); 
    refImg = data.ops.meanImg;

    %select ROIs for different planes
    if j==1
        rois = roisCoord(1:neuronEachPlane(j));
    else
        rois = roisCoord(neuronEachPlane(j-1)+1:neuronPlane(end)); 
    end
    
    %     roiName = [ filePath '/' roiFolderName '/' fileName '_roi' int2str(j),suffix,'.zip'];
    %roiName = [ filePath '/' roiFolderName '/' fileName '_roi' int2str(j) '.zip'];
    %roiName = [filePath '/' roiFolderName '/' fileName '_ot_' num2str(whichPlane-1,'%03d') '_rois.zip'];
    %rois = ReadImageJROI(roiName);
    %cellThisPlane = length(rois);
    %startingIndex = sum(cellPerPlane);
    %cellPerPlane(j) = cellThisPlane;
    significantNeurons=[];
    for i = 1:neuronEachPlane(j)
        tic;
        
        %create plots for single neurons but not show them
        figure('visible','off');
        cellIndex = currentNeuron(j) + i;
        % REMEMBER TO CHANGE THIS LINE
        subplot(2,2,1)
        imagesc(refImg(:,:));colormap gray;hold on;

%         x=roisCoord{i+currentNeuron(j)}.xpix; %freehand rois have the outlines in x-y coordinates
%         y=roisCoord{i+currentNeuron(j)}.ypix; %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
        if strcmp(file,'suite2p')
            x=rois{i+currentNeuron(j)}.xpix; %freehand rois have the outlines in x-y coordinates
            y=rois{i+currentNeuron(j)}.ypix; %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
        elseif strcmp(file,'handroi')
            x = rois{1,i}.mnCoordinates(:,1); %freehand rois have the outlines in x-y coordinates
            y = rois{1,i}.mnCoordinates(:,2);
        end
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
        disp(signifTonePoint);
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
        
%         if ~isempty(signifTonePoint)
%             saveas(gcf,[suite2ppath '/' tuningFolderName ...
%                 '/singleNeuron/significant/Neuron' num2str(cellIndex,'%03d') '.png']);
%             saveas(gcf,[suite2ppath '/' tuningFolderName ...
%                 '/singleNeuron/significant/Neuron' num2str(cellIndex,'%03d') '.m']);
%             close all;
%             timeElapsed = toc;
%         else
            saveas(gcf,[suite2ppath '/' tuningFolderName ...
                '/singleNeuron/Neuron' num2str(cellIndex,'%03d') '.png']);
            saveas(gcf,[suite2ppath '/' tuningFolderName ...
                '/singleNeuron/Neuron' num2str(cellIndex,'%03d') '.m']);
%             close all;
            timeElapsed = toc;
%         end
        
        disp(['Cell #' int2str(cellIndex) ' Plane #' int2str(j) ' Time=' num2str(timeElapsed,'%03f') ' secs'])

    end
end



end