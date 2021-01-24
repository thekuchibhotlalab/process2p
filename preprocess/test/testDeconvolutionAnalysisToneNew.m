%% load the data from the whole traces
function testDeconvolutionAnalysisToneNew(mouse,day)
configpath = 'C:\Users\zzhu34\Documents\tempdata\excitData\config\mouse\';
datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mouse '\allSessions\'];
behavpath = ['C:\Users\zzhu34\Documents\tempdata\excitData\' mouse '\behavior'];
load(['C:\Users\zzhu34\Documents\tempdata\excitData\' mouse '\roi\ishere_plane0.mat'],'ishere');
imgSavePath = [datapath '\fig\']; mkdir(imgSavePath);
if iscell(ishere); ishere = ishere{1}; end
peakDetection = false; peakFrames = 1:6;

[nFrames, ~] = func_readnFrames(mouse,'root',configpath);
nFrames = [0 0; nFrames];
configTable =  readtable([configpath '\' mouse '_config.csv']);

filenameList = {[mouse '_calman_ar2_foo90_pars_allday_s.mat'],[mouse '_calman_ar2_foo90_pars_day_s.mat'],...
    [mouse '_calman_ar2_foo95_optb_nosmin_s.mat']};

for i = 1:length(filenameList); shortFilename{i} = filenameList{i}(24:end-6); shortFilename{i}(strfind(shortFilename{i},'_'))=' ';end
selectNeuron = ishere(:,day+1)==1;
tPeakMean = []; fPeakMean = [];
tPeakTrial = []; fPeakTrial = [];
tPeak =[]; fPeak = [];
sessionIdx = find(configTable.Day == (day) & ~strcmp(configTable.BehavType,'Baseline')); 
for i = 1:length(filenameList)
    disp([filenameList{i} ' started!'])
    load([datapath '\' filenameList{i}],'s');
    %tPeakMeanTemp = []; fPeakMeanTemp = [];
    tPeakTrialTemp = []; fPeakTrialTemp = [];
    tPeakTemp = []; fPeakTemp= [];
    for j = 1:length(sessionIdx)
        tempSession = sessionIdx(j);
        tempS = [nan(size(s{tempSession},1),1) s{tempSession}]; tempS = tempS(selectNeuron,:);      
        behavMatrix = importdata([behavpath '\' configTable.BehavFile{tempSession}]);
        [tempTPM, tempFPM,tempT,tempF,tempTPT,tempFPT] = getPeakAct(behavMatrix,tempS,peakDetection,peakFrames);
        %tPeakMeanTemp = cat(2,tPeakMeanTemp,tempTPM);fPeakMeanTemp = cat(2,fPeakMeanTemp,tempFPM) ;
        tPeakTrialTemp = cat(2,tPeakTrialTemp,tempTPT);fPeakTrialTemp = cat(2,fPeakTrialTemp,tempFPT);
        tPeakTemp = cat(3,tPeakTemp,tempT);fPeakTemp = cat(3,fPeakTemp,tempF);
    end
    tPeakMean = cat(2,tPeakMean,mean(tPeakTrialTemp,2));fPeakMean = cat(2,fPeakMean,mean(fPeakTrialTemp,2)) ;
    tPeakTrial = cat(3,tPeakTrial,tPeakTrialTemp);fPeakTrial = cat(3,fPeakTrial,fPeakTrialTemp);
    tPeak = cat(4,tPeak,tPeakTemp);fPeak = cat(4,fPeak,fPeakTemp);
end

% load and process dff data
disp('Dff started!');
try
    load(['C:\Users\zzhu34\Documents\tempdata\excitData\' mouse '\TC\' mouse '_TC_plane0.mat'],'tempTC');
    TC = tempTC;
catch 
    load(['C:\Users\zzhu34\Documents\tempdata\excitData\' mouse '\TC\' mouse '_TC_plane0_final.mat'],'TC');
end
%tPeakMeanDff = []; fPeakMeanDff = [];
tPeakTrialDff = []; fPeakTrialDff = [];
tPeakDff = []; fPeakDff= [];
peakDetection = true; peakFrames = 3:12;

selectDff = cell(1,length(sessionIdx));
for i = 1:length(sessionIdx)
    sessionTC = tempTC(:,(nFrames(sessionIdx(i),1)+1):nFrames(sessionIdx(i)+1,1));
    sessionDff = (sessionTC - repmat(median(sessionTC,2),1,size(sessionTC,2))) ./ repmat(median(sessionTC,2),1,size(sessionTC,2));
    sessionDff = sessionDff - smoothdata(sessionDff,2,'movmedian',3000);
    selectDff{i} = sessionDff(selectNeuron,2:end);
end


for j = 1:length(sessionIdx)
    tempSession = sessionIdx(j);
    tempS = [nan(size(selectDff{j},1),1) selectDff{j}]; 
    behavMatrix = importdata([behavpath '\' configTable.BehavFile{tempSession}]);
    [tempTPM, tempFPM,tempT,tempF,tempTPT,tempFPT] = getPeakAct(behavMatrix,tempS,peakDetection,peakFrames);
    %tPeakMeanDff = cat(2,tPeakMeanDff,tempTPM);fPeakMeanDff = cat(2,fPeakMeanDff,tempFPM) ;
    tPeakTrialDff = cat(2,tPeakTrialDff,tempTPT);fPeakTrialDff = cat(2,fPeakTrialDff,tempFPT);
    tPeakDff = cat(3,tPeakDff,tempT);fPeakDff = cat(3,fPeakDff,tempF);
end
tPeakMeanDff = mean(tPeakTrialDff,2);fPeakMeanDff = mean(fPeakTrialDff,2) ;

%% plot 1 - the averaged activity over days in tone-evoked period
f1 = figure; subplot(1,2,1);hold on; set(f1, 'Units', 'Normalized', 'OuterPosition',  [0.03, 0.4, 0.6, 0.5]);
dffFlat = [tPeakMeanDff(:)' fPeakMeanDff(:)'];
spikeFlat =  cat(1,tPeakMean,fPeakMean);
nbins = 10;
[N,edges,bin] = histcounts(dffFlat,linspace(0,prctile(dffFlat,95),nbins+1));
tempX = []; tempMean = []; tempSEM = [];
for j = 1:nbins
    tempX(j) = mean(edges(j:j+1));
    for i = 1:length(filenameList)
        tempMean(j,i) = nanmean(spikeFlat(bin==j,i)); tempSEM(j,i) = nanstd(spikeFlat(bin==j,i))/sqrt(N(j));     
    end
end
for i = 1:length(filenameList); plot(tempX,tempMean(:,i),'Color',matlabColors(i)); end
tempLegend = [];
hold on;
for i = 1:length(filenameList)
    f = fillErrorbarPlot(tempX,tempMean(:,i)', tempSEM(:,i)',matlabColors(i),'LineStyle','none');
    f.FaceAlpha = 0.1;
    scatter(dffFlat,spikeFlat(:,i),5,matlabColors(i), 'filled' , 'o', 'MarkerFaceAlpha', 0.9); 
end
ylabel('spike'); xlabel('dff'); 
xlim([0 prctile(dffFlat,98)]);
ylim([0 prctile(spikeFlat(:),99)]);
legend(shortFilename{:},'Location','Best');
title('Tone-evoked activity')

tPSTH = []; fPSTH = [];
binEdges = linspace(prctile(tPeakMeanDff,1),prctile(tPeakMeanDff,99),10); 
binSelection = 1:5;[N, edges, tbin] = histcounts(tPeakMeanDff,binEdges);
for i = 1:length(filenameList)
    tPSTH(:,:,i) = nanmean(tPeak(:,:,:,i),3);
    fPSTH(:,:,i) = nanmean(fPeak(:,:,:,i),3);
end

subplot(1,2,2); for i = 1:length(filenameList)
    plot(smoothdata(nanmean(tPSTH(:,:,i),1) - nanmean(nanmean(tPSTH(:,1:5,i),1),2),'gaussian',3)); hold on; end
saveas(gcf, [imgSavePath 'f1_day' int2str(day) '.png'])
%% plot 2.1 - PSTH of different groups of neurons by dff
figure;
for i = 1:length(filenameList)
    subplot(2,2,i);
    for j = 1:length(binSelection)
        tempAct = nanmean(tPSTH(tbin == binSelection(j) ,:,i),1) - nanmean(nanmean(tPSTH(tbin == binSelection(j),1:3,i),1),2);
        plot(tempAct); hold on; 
    end
    if i == 1; legend('dff 0.01-0.02','dff 0.03-0.04','dff 0.05-0.06','dff 0.07-0.08','dff 0.09-0.1');end
    title(shortFilename{i})
end
subplot(2,2,length(filenameList)+1);
tPSTHDff = nanmean(tPeakDff,3);
for j = 1:length(binSelection)
    tempAct = nanmean(tPSTHDff(tbin == binSelection(j) ,:),1) - nanmean(nanmean(tPSTHDff(tbin == binSelection(j),1:3),1),2);
    plot(tempAct); hold on; 
end
title('dff')
saveas(gcf, [imgSavePath 'f2.1_day' int2str(day) '.png'])
%% plot 2.2 - PSTH of different groups of neurons by percentile
binSelection = floor([.05 .2 .4 .6 .8 .95] * length(tPeakMeanDff));
[~, tSortIdx] = sort(tPeakMeanDff,'ascend');
tPSTH = []; fPSTH = [];
for i = 1:length(filenameList)
    tPSTH(:,:,i) = nanmean(tPeak(:,:,:,i),3);
    fPSTH(:,:,i) = nanmean(fPeak(:,:,:,i),3);
end
tPSTHDff = nanmean(tPeakDff,3); fPSTHDff = nanmean(fPeakDff,3);

figure;
for i = 1:length(filenameList)
    subplot_tight(2,2,i);
    for j = 1:length(binSelection)-1
        tempIdx= tSortIdx(binSelection(j): binSelection(j+1)-1);
        tempAct = nanmean(tPSTH(tempIdx,:,i),1) - nanmean(nanmean(tPSTH(tempIdx,1:3,i),1),2);
        plot(tempAct); hold on; 
    end
    if i == 1; legend('5-20 perc','20-40 perc','40-60 perc','60-80 perc','80-95 perc');end
    title(shortFilename{i})
end
subplot_tight(2,2,length(filenameList)+1); title('dff')
for j = 1:length(binSelection)-1
    tempIdx= tSortIdx(binSelection(j): binSelection(j+1)-1);
    tempAct = nanmean(tPSTHDff(tempIdx ,:),1) - nanmean(nanmean(tPSTHDff(tempIdx,1:3),1),2);
    plot(tempAct); hold on; 
end
saveas(gcf, [imgSavePath 'f2.2_day' int2str(day) '.png'])
%% plot 2.3 - Correlation of peak activity 

actMatT = cat(2,tPeakMean, tPeakMeanDff);
actMatF = cat(2,fPeakMean, fPeakMeanDff);
actMat = cat(1,actMatT,actMatF);
f3 = figure; set(f3, 'Units', 'Normalized', 'OuterPosition',  [0.03, 0.4, 0.6, 0.5]);
subplot(1,2,1);imagesc(corr(actMat)); colorbar; caxis([0.6 1]); title('Correlation of Mean Activity')

actMatT = cat(2,reshape(tPeakTrial,size(tPeakTrial,1)*size(tPeakTrial,2),size(tPeakTrial,3)), tPeakTrialDff(:));
actMatF = cat(2,reshape(fPeakTrial,size(fPeakTrial,1)*size(fPeakTrial,2),size(fPeakTrial,3)), fPeakTrialDff(:));
actMat = cat(1,actMatT,actMatF);
subplot(1,2,2); imagesc(corr(actMat)); colorbar; caxis([0.3 1]); title('Correlation of Trial-by-Trial Activity ')
saveas(gcf, [imgSavePath 'f2.3_day' int2str(day) '.png'])
%% plot 3 - get number of responsive neurons 
tAct = tPeak;fAct = fPeak;
peakDetection = false;peakFrames = 1:6;smoothWindow = 0;
for i = 1:length(filenameList)
    for j = 1:size(tAct,1)
        tempAct = squeeze(tAct(j,:,:,i)); [h, tempAuc, tempIdx] = testResponsive(tempAct,peakDetection,peakFrames,smoothWindow);
        ttestT(j,i) = h;rocAucT(j,i) = tempAuc; peakIdxT(j,i) = tempIdx;
        tempAct = squeeze(fAct(j,:,:,i)); [h, tempAuc] = testResponsive(tempAct,peakDetection,peakFrames,smoothWindow);
        ttestF(j,i) = h;rocAucF(j,i) = tempAuc;
    end
end
peakDetection = true;peakFrames = 3:12;smoothWindow = 3;
for j = 1:size(tPeakDff,1)
    tempAct = squeeze(tPeakDff(j,:,:)); [h, tempAuc, tempIdx] = testResponsive(tempAct,peakDetection,peakFrames,smoothWindow);
    ttestTDff(j) = h; rocAucTDff(j) = tempAuc; peakIdxT(j,i) = tempIdx;
    tempAct = squeeze(fPeakDff(j,:,:)); [h, tempAuc] = testResponsive(tempAct,peakDetection,peakFrames,smoothWindow);
    ttestFDff(j) = h; rocAucFDff(j) = tempAuc;
end

%% plot 3.1 - pie plots
figure;
for i = 1:length(filenameList)
    h = subplot(1,3,i);
    bothFlag = ((ttestF(:,i)==1)' & ttestFDff==1);
    dffFlag = ((ttestF(:,i)==0)' & ttestFDff==1);
    spkFlag = ((ttestF(:,i)==1)' & ttestFDff==0);
    fn_plotPieMultiPanel({{bothFlag,dffFlag,spkFlag}}, 'legend', {{'both','dff only','spike only'}},...
        'title',shortFilename(i),'axes',h);
end
saveas(gcf, [imgSavePath 'f3.1_day' int2str(day) '.png'])
%% plot 3.2 - pie plots
figure; subplot(2,2,1);hold on;
for i = 1:length(filenameList)
    bins = 0.4:0.005:1.0; binsPlot = bins(1:end-1) + 0.005/2;
    %[N, ~,~] = histcounts(rocAucT(:,i), 0.4:0.005:0.8,'Normalization','probability');
    cdfplot([rocAucT(:,i);rocAucF(:,i)]); xlim([0.5 1.0])
end
h = cdfplot([rocAucTDff rocAucFDff] ); xlim([0.5 1.0])
set( h, 'Color', [0 0 0], 'LineWidth', 1.5); legend({shortFilename{:},'dff'},'Location','Best')
title('Distribution of ROC value');xlabel('AUC value of ROC')

subplot(2,2,2);  hold on;
respFlag = ttestTDff == 1 | ttestFDff == 1;
h = cdfplot([rocAucT(respFlag,3);rocAucF(respFlag,3)]); set( h, 'Color', matlabColors(1), 'LineWidth', 1.5);
h = cdfplot([rocAucT(~respFlag,3);rocAucF(~respFlag,3)]);set( h, 'Color', matlabColors(1),'LineStyle','--', 'LineWidth', 1.5);
h = cdfplot([rocAucTDff(respFlag) rocAucFDff(respFlag)]); set( h, 'Color', [0 0 0], 'LineWidth', 1.5);
h = cdfplot([rocAucTDff(~respFlag) rocAucFDff(~respFlag)]); set( h, 'Color', [0 0 0], 'LineStyle','--', 'LineWidth', 1.5);
title('Target & Foil ROC value'); legend('spk-resp','spk-not resp','dff-resp','dff-not resp','Location','Best')
xlabel('AUC value of ROC');

subplot(2,2,3); hold on; 
respFlag = (ttestTDff' == 1 & sum(ttestT,2) == 3); notRespFlag = (ttestTDff' == 1 & sum(ttestT,2) == 0);
h = cdfplot(mean(rocAucT(respFlag,:),2)); set( h, 'Color', matlabColors(1), 'LineWidth', 1.5);
h = cdfplot(mean(rocAucT(notRespFlag,:),2));set( h, 'Color', matlabColors(1),'LineStyle','--', 'LineWidth', 1.5);
h = cdfplot(rocAucTDff(respFlag)); set( h, 'Color', [0 0 0], 'LineWidth', 1.5);
h = cdfplot(rocAucTDff(notRespFlag)); set( h, 'Color', [0 0 0], 'LineStyle','--', 'LineWidth', 1.5);
legend('spk-resp','spk-not resp','dff-resp','dff-not resp','Location','Best')
xlabel('AUC value of ROC');title('Targ Resp Neuron ROC')

subplot(2,2,4); hold on; 
respFlag = (ttestFDff' == 1 & sum(ttestF,2) == 3); notRespFlag = (ttestFDff' == 1 & sum(ttestF,2) == 0);
h = cdfplot(mean(rocAucF(respFlag,:),2)); set( h, 'Color', matlabColors(1), 'LineWidth', 1.5);
h = cdfplot(mean(rocAucF(notRespFlag,:),2));set( h, 'Color', matlabColors(1),'LineStyle','--', 'LineWidth', 1.5);
h = cdfplot(rocAucFDff(respFlag)); set( h, 'Color', [0 0 0], 'LineWidth', 1.5);
h = cdfplot(rocAucFDff(notRespFlag)); set( h, 'Color', [0 0 0], 'LineStyle','--', 'LineWidth', 1.5);
legend('spk-resp','spk-not resp','dff-resp','dff-not resp','Location','Best')
xlabel('AUC value of ROC');title('Foil Resp Neuron ROC')
saveas(gcf, [imgSavePath 'f3.2_day' int2str(day) '.png'])
%% plot 4 - calculate the SI of all neurons
spikeSI = abs(tPeakMean- fPeakMean) ./ (abs(tPeakMean) +  abs(fPeakMean));

spikeSI(isnan(spikeSI)) = 0;
spikeSIDff = abs(tPeakMeanDff- fPeakMeanDff) ./ (abs(tPeakMeanDff) +  abs(fPeakMeanDff));

f4 = figure; set(f4, 'Units', 'Normalized', 'OuterPosition',  [0.03, 0.4, 0.9, 0.5]);
subplot(1,3,1); hold on;
for i = 1:length(filenameList); cdfplot(spikeSI(:,i));end; h = cdfplot(spikeSIDff);
set( h, 'Color', [0 0 0], 'LineWidth', 1.5); legend({shortFilename{:},'dff'},'Location','Best');
title('Distribution of SI')
subplot(1,3,2); hold on;
for i = 1:length(filenameList); cdfplot(spikeSI(:,i)-spikeSIDff);end
legend(shortFilename,'Location','Best');title('Distribution of deltaSI')
subplot(1,3,3); 
tempMat = cat(2,spikeSI,spikeSIDff);imagesc(corr(tempMat));
title('Correlation of SI'); colorbar; caxis([0.4 1])
saveas(gcf, [imgSavePath 'f4_day' int2str(day) '.png'])
%% plot 5 - PSTH plotted for all neurons
tPSTH = []; fPSTH = [];
tPSTHDff = nanmean(tPeakDff,3); fPSTHDff = nanmean(fPeakDff,3);
f5 = figure; set(f5, 'Units', 'Normalized', 'OuterPosition',  [0.02, 0.3, 0.95, 0.6]);
subplot_tight(1+length(filenameList),1,1);
dffTrace = [ reshape(tPSTHDff',numel(tPSTHDff),1)' reshape(fPSTHDff',numel(fPSTHDff),1)']; plot(dffTrace);
xlim([1 length(dffTrace)]); title('dff')
for i = 1:length(filenameList)
    tPSTH(:,:,i) = nanmean(tPeak(:,:,:,i),3);
    fPSTH(:,:,i) = nanmean(fPeak(:,:,:,i),3);
    subplot_tight(1+length(filenameList),1,1+i); hold on;
    tempTrace = [reshape(tPSTH(:,:,i)',numel(tPSTH(:,:,i)),1)' reshape(fPSTH(:,:,i)',numel(fPSTH(:,:,i)),1)'];
    tempReg = [dffTrace;ones(size(dffTrace))]; coeff = tempReg' \ tempTrace'; 
    
    plot(tempTrace);xlim([1 length(tempTrace)]); 
    plot(tempReg'*coeff,'Color',matlabColors(2) * 0.5 + [1 1 1]*0.5);
    ylim([prctile(tempTrace,0.1) prctile(tempTrace,99.8)])
    title(shortFilename{i});
end
saveas(gcf, [imgSavePath 'f5_day' int2str(day) '.png'])
close all;
end



%% functions 
function [h, tempAuc, peakIdx] = testResponsive(tempAct,peakDetection,peakFrames,smoothWindow)
    tempBaseline = tempAct(1:3,:); tempBaselineMean = mean(tempBaseline,1);
    if smoothWindow ~= 0; tempMeanAct = smoothdata(mean(tempAct,2),'gaussian',smoothWindow);
    else; tempMeanAct = mean(tempAct,2); end
    toneFrame = 5; % number of pretone frames
    selectToneFrame = toneFrame + peakFrames;
    [~,tempPeakIdx] = max(tempMeanAct(selectToneFrame));
    peakIdx = tempPeakIdx + selectToneFrame(1)-1; 
    if peakDetection; tempPeak = tempAct(peakIdx,:); 
    else; tempPeak = mean(tempAct(selectToneFrame,:),1); end
    [h,p] = ttest(tempPeak,tempBaselineMean,'tail','right');
    % do an roc for each tone 
    rocAct = [tempBaseline(:)' tempPeak]; 
    rocAct(rocAct<1e-4) = randn(1,sum(rocAct<1e-4))*1e-4;% introduce small noise to spks to avoid domination of 0s
    rocLabel = [zeros(1,numel(tempBaseline)) ones(1,length(tempPeak))];
    [tpr, fpr, threshold] = roc(rocLabel, rocAct);
    tempAuc = trapz([0 fpr 1],[0 tpr 1]);
end

function [tPeakMean, fPeakMean, tAct, fAct,tPeakTrial,fPeakTrial] = getPeakAct(behavMatrix,act,peakDetection,peakFrames)
    tFrame = floor(behavMatrix(behavMatrix(:,4)==1 | behavMatrix(:,4)==2,12)/2);
    fFrame = floor(behavMatrix(behavMatrix(:,4)==3 | behavMatrix(:,4)==4,12)/2);
    selectFrame = -5 : 30; toneFrame = abs(selectFrame(1)); selectToneFrame = toneFrame + peakFrames;
    tAct = []; fAct = [];
    for k = 1:length(tFrame)
        tAct(:,:,k) = act(:,tFrame(k) + selectFrame);
    end
    for k = 1:length(fFrame)
        fAct(:,:,k) = act(:,fFrame(k) + selectFrame);
    end
    if peakDetection
        tActMean = mean(tAct,3); tPeakMean = max(smoothdata(tActMean(:,selectToneFrame),2,'gaussian',3),[],2);
        fActMean = mean(fAct,3); fPeakMean = max(smoothdata(fActMean(:,selectToneFrame),2,'gaussian',3),[],2);
        tPeakTrial = squeeze(max(smoothdata(tAct(:,selectToneFrame,:),2,'gaussian',3),[],2));
        fPeakTrial = squeeze(max(smoothdata(fAct(:,selectToneFrame,:),2,'gaussian',3),[],2));
        
    else
        tPeakTrial = squeeze(mean(tAct(:,selectToneFrame,:),2));tPeakMean = mean(tPeakTrial,2);
        fPeakTrial = squeeze(mean(fAct(:,selectToneFrame,:),2));fPeakMean = mean(fPeakTrial,2);
    end
end

function f = fillErrorbarPlot(xdata,yMean, ySEM,varargin)
    x = [xdata, fliplr(xdata)];
    y = [yMean+ySEM, fliplr(yMean-ySEM)];
    f = fill(x,y,varargin{:});
end