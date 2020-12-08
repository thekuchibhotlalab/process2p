clear;
datapath ='C:\Users\zzhu34\Documents\tempdata\deconv_test_cd036\smooth3000\behavior';
behavpath = 'C:\Users\zzhu34\Documents\tempdata\excitData\cd036\behavior';
%filenameList = {'cd036_confoo_smin.mat','cd036_foo_smin.mat','cd036_thre_smin.mat'};
%filenameList = {'cd036_calman_confoo_optb.mat','cd036_calman_foo90_optb_nosmin.mat',...
%    'cd036_calman_foo95_optb_nosmin.mat','cd036_calman_foo95_optb.mat','cd036_calman_foo99_optb.mat','cd036_calman_thre95_optb.mat'};
filenameList = {'cd036_calman_confoo_optb.mat','cd036_calman_foo90_optb_nosmin.mat',...
    'cd036_calman_foo95_optb_nosmin.mat','cd036_calman_foo95_optb.mat','cd036_calman_foo99_optb.mat','cd036_calman_ar2_foo90_optb_nosmin.mat'};
%filenameList = {'cd036_foo_lambda0.mat','cd036_foo_lambda05.mat','cd036_foo_lambda10.mat'};
shortFilename = [];
global peakDetection peakFrames;
peakDetection = true; peakFrames = 1:6;
for i = 1:length(filenameList); shortFilename{i} = filenameList{i}(14:end-4); shortFilename{i}(strfind(shortFilename{i},'_'))=' ';end
for i = 1:length(filenameList)
    disp([filenameList{i} ' started!'])
    load([datapath '\' filenameList{i}],'s','c','day');
    for j = 1:length(day)
        behavMatrix = importdata([behavpath '\cd036_' int2str(day(j)-1) 'v1.txt']);
        [tempTPM, tempFPM,tempT,tempF] = getPeakAct(behavMatrix,s{j});
        tPeakMean(:,j,i) = tempTPM;fPeakMean(:,j,i) = tempFPM;
        tPeak(:,:,:,j,i) = tempT;fPeak(:,:,:,j,i) = tempF;
    end
end
load([datapath '\selectDff.mat'],'selectDff');
disp('selected dff started!')
for j = 1:length(day)
    behavMatrix = importdata([behavpath '\cd036_' int2str(day(j)-1) 'v1.txt']);
    [tempTP, tempFP, tempT, tempF] = getPeakAct(behavMatrix,selectDff{j});
    tPeakMeanDff(:,j) = tempTP;fPeakMeanDff(:,j) = tempFP;
    tPeakDff(:,:,:,j) = tempT;fPeakDff(:,:,:,j) = tempF;
end
%% plot 1 - the averaged activity over days in tone-evoked period
figure; hold on;
dffFlat = [tPeakMeanDff(:)' fPeakMeanDff(:)'];
for i = 1:length(filenameList); spikeFlat(:,i) = [reshape(tPeakMean(:,:,i),1,[]) reshape(fPeakMean(:,:,i),1,[])];end
[N,edges,bin] = histcounts(dffFlat,linspace(0,prctile(dffFlat,98),21));
for j = 1:20 
    tempX(j) = mean(edges(j:j+1));
    for i = 1:length(filenameList)
        tempMean(j,i) = mean(spikeFlat(bin==j,i)); tempSEM(j,i) = std(spikeFlat(bin==j,i))/sqrt(N(j));     
    end
end
for i = 1:length(filenameList); plot(tempX,tempMean(:,i),'Color',matlabColors(i)); end
tempLegend = [];
for i = 1:length(filenameList)
    f = fillErrorbarPlot(tempX,tempMean(:,i)', tempSEM(:,i)',matlabColors(i),'LineStyle','none');
    f.FaceAlpha = 0.1;
    scatter(dffFlat,spikeFlat(:,i),5,matlabColors(i), 'filled' , 'o', 'MarkerFaceAlpha', 0.1); 
end
ylabel('spike'); xlabel('dff'); 
xlim([0 prctile(dffFlat,98)]);
ylim([0 prctile(spikeFlat(:),99)]);
legend(shortFilename{:},'Location','Best');
title('Tone-evoked activity')
%% plot 1.1 - the averaged activity over days in tone-evoked period
figure; 
for j = 1:length(filenameList); subplot_tight(2,3,j,[0.15,0.06]);hold on; title(shortFilename{j})
for i= 1:7
    scatter([tPeakMeanDff(:,i)' fPeakMeanDff(:,i)'],...
        [reshape(tPeakMean(:,i,j),1,[]) reshape(fPeakMean(:,i,j),1,[])],...
        5, 'filled' , 'o','MarkerFaceAlpha', 0.4);
end
ylabel('spike'); xlabel('dff'); 
xlim([0 prctile([tPeakMeanDff(:)' fPeakMeanDff(:)'],99)]);
ylim([0 prctile([reshape(tPeakMean(:,:,j),1,[]) reshape(fPeakMean(:,:,j),1,[])],99)]);
end

%% plot 2.1 - PSTH of different groups of neurons by dff
binEdges = 0:0.01:0.15; binSelection = [2 4 6 8 10];
tPeakMeanDffAvg = mean(tPeakMeanDff,2); [N, edges, tbin] = histcounts(tPeakMeanDffAvg,binEdges);
for i = 1:length(filenameList)
    tPSTH(:,:,i) = mean(reshape(tPeak(:,:,:,:,i),size(tPeak,1), size(tPeak,2),[]),3);
    fPSTH(:,:,i) = mean(reshape(fPeak(:,:,:,:,i),size(fPeak,1), size(fPeak,2),[]),3);
end
tDffPSTH = mean(reshape(tPeakDff,size(tPeakDff,1), size(tPeakDff,2),[]),3);
fDffPSTH = mean(reshape(fPeakDff,size(fPeakDff,1), size(fPeakDff,2),[]),3);
figure; for i = 1:length(binSelection)
    plot(smoothdata(mean(tPSTH(:,:,i),1) - mean(mean(tPSTH(:,1:5,i),1),2),'gaussian',3)); hold on; end
figure;
for i = 1:length(filenameList)
    subplot_tight(3,3,i);
    for j = 1:length(binSelection)
        tempAct = mean(tPSTH(tbin == binSelection(j) ,:,i),1) - mean(mean(tPSTH(tbin == binSelection(j),1:3,i),1),2);
        plot(tempAct); hold on; 
    end
    if i == 1; legend('dff 0.01-0.02','dff 0.03-0.04','dff 0.05-0.06','dff 0.07-0.08','dff 0.09-0.1');end
    title(shortFilename{i})
end
subplot_tight(3,3,length(filenameList)+1); title('dff')
for j = 1:length(binSelection)
    tempAct = mean(tDffPSTH(tbin == binSelection(j) ,:),1) - mean(mean(tDffPSTH(tbin == binSelection(j),1:3),1),2);
    plot(tempAct); hold on; 
end

%% plot 2.2 - PSTH of different groups of neurons by percentile
binEdges = 0:0.01:0.15; binSelection = floor([.05 .2 .4 .6 .8 .95] * length(tPeakMeanDffAvg));
tPeakMeanDffAvg = mean(tPeakMeanDff,2); [~, tSortIdx] = sort(tPeakMeanDffAvg,'ascend');
for i = 1:length(filenameList)
    tPSTH(:,:,i) = mean(reshape(tPeak(:,:,:,:,i),size(tPeak,1), size(tPeak,2),[]),3);
    fPSTH(:,:,i) = mean(reshape(fPeak(:,:,:,:,i),size(fPeak,1), size(fPeak,2),[]),3);
end
tDffPSTH = mean(reshape(tPeakDff,size(tPeakDff,1), size(tPeakDff,2),[]),3);
fDffPSTH = mean(reshape(fPeakDff,size(fPeakDff,1), size(fPeakDff,2),[]),3);
figure;
for i = 1:length(filenameList)
    subplot_tight(3,3,i);
    for j = 1:length(binSelection)-1
        tempIdx= tSortIdx(binSelection(j): binSelection(j+1)-1);
        tempAct = mean(tPSTH( tempIdx,:,i),1) - mean(mean(tPSTH(tempIdx,1:3,i),1),2);
        plot(tempAct); hold on; 
    end
    if i == 1; legend('5-20 perc','20-40 perc','40-60 perc','60-80 perc','80-95 perc');end
    title(shortFilename{i})
end
subplot_tight(3,3,length(filenameList)+1); title('dff')
for j = 1:length(binSelection)-1
    tempIdx= tSortIdx(binSelection(j): binSelection(j+1)-1);
    tempAct = mean(tDffPSTH(tempIdx ,:),1) - mean(mean(tDffPSTH(tempIdx,1:3),1),2);
    plot(tempAct); hold on; 
end

%% plot 3 - get number of responsive neurons 
tAct = reshape(tPeak,size(tPeak,1),size(tPeak,2),size(tPeak,3)*size(tPeak,4),size(tPeak,5));
fAct = reshape(fPeak,size(fPeak,1),size(fPeak,2),size(fPeak,3)*size(fPeak,4),size(fPeak,5));
for i = 1:length(filenameList)
    for j = 1:size(tAct,1)
        tempAct = squeeze(tAct(j,:,:,i)); [h, tempAuc, tempIdx] = testResponsive(tempAct);
        ttestT(j,i) = h;rocAucT(j,i) = tempAuc; peakIdxT(j,i) = tempIdx;
        tempAct = squeeze(fAct(j,:,:,i)); [h, tempAuc] = testResponsive(tempAct);
        ttestF(j,i) = h;rocAucF(j,i) = tempAuc;
    end
end
tActDff = reshape(tPeakDff,size(tPeakDff,1),size(tPeakDff,2),[]);fActDff = reshape(fPeakDff,size(fPeakDff,1),size(fPeakDff,2),[]);
for j = 1:size(tActDff,1)
    tempAct = squeeze(tActDff(j,:,:)); [h, tempAuc, tempIdx] = testResponsive(tempAct);
    ttestTDff(j) = h; rocAucTDff(j) = tempAuc; peakIdxT(j,i) = tempIdx;
    tempAct = squeeze(fActDff(j,:,:)); [h, tempAuc] = testResponsive(tempAct);
    ttestFDff(j) = h; rocAucFDff(j) = tempAuc;
end

%% plot 3.1 - 
figure;
for i = 1:length(filenameList)
    h = subplot(2,3,i);
    bothFlag = ((ttestT(:,i)==1)' & ttestTDff==1);
    dffFlag = ((ttestT(:,i)==0)' & ttestTDff==1);
    spkFlag = ((ttestT(:,i)==1)' & ttestTDff==0);
    func_piePlot({{bothFlag,dffFlag,spkFlag}}, 'legend', {{'both','dff only','spike only'}},...
        'title',shortFilename(i),'axes',h);
end
figure; hold on;
for i = 1:length(filenameList)
    bins = 0.4:0.005:0.8; binsPlot = bins(1:end-1) + 0.005/2;
    %[N, ~,~] = histcounts(rocAucT(:,i), 0.4:0.005:0.8,'Normalization','probability');
    cdfplot(rocAucT(:,i)); xlim([0.5 0.8])
end
h = cdfplot(rocAucTDff); xlim([0.5 0.8])
set( h, 'Color', [0 0 0], 'LineWidth', 1.5); legend({shortFilename{:},'dff'},'Location','Best')
title('Distribution of ROC value')
%% plot 4 - calculate the SI of all neurons
spikeSI = abs(tPeakMean- fPeakMean) ./ (abs(tPeakMean) +  abs(fPeakMean));
spikeSI = reshape(spikeSI,size(spikeSI,1)*size(spikeSI,2),[]);
spikeSI(isnan(spikeSI)) = 0;
spikeSIDff = abs(tPeakMeanDff- fPeakMeanDff) ./ (abs(tPeakMeanDff) +  abs(fPeakMeanDff));
spikeSIDff = spikeSIDff(:);
figure; hold on;
for i = 1:length(filenameList); cdfplot(spikeSI(:,i));end; h = cdfplot(spikeSIDff);
set( h, 'Color', [0 0 0], 'LineWidth', 1.5); legend({shortFilename{:},'dff'},'Location','Best');
title('Distribution of SI')
%% plot 4.1 - the histogram of SI that 
%% functions 
function [h, tempAuc, peakIdx] = testResponsive(tempAct)
    global peakDetection peakFrames;
    tempBaseline = tempAct(1:5,:); tempBaselineMean = mean(tempBaseline,1);
    tempMeanAct = smoothdata(mean(tempAct,2),'gaussian',3);
    toneFrame = 6; selectToneFrame = toneFrame + peakFrames;
    [~,tempPeak] = max(tempMeanAct(selectToneFrame));
    tempPeak = tempPeak + selectToneFrame(1)-1;  peakIdx = tempPeak;
    if peakDetection; tempPeak = tempAct(tempPeak,:); 
    else; tempPeak = mean(tempAct(selectToneFrame,:),1); end
    [h,p] = ttest(tempPeak,tempBaselineMean);  
    % do an roc for each tone 
    rocAct = [tempBaseline(:)' tempPeak];
    rocLabel = [zeros(1,numel(tempBaseline)) ones(1,length(tempPeak))];
    [tpr, fpr, threshold] = roc(rocLabel, rocAct);
    tempAuc = trapz([0 fpr 1],[0 tpr 1]);
end

function [tPeakMean, fPeakMean, tAct, fAct] = getPeakAct(behavMatrix,act)
    global peakDetection peakFrames;
    tFrame = ceil(behavMatrix(behavMatrix(:,4)==1 | behavMatrix(:,4)==2,12)/2);
    fFrame = ceil(behavMatrix(behavMatrix(:,4)==3 | behavMatrix(:,4)==4,12)/2);
    selectFrame = -5 : 10; toneFrame = abs(selectFrame(1))+1; selectToneFrame = toneFrame + peakFrames;
    tAct = []; fAct = [];
    for k = 1:length(tFrame)
        tAct(:,:,k) = act(:,tFrame(k) + selectFrame);
        fAct(:,:,k) = act(:,fFrame(k) + selectFrame);
    end
    if peakDetection
        tActMean = mean(tAct,3); tPeakMean = max(smoothdata(tActMean(:,selectToneFrame,:),2,'gaussian',3),[],2);
        fActMean = mean(fAct,3); fPeakMean = max(smoothdata(fActMean(:,selectToneFrame,:),2,'gaussian',3),[],2);
        tPeakTrial = squeeze(max(smoothdata(tAct(:,selectToneFrame,:),2,'gaussian',3),[],2));
        fPeakTrial = squeeze(max(smoothdata(fAct(:,selectToneFrame,:),2,'gaussian',3),[],2));
    else
        tActMean = mean(tAct,3); tPeakMean = squeeze(mean(tActMean(:,selectToneFrame,:),2));
        fActMean = mean(fAct,3); fPeakMean = squeeze(mean(fActMean(:,selectToneFrame,:),2));
    end
    
    
end

function f = fillErrorbarPlot(xdata,yMean, ySEM,varargin)
    x = [xdata, fliplr(xdata)];
    y = [yMean+ySEM, fliplr(yMean-ySEM)];
    f = fill(x,y,varargin{:});
end