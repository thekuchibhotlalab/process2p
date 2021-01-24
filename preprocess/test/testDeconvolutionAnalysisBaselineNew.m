function testDeconvolutionAnalysisBaselineNew(mouse,selectday)

configpath = 'C:\Users\zzhu34\Documents\tempdata\excitData\config\mouse\';
datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mouse '\allSessions\'];
imgSavePath = [datapath '\fig_baseline\']; mkdir(imgSavePath);
load(['C:\Users\zzhu34\Documents\tempdata\excitData\' mouse '\roi\ishere_plane0.mat'],'ishere');
if iscell(ishere); ishere = ishere{1}; end

filenameList = {[mouse '_calman_ar2_foo90_pars_allday_s.mat'],[mouse '_calman_ar2_foo90_pars_day_s.mat'],...
    [mouse '_calman_ar2_foo95_optb_nosmin_s.mat']};
for i = 1:length(filenameList); shortFilename{i} = filenameList{i}(24:end-6); shortFilename{i}(strfind(shortFilename{i},'_'))=' ';end

[nFrames, ~] = func_readnFrames(mouse,'root',configpath);
nFrames = [0 0; nFrames];
configTable =  readtable([configpath '\' mouse '_config.csv']);

plotNeuron = true;

selectNeuron = ishere(:,selectday+1)==1; alls = cell(1,length(filenameList));
sessionIdx = find(configTable.Day == (selectday) &  strcmp(configTable.BehavType,'Baseline')); 
for i = 1:length(filenameList)
    disp([filenameList{i} ' started!'])
    load([datapath '\' filenameList{i}],'s');
    tempS =  s{sessionIdx}(selectNeuron,:);
    alls{i} = tempS;
end

disp('Dff started!');
try
    load(['C:\Users\zzhu34\Documents\tempdata\excitData\' mouse '\TC\' mouse '_TC_plane0.mat'],'tempTC');
    TC = tempTC;
catch 
    load(['C:\Users\zzhu34\Documents\tempdata\excitData\' mouse '\TC\' mouse '_TC_plane0_final.mat'],'TC');
end
sessionTC = tempTC(:,(nFrames(sessionIdx,1)+1):nFrames(sessionIdx+1,1));
sessionDff = (sessionTC - repmat(median(sessionTC,2),1,size(sessionTC,2))) ./ repmat(median(sessionTC,2),1,size(sessionTC,2));
sessionDff = sessionDff - smoothdata(sessionDff,2,'movmedian',3000);
selectDff = sessionDff(selectNeuron,2:end); 

% for cd036, also plot the ar1 models
if strcmp(mouse,'cd036') 
    extraAlls = [];
    extraDatapath ='C:\Users\zzhu34\Documents\tempdata\deconv_test_cd036\smooth3000\baseline';
    extraFilenameList = {'cd036_calman_ar2_confoo_optb.mat','cd036_calman_ar2_foo95_optb_nosmin.mat'}; 
    extraShortFilename = []; for i = 1:length(extraFilenameList)
        extraShortFilename{i} = extraFilenameList{i}(14:end-4); extraShortFilename{i}(strfind(extraShortFilename{i},'_'))=' ';end
    for i = 1:length(extraFilenameList)
        disp([extraFilenameList{i} ' started!'])
        load([extraDatapath '\' extraFilenameList{i}],'s','c','day');
        if any((selectday+1)==day); extraAlls{i} = s{(selectday+1)==day};end
    end
end
%% plot 1 - the averaged activity over days in tone-evoked period
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition',  [0.03, 0.4, 0.6, 0.5]);
spikeCount = []; ISIcount = [];
for i = 1:length(alls)
    for j = 1:size(alls{i},1)

        spikes = find(alls{i}(j,:)>1e-3) ;
        ISI = diff(spikes);
        [N,edges,bins] = histcounts(ISI,1:100,'Normalization','probability');
        ISIcount{i}(:,j) = N;
        spikeCount{i}(j) = length(spikes);

    end
    spikeCountAvg = spikeCount{i};ISIcountAvg = nanmean(ISIcount{i},2);
    subplot(1,2,1); cdfplot(spikeCountAvg); hold on;
    subplot(1,2,2); plot(ISIcountAvg); hold on;
end
tempNames =shortFilename ;
if strcmp(mouse,'cd036') && ~isempty(extraAlls)
    spikeCount = []; ISIcount = [];
    for i = 1:length(extraAlls)
        for j = 1:size(extraAlls{i},1)

            spikes = find(extraAlls{i}(j,:)>1e-3) ;
            ISI = diff(spikes);
            [N,edges,bins] = histcounts(ISI,1:100,'Normalization','probability');
            ISIcount{i}(:,j) = N;
            spikeCount{i}(j) = length(spikes);

        end
    spikeCountAvg = spikeCount{i};ISIcountAvg = nanmean(ISIcount{i},2);
    subplot(1,2,1); cdfplot(spikeCountAvg); hold on;
    subplot(1,2,2); plot(ISIcountAvg); hold on;
    end
    tempNames =[tempNames extraShortFilename];
end
subplot(1,2,1); title('spike count distribution'); legend(tempNames{:},'Location','Best');
subplot(1,2,2); title('spike ISI distribution'); legend(tempNames{:},'Location','Best');
saveas(gcf,[ imgSavePath '\fig1_day_' int2str(selectday) '.png']);
%% plot 2 - the averaged activity over days in tone-evoked period
corrFlat = [];allsFlat = [];
for i = 1:length(alls)
    allsFlat(:,i) = alls{i}(:);

    corrMat = corr(alls{i}(:,:)');
    upperTriFlag = triu(ones(size(corrMat)),1);
    corrFlat(:,i) = corrMat(logical(upperTriFlag));
    
end
tempNames = shortFilename;
corrFlat(sum(sum(isnan(corrFlat),2),3)>0,:,:) = [];
corrCorr = []; corrCorr = corr(corrFlat); 
totalCorr = corr(allsFlat);
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition',  [0.03, 0.4, 0.8, 0.6]);
subplot(1,2,1);imagesc(totalCorr); colorbar;
title('correlation of spike activity between different methods')
xticks(1:length(alls))
xticklabels(tempNames)
subplot(1,2,2);imagesc(mean(corrCorr,3)); colorbar;% caxis([0 1]);
title('correlation of correlation coefficients between different methods')
xticks(1:length(alls))
xticklabels(tempNames)
saveas(gcf,[ imgSavePath '\fig2_day_' int2str(selectday) '.png']);
close all;
%% plot 3 - save examples of every cells
neuronLim = 60;
if plotNeuron
    selectFrame = 1:2000; 
    for i = 1:neuronLim
        bsum = 0;
        for j = 1:length(filenameList); load([datapath '\' filenameList{j}(1:end-5) 'p.mat'],'p');
            if size(p,1) == 1; bsum = bsum + p{sessionIdx}{i}.b;else;bsum = bsum + p{sessionIdx,i}.b;end
        end
        figdays = figure('visible','off');
        set(figdays, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.04, 0.85, 0.96]);
        
        totalSubPlot = length(filenameList)+1; if strcmp(mouse,'cd036'); totalSubPlot = totalSubPlot + length(extraAlls); end
        subplot_tight(totalSubPlot,1,1);
        plot(selectDff(i,selectFrame)); hold on; 
        plot(ones(1,length(selectFrame))*bsum/length(filenameList),'Color',[0.8 0.8 0.8], 'LineWidth',2); title('dff')
        %tempTitles = {'foopsi','constrained foopsi', 'thresholded'};
        for j = 1:length(filenameList)
            subplot_tight(totalSubPlot,1,1+j); plot(alls{j}(i,selectFrame))
            if j~=length(filenameList); set(gca,'xtick',[]); end ; title(shortFilename{j})
        end
        if strcmp(mouse,'cd036') && ~isempty(extraAlls)
            for j = 1:length(extraAlls)
                subplot_tight(totalSubPlot,1,length(filenameList)+1+j);
                plot(extraAlls{j}(i,selectFrame));
                if j~=length(extraAlls); set(gca,'xtick',[]);end;  title(extraShortFilename{j})
            end
        end
        
        saveas(figdays,[ imgSavePath '\cell_' int2str(i) 'day_' int2str(selectday) '.png']);
        close (figdays);
    end
end
end

%% functions 
function [tPeakMean, fPeakMean] = getPeakAct(behavMatrix,act)
    tFrame = ceil(behavMatrix(behavMatrix(:,4)==1 | behavMatrix(:,4)==2,12)/2);
    fFrame = ceil(behavMatrix(behavMatrix(:,4)==3 | behavMatrix(:,4)==4,12)/2);
    selectFrame = -5 : 10;
    tAct = []; fAct = [];
    for k = 1:length(tFrame)
        tAct(:,:,k) = act(:,tFrame(k) + selectFrame);
        fAct(:,:,k) = act(:,fFrame(k) + selectFrame);
    end
    tActMean = mean(tAct,3); tPeakMean = max(smoothdata(tActMean,2,'gaussian',3),[],2);
    fActMean = mean(fAct,3); fPeakMean = max(smoothdata(fActMean,2,'gaussian',3),[],2);
    tPeakTrial = squeeze(max(smoothdata(tAct,2,'gaussian',3),[],2));
    fPeakTrial = squeeze(max(smoothdata(fAct,2,'gaussian',3),[],2));
end

function f = fillErrorbarPlot(xdata,yMean, ySEM,varargin)
    x = [xdata, fliplr(xdata)];
    y = [yMean+ySEM, fliplr(yMean-ySEM)];
    f = fill(x,y,varargin{:});
end