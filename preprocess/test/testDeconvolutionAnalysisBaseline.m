clear;
datapath ='C:\Users\zzhu34\Documents\tempdata\deconv_test_cd036\smooth3000\baseline';
%filenameList = {'cd036_confoo_smin.mat','cd036_foo_smin.mat','cd036_thre_smin.mat'};
%filenameList = {'cd036_calman_confoo_optb.mat','cd036_calman_foo90_optb_nosmin.mat',...
%    'cd036_calman_foo95_optb_nosmin.mat','cd036_calman_foo95_optb.mat','cd036_calman_foo99_optb.mat','cd036_calman_thre95_optb.mat'};
filenameList = {'cd036_calman_ar2_confoo_optb.mat','cd036_calman_ar2_foo90_optb_nosmin.mat','cd036_calman_foo90_optb_nosmin.mat',...
    'cd036_calman_ar2_foo95_optb_nosmin.mat','cd036_calman_foo95_optb_nosmin.mat'};
shortFilename = []; for i = 1:length(filenameList)
    shortFilename{i} = filenameList{i}(14:end-4); shortFilename{i}(strfind(shortFilename{i},'_'))=' ';end
for i = 1:length(filenameList)
    disp([filenameList{i} ' started!'])
    load([datapath '\' filenameList{i}],'s','c','day');
    alls{i} = reshape(cell2mat(s),size(s{1},1),size(s{1},2),length(day));
    for j = 1:length(day)
        spikeMean(:,j,i) = mean(s{j},2);
        cMean(:,j,i) = mean(c{j},2);
    end
end

plotNeuron = true;

load([datapath '\selectDff.mat'],'selectDff');
disp('selected dff started!')
for j = 1:length(day)
    dffMean(:,j) = mean(selectDff{j},2);
end
%% plot 1 - the averaged activity over days in tone-evoked period
figure; hold on;
dffFlat = dffMean(:)';tempX =[];
for i = 1:length(filenameList); spikeFlat(:,i) = reshape(spikeMean(:,:,i),1,[]) ; end
[N,edges,bin] = histcounts(dffFlat,linspace(0,prctile(dffFlat,98),21));
for j = 1:20 
    tempX(j) = mean(edges(j:j+1));
    for i = 1:length(filenameList)
        tempMean(j,i) = mean(spikeFlat(bin==j,i)); tempSEM(j,i) = std(spikeFlat(bin==j,i))/sqrt(N(j));     
    end
end
for i = 1:length(filenameList); plot(tempX,tempMean(:,i),'Color',matlabColors(i)); end
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
%% plot 1.1 - the averaged activity over days in tone-evoked period, c and s
figure; hold on;tempX =[];
for i = 1:length(filenameList)
    spikeFlat(:,i) = reshape(spikeMean(:,:,i),1,[]) ; 
    cFlat(:,i) = reshape(cMean(:,:,i),1,[]) ;
end
for i = 1:length(filenameList)
    [N,edges,bin] = histcounts(cFlat(:,i),linspace(0,prctile(cFlat(:,i),98),21));
    for j = 1:20 
        tempX(j,i) = mean(edges(j:j+1));
        tempMean(j,i) = mean(spikeFlat(bin==j,i)); tempSEM(j,i) = std(spikeFlat(bin==j,i))/sqrt(N(j));      
    end
end
for i = 1:length(filenameList); plot(tempX(:,i),tempMean(:,i),'Color',matlabColors(i)); end
for i = 1:length(filenameList)
    f = fillErrorbarPlot(tempX(:,i)',tempMean(:,i)', tempSEM(:,i)',matlabColors(i),'LineStyle','none');
    f.FaceAlpha = 0.1;
    scatter(cFlat(:,i)',spikeFlat(:,i),5,matlabColors(i), 'filled' , 'o', 'MarkerFaceAlpha', 0.1); 
end
ylabel('spike'); xlabel('dff'); 
xlim([0 prctile(cFlat(:),98)]);
ylim([0 prctile(spikeFlat(:),99)]);
legend(shortFilename{:},'Location','Best');
title('Tone-evoked activity')
%% plot 1.2 - the averaged activity over days in tone-evoked period
figure; 
for j = 1:length(filenameList); subplot_tight(1,length(filenameList),j,[0.15,0.06]);hold on; title(shortFilename{j})
for i= 1:7   
    scatter(dffMean(:,i)',reshape(spikeMean(:,i,j),1,[]),...
        5, 'filled' , 'o','MarkerFaceAlpha', 0.4);
end
ylabel('spike'); xlabel('dff'); 
xlim([0 prctile(dffMean(:,i)',99)]);
ylim([0 prctile(reshape(spikeMean(:,i,j),1,[]),99)]);
end

%% plot 2 - the averaged activity over days in tone-evoked period
figure; 
spikeCount = []; ISIcount = [];
for i = 1:length(alls)
    for j = 1:size(alls{i},1)
        for k = 1:length(day)
            spikes = find(alls{i}(j,:,k)>1e-5) ;
            ISI = diff(spikes);
            [N,edges,bins] = histcounts(ISI,1:100,'Normalization','probability');
            ISIcount{i}(:,j,k) = N;
            spikeCount{i}(j,k) = length(spikes);
        end
    end
    spikeCountAvg = mean(spikeCount{i},2);ISIcountAvg = nanmean(nanmean(ISIcount{i},3),2);
    subplot(1,2,1); cdfplot(spikeCountAvg); hold on;
    subplot(1,2,2); plot(ISIcountAvg); hold on;
end
subplot(1,2,1); title('spike count distribution'); legend(shortFilename{:},'Location','Best');
subplot(1,2,2); title('spike ISI distribution'); legend(shortFilename{:},'Location','Best');
%% plot 3 - the averaged activity over days in tone-evoked period
corrFlat = [];allsFlat = [];
for i = 1:length(alls)
    allsFlat(:,i) = alls{i}(:);
    for j = 1:length(day)
        corrMat = corr(alls{i}(:,:,j)');
        upperTriFlag = triu(ones(size(corrMat)),1);
        corrFlat(:,i,j) = corrMat(logical(upperTriFlag));
    end
end
corrFlat(sum(sum(isnan(corrFlat),2),3)>0,:,:) = [];
corrCorr = []; for j = 1:length(day); corrCorr(:,:,j) = corr(corrFlat(:,:,j)); end
totalCorr = corr(allsFlat);
figure; subplot(1,2,1);imagesc(totalCorr); colorbar;
title('correlation of spike activity between different methods')
xticks(1:length(alls))
xticklabels(shortFilename)
subplot(1,2,2);imagesc(mean(corrCorr,3)); colorbar;% caxis([0 1]);
title('correlation of correlation coefficients between different methods')
xticks(1:length(alls))
xticklabels(shortFilename)
%% plot 4 - save examples of every cells

if plotNeuron
    selectFrame = 1:2000; selectDay=1;
    for i = 1:size(alls{1},1)
        bsum = 0;
        for j = 1:length(filenameList); load([datapath '\' filenameList{1}],'p');
            bsum = bsum + p{selectDay,i}.b; end
        figdays = figure('visible','off');
        set(figdays, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.04, 0.85, 0.96]);
        subplot_tight(length(filenameList)+1,1,1);
        plot(selectDff{selectDay}(i,selectFrame)); hold on; 
        plot(ones(1,length(selectFrame))*bsum/length(filenameList),'Color',[0.8 0.8 0.8], 'LineWidth',2); title('dff')
        %tempTitles = {'foopsi','constrained foopsi', 'thresholded'};
        for j = 1:length(filenameList)
            subplot_tight(length(filenameList)+1,1,1+j); plot(alls{j}(i,selectFrame,selectDay))
            if j~=length(filenameList); set(gca,'xtick',[]); end ; title(shortFilename{j})
        end
        saveas(figdays,[ 'cell' int2str(i) '.png']);
        close (figdays);
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