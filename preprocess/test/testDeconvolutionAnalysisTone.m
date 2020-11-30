clear;
datapath ='C:\Users\zzhu34\Documents\tempdata\deconv_test_cd036\smooth3000\behavior';
behavpath = 'C:\Users\zzhu34\Documents\tempdata\excitData\cd036\behavior';
%filenameList = {'cd036_confoo_smin.mat','cd036_foo_smin.mat','cd036_thre_smin.mat'};
filenameList = {'cd036_calman_thre90_optb_nosmin.mat','cd036_calman_thre95_optb_nosmin.mat','cd036_calman_thre99_optb_nosmin.mat'};
%filenameList = {'cd036_foo_lambda0.mat','cd036_foo_lambda05.mat','cd036_foo_lambda10.mat'};
shortFilename = [];
for i = 1:length(filenameList); shortFilename{i} = filenameList{i}(7:end-4); shortFilename{i}(strfind(shortFilename{i},'_'))=' ';end
for i = 1:length(filenameList)
    disp([filenameList{i} ' started!'])
    load([datapath '\' filenameList{i}],'s','c','day');
    for j = 1:length(day)
        behavMatrix = importdata([behavpath '\cd036_' int2str(day(j)-1) 'v1.txt']);
        [tempTPM, tempFPM] = getPeakAct(behavMatrix,s{j});
        tPeakMean(:,j,i) = tempTPM;fPeakMean(:,j,i) = tempFPM;
    end
end
load([datapath '\selectDff.mat'],'selectDff');
disp('selected dff started!')
for j = 1:length(day)
    behavMatrix = importdata([behavpath '\cd036_' int2str(day(j)-1) 'v1.txt']);
    [tempTP, tempFP] = getPeakAct(behavMatrix,selectDff{j});
    tPeakMeanDff(:,j) = tempTP;fPeakMeanDff(:,j) = tempFP;
end
%% plot 1 - the averaged activity over days in tone-evoked period
figure; hold on;
dffFlat = [tPeakMeanDff(:)' fPeakMeanDff(:)'];
for i = 1:3; spikeFlat(:,i) = [reshape(tPeakMean(:,:,i),1,[]) reshape(fPeakMean(:,:,i),1,[])];end
[N,edges,bin] = histcounts(dffFlat,linspace(0,prctile(dffFlat,98),21));
for j = 1:20 
    tempX(j) = mean(edges(j:j+1));
    for i = 1:3
        tempMean(j,i) = mean(spikeFlat(bin==j,i)); tempSEM(j,i) = std(spikeFlat(bin==j,i))/sqrt(N(j));     
    end
end
for i = 1:3; plot(tempX,tempMean(:,i),'Color',matlabColors(i)); end
tempLegend = [];
for i = 1:3
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
for j = 1:3; subplot_tight(1,3,j,[0.15,0.06]);hold on; title(shortFilename{j})
for i= 1:7   
    scatter([tPeakMeanDff(:,i)' fPeakMeanDff(:,i)'],...
        [reshape(tPeakMean(:,i,j),1,[]) reshape(fPeakMean(:,i,j),1,[])],...
        5, 'filled' , 'o','MarkerFaceAlpha', 0.4);
end
ylabel('spike'); xlabel('dff'); 
xlim([0 prctile([tPeakMeanDff(:)' fPeakMeanDff(:)'],99)]);
ylim([0 prctile([reshape(tPeakMean(:,:,j),1,[]) reshape(fPeakMean(:,:,j),1,[])],99)]);
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