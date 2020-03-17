%% load data
clear; 
mouse = 'cd042';
switch mouse
    case 'cd017'
        preName = 'cd017_000_001';
        postName = 'cd017_016_007';
        target = 9514;
        foil = 16000;
    case 'cd036'
        preName = 'cd036_000_003';
        postName = 'cd036_017_001';
        target = 13454;
        foil = 22627;
    case 'cd037'
        preName = 'cd037_000_001';
        postName = 'cd037_017_003';
        target = 11314;
        foil = 6727;
    case 'cd042'
        preName = 'cd042_000_003';
        postName = 'cd042_016_001';
        target = 11314;
        foil = 6727;
    case 'cd044'
        preName = 'cd044_000_001';
        postName = 'cd044_016_001';
        target = 9514;
        foil = 16000;
        
end
preTuning = load(['D:\labData\excitatory\tuning\' preName '_Tuning\population\tuning.mat']);
postTuning = load(['D:\labData\excitatory\tuning\' postName '_Tuning\population\tuning.mat']);
savePath = 'D:\labData\excitatory\tuning\masterData\';
%ishere1 = load('T:\LabData6\ziyi\cd017\roi\ishere_plane0.mat'); ishere1 = ishere1.ishere{1};
%ishere2 = load('T:\LabData6\ziyi\cd017\roi\ishere_plane1.mat'); ishere2 = ishere2.ishere{2};
%ishere = [ishere1; ishere2];
%ishereFlag = sum(ishere(:,[1 17]),2)==2 ;

%% declare useful parameters and maps
%responsiveFlagPre= sum(preTuning.signrankToneH,1)>0;
%responsiveFlagPost= sum(postTuning.signrankToneH,1)>0;
pretoneFrames = 10;
rbmap = imresize(redbluecmap, [256, 3]);  % original color map contain just 11 colors, this increase it to 64
rbmap = min(max(rbmap, 0), 1);
nPlanes = 2;
nNeuron = size(preTuning.TCpretone_reorderCorr,4);

tones = [4000,4757,5657,6727,8000,...
    9514,11314,13454,16000,19027,22627,...
    26909,32000,38055,45255,53817,64000];
targIdx = find(tones==target);
foilIdx = find(tones==foil);

%% get useful data
responsiveFlagPre= preTuning.anovaPeakCorr(2,:) > 0 & sum(preTuning.anovaSignifToneCorr,1) > 0;
responsiveFlagPost= postTuning.anovaPeakCorr(2,:) > 0 & sum(postTuning.anovaSignifToneCorr,1) > 0;
respCellFlag = responsiveFlagPre | responsiveFlagPost;
img = preTuning.refImg;
rois = preTuning.roisBound;
neuronEachPlane = preTuning.neuronEachPlane;

preTuning.TCpretone_zscoreCorr = zscoreTC(preTuning.TCpretone_reorderCorr);
postTuning.TCpretone_zscoreCorr = zscoreTC(postTuning.TCpretone_reorderCorr);
[peakActPre,~] = getPeakAct(preTuning.TCpretone_zscoreCorr,pretoneFrames,preTuning.peakFrames);
[peakActPost,~] = getPeakAct(postTuning.TCpretone_zscoreCorr,pretoneFrames,postTuning.peakFrames);

%% plot1: 
colormapIndex = round(linspace(1,64,17));
plotMapIdx(img, rois, neuronEachPlane,jet,colormapIndex,...
   {preTuning.tuningPeak, postTuning.tuningPeak}, ...
    {responsiveFlagPre, responsiveFlagPost});
saveas(gcf,[savePath mouse '_figurePeak1.png']);
%% plot1.5:
tuningChange = preTuning.tuningPeak - postTuning.tuningPeak;
tempMax = prctile(tuningChange(responsiveFlagPre & responsiveFlagPost),95);
tempMin = prctile(tuningChange(responsiveFlagPre & responsiveFlagPost),5);
changeIdx = normalizeIndex(tuningChange, tempMax, tempMin);

plotMapIdx(img, rois, neuronEachPlane,rbmap,1:256,...
   changeIdx,responsiveFlagPre & responsiveFlagPost);
saveas(gcf,[savePath mouse '_figurePeak2.png']);
%% plot1.9: get different kinds of peak shifts

tPre = (targIdx == preTuning.tuningPeak);
fPre = (foilIdx == preTuning.tuningPeak);

if targIdx < foilIdx
    middlePre = (targIdx < preTuning.tuningPeak & foilIdx > preTuning.tuningPeak);
    tSidePre = (targIdx > preTuning.tuningPeak & foilIdx > preTuning.tuningPeak);
    fSidePre = (targIdx < preTuning.tuningPeak & foilIdx < preTuning.tuningPeak);
else
    middlePre = (targIdx > preTuning.tuningPeak & foilIdx < preTuning.tuningPeak);
    tSidePre = (targIdx < preTuning.tuningPeak & foilIdx < preTuning.tuningPeak);
    fSidePre = (targIdx > preTuning.tuningPeak & foilIdx > preTuning.tuningPeak);
end

prepostTuning.peak.pre = preTuning.tuningPeak;
prepostTuning.peak.post = postTuning.tuningPeak;
prepostTuning.peak.tuningChange = tuningChange;
prepostTuning.peak.tuningUp = sum(tuningChange(respCellFlag)>0);
prepostTuning.peak.tuningDown = sum(tuningChange(respCellFlag)<0);
prepostTuning.peak.tuningStay = sum(tuningChange(respCellFlag)==0);

prepostTuning.peak.dist.tPreDist = abs(targIdx - preTuning.tuningPeak);
prepostTuning.peak.dist.fPreDist = abs(foilIdx - preTuning.tuningPeak);
prepostTuning.peak.dist.tPostDist = abs(targIdx - postTuning.tuningPeak);
prepostTuning.peak.dist.fPostDist = abs(foilIdx - postTuning.tuningPeak);

if targIdx < foilIdx
    prepostTuning.peak.tPre.stay = sum(tuningChange(tPre & respCellFlag)==0);
    prepostTuning.peak.tPre.switch = sum(tuningChange(tPre & respCellFlag)>0);
    prepostTuning.peak.tPre.away = sum(tuningChange(tPre & respCellFlag)<0);

    prepostTuning.peak.fPre.stay = sum(tuningChange(fPre & respCellFlag)==0);
    prepostTuning.peak.fPre.switch = sum(tuningChange(fPre & respCellFlag)<0);
    prepostTuning.peak.fPre.away = sum(tuningChange(fPre & respCellFlag)>0);

    prepostTuning.peak.middlePre.stay = sum(tuningChange(middlePre & respCellFlag)==0);
    prepostTuning.peak.middlePre.T = sum(tuningChange(middlePre & respCellFlag)<0);
    prepostTuning.peak.middlePre.F = sum(tuningChange(middlePre & respCellFlag)>0);

    prepostTuning.peak.tSidePre.stay = sum(tuningChange(tSidePre & respCellFlag)<=0);
    prepostTuning.peak.tSidePre.closeT = sum(tuningChange(tSidePre & respCellFlag)>0 & (postTuning.tuningPeak(tSidePre & respCellFlag) < mean([targIdx foilIdx]) )  );
    prepostTuning.peak.tSidePre.closeF = sum(tuningChange(tSidePre & respCellFlag)>0 & (postTuning.tuningPeak(tSidePre & respCellFlag) > mean([targIdx foilIdx]) )  );

    prepostTuning.peak.fSidePre.stay = sum(tuningChange(fSidePre & respCellFlag)>=0);
    prepostTuning.peak.fSidePre.closeT = sum(tuningChange(fSidePre & respCellFlag)<0 & (postTuning.tuningPeak(fSidePre & respCellFlag) < mean([targIdx foilIdx]) )  );
    prepostTuning.peak.fSidePre.closeF = sum(tuningChange(fSidePre & respCellFlag)<0 & (postTuning.tuningPeak(fSidePre & respCellFlag) > mean([targIdx foilIdx]) )  );

else
    prepostTuning.peak.tPre.stay = sum(tuningChange(tPre & respCellFlag)==0);
    prepostTuning.peak.tPre.switch = sum(tuningChange(tPre & respCellFlag)<0);
    prepostTuning.peak.tPre.away = sum(tuningChange(tPre & respCellFlag)>0);

    prepostTuning.peak.fPre.stay = sum(tuningChange(fPre & respCellFlag)==0);
    prepostTuning.peak.fPre.switch = sum(tuningChange(fPre & respCellFlag)>0);
    prepostTuning.peak.fPre.away = sum(tuningChange(fPre & respCellFlag)<0);

    prepostTuning.peak.middlePre.stay = sum(tuningChange(middlePre & respCellFlag)==0);
    prepostTuning.peak.middlePre.T = sum(tuningChange(middlePre & respCellFlag)>0);
    prepostTuning.peak.middlePre.F = sum(tuningChange(middlePre & respCellFlag)<0);

    prepostTuning.peak.tSidePre.stay = sum(tuningChange(tSidePre & respCellFlag)>=0);
    prepostTuning.peak.tSidePre.closeT = sum(tuningChange(tSidePre & respCellFlag)<0 & (postTuning.tuningPeak(tSidePre & respCellFlag) > mean([targIdx foilIdx]) )  );
    prepostTuning.peak.tSidePre.closeF = sum(tuningChange(tSidePre & respCellFlag)<0 & (postTuning.tuningPeak(tSidePre & respCellFlag) < mean([targIdx foilIdx]) )  );

    prepostTuning.peak.fSidePre.stay = sum(tuningChange(fSidePre & respCellFlag)<=0);
    prepostTuning.peak.fSidePre.closeT = sum(tuningChange(fSidePre & respCellFlag)>0 & (postTuning.tuningPeak(fSidePre & respCellFlag) > mean([targIdx foilIdx]) )  );
    prepostTuning.peak.fSidePre.closeF = sum(tuningChange(fSidePre & respCellFlag)>0 & (postTuning.tuningPeak(fSidePre & respCellFlag) < mean([targIdx foilIdx]) )  );

end
prepostTuning.peak.nResp = sum(respCellFlag);
%% plot2: change in T or F activity
targDiff = squeeze((peakActPre(targIdx,:,:)- peakActPost(targIdx,:,:)));
foilDiff = squeeze((peakActPre(foilIdx,:,:)- peakActPost(foilIdx,:,:)));
tempMax = prctile([mean(targDiff) mean(foilDiff)],90);
tempMin = prctile([mean(targDiff) mean(foilDiff)],10);

targDiffIdx = normalizeIndex(mean(targDiff,1), tempMax, tempMin);
foilDiffIdx = normalizeIndex(mean(foilDiff,1), tempMax, tempMin);

plotMapIdx(img, rois, neuronEachPlane,rbmap,1:256,...
   {targDiffIdx, foilDiffIdx}, ...
    {responsiveFlagPre & responsiveFlagPost, responsiveFlagPre & responsiveFlagPost});
saveas(gcf,[savePath mouse '_figureActChange1.png']);
%% plot3: T/F selectivity before and after learning
tempT = squeeze(mean(peakActPre(targIdx,:,:),2));
tempF = squeeze(mean(peakActPre(foilIdx,:,:),2));
tfPre = (abs(tempT) - abs(tempF)) ./ (abs(tempT) + abs(tempF));

tempT = squeeze(mean(peakActPost(targIdx,:,:),2));
tempF = squeeze(mean(peakActPost(foilIdx,:,:),2));
tfPost = (abs(tempT) - abs(tempF)) ./ (abs(tempT) + abs(tempF));

tempMax = prctile([tfPre;tfPost],90);
tempMin = prctile([tfPre;tfPost],10);

targDiffIdx = normalizeIndex(tfPre, tempMax, tempMin);
foilDiffIdx = normalizeIndex(tfPost, tempMax, tempMin);

plotMapIdx(img, rois, neuronEachPlane,rbmap,1:256,...
   {targDiffIdx, foilDiffIdx}, ...
    {responsiveFlagPre & responsiveFlagPost, responsiveFlagPre & responsiveFlagPost});
saveas(gcf,[savePath mouse '_figureSelect1.png']);
%% plot4: change in T/F selectivity
tfDiff = tfPre - tfPost;
prepostTuning.tfSelect.tfDiff = tfDiff;
prepostTuning.tfSelect.tfPre = tfPre;
prepostTuning.tfSelect.tfPost = tfPost;

tempMax = prctile(tfDiff,90);
tempMin = prctile(tfDiff,10);

tfDiffIdx = normalizeIndex(tfDiff, tempMax, tempMin);

plotMapIdx(img, rois, neuronEachPlane,rbmap,1:256,tfDiffIdx,...
    responsiveFlagPre & responsiveFlagPost);
saveas(gcf,[savePath mouse '_figureSelect2.png']);
%% plot4: visualize decoder weight 

tempT = squeeze(peakActPre(targIdx,:,respCellFlag));
tempF = squeeze(peakActPre(foilIdx,:,respCellFlag));

decoderPre = runDecoder(tempT, tempF,0.6);
[h,~] = ttest(decoderPre.cellAcc(:,:,2)',decoderPre.cellAccShuf(:,:,2)','tail','right');
signifFlagPre = logical(h);
meanAccPre = mean(decoderPre.cellAcc(signifFlagPre,:,2),2);
cellWPre = mean(decoderPre.cellW(signifFlagPre,:),2);
[~,tempidx] = sort(abs(mean(cellWPre,2)),'descend');
%figure; plot(abs(mean(cellW(tempidx,:),2))); hold on; plot(std(cellW(tempidx,:),0,2));

tempT = squeeze(peakActPost(targIdx,:,respCellFlag));
tempF = squeeze(peakActPost(foilIdx,:,respCellFlag));
decoderPost = runDecoder(tempT, tempF,0.6);
[h,~] = ttest(decoderPost.cellAcc(:,:,2)',decoderPost.cellAccShuf(:,:,2)','tail','right');
signifFlagPost = logical(h);
meanAccPost = mean(decoderPost.cellAcc(signifFlagPost,:,2),2);
cellWPost = mean(decoderPost.cellW(signifFlagPost,:),2);

prepostTuning.decode.cellWPre = cellWPre;
prepostTuning.decode.cellWPost = cellWPost;

tempMax = prctile([cellWPre;cellWPost],90); tempMin = prctile([cellWPre;cellWPost],10);

decoderIdxPre = normalizeIndex(cellWPre, tempMax, tempMin);
tempidx = find(respCellFlag); tempidx = tempidx(signifFlagPre);
preFlag = zeros(1,nNeuron);preFlag(tempidx) = 1;
preIdx = zeros(1,nNeuron);preIdx(tempidx) = decoderIdxPre;

decoderIdxPost = normalizeIndex(cellWPost, tempMax, tempMin);
tempidx = find(respCellFlag); tempidx = tempidx(signifFlagPost);
postFlag = zeros(1,nNeuron);postFlag(tempidx) = 1;
postIdx = zeros(1,nNeuron);postIdx(tempidx) = decoderIdxPost;

plotMapIdx(img, rois, neuronEachPlane,rbmap,1:256,{preIdx postIdx},...
    {preFlag, postFlag});
saveas(gcf,[savePath mouse '_figureDecode1.png']);
%% plot: change in decoder weight
% here we take all the neurons that are selective either in pre or post
decodeFlag = (signifFlagPre | signifFlagPost);
cellWPre = mean(decoderPre.cellW(decodeFlag,:),2);
cellWPost = mean(decoderPost.cellW(decodeFlag,:),2);

prepostTuning.decode.allCellWPre = cellWPre;
prepostTuning.decode.allCellWPost = cellWPost;

cellPrePost = cellWPre - cellWPost;
tempMax = prctile([cellPrePost],90); tempMin = prctile([cellPrePost],10);
idxPrePost = normalizeIndex(cellPrePost, tempMax, tempMin);

tempidx = find(respCellFlag); tempidx = tempidx(decodeFlag);
prePostFlag = zeros(1,nNeuron);prePostFlag(tempidx) = 1;
prePostIdx = zeros(1,nNeuron);prePostIdx(tempidx) = idxPrePost;

plotMapIdx(img, rois, neuronEachPlane,rbmap,1:256,prePostIdx,...
    prePostFlag);
saveas(gcf,[savePath mouse '_figureDecode2.png']);
%% plot: change in decoder weight based on their initial decoding weight
% restrict to discriminating neurons
%tempT = squeeze(peakActPre(targIdx,:,:));
%tempF = squeeze(peakActPre(foilIdx,:,:));
%[tPreferIdx,~] = ttest(tempT, tempF, 'tail','left','alpha',0.025);
%[fPreferIdx,~] = ttest(tempT, tempF, 'tail','right','alpha',0.025);
tFlagPre = mean(decoderPre.cellW(:,:),2) > 0 & signifFlagPre'; tFlagPre = tFlagPre(decodeFlag);
fFlagPre = mean(decoderPre.cellW(:,:),2) < 0 & signifFlagPre'; fFlagPre = fFlagPre(decodeFlag);
tFlagPost = mean(decoderPost.cellW(:,:),2) > 0 & signifFlagPost'; tFlagPost = tFlagPost(decodeFlag);
fFlagPost = mean(decoderPost.cellW(:,:),2) < 0 & signifFlagPost'; fFlagPost = fFlagPost(decodeFlag);

prepostTuning.decode.prePost.tFlagPre = tFlagPre;
prepostTuning.decode.prePost.fFlagPre = fFlagPre;
prepostTuning.decode.prePost.tFlagPost = tFlagPost;
prepostTuning.decode.prePost.fFlagPost = fFlagPost;

rep = 100;
shufData = [];
shufDataPre = [];
for i = 1:100
    randidx1 = randperm(sum(respCellFlag)-sum(decodeFlag));
    temp = mean(decoderPre.cellW(~decodeFlag,:),2); tempPre = temp(randidx1);
    randidx2 = randperm(sum(respCellFlag)-sum(decodeFlag));
    temp = mean(decoderPost.cellW(~decodeFlag,:),2); tempPost = temp(randidx2);
    shufData = [shufData;tempPre-tempPost];
    shufDataPre = [shufDataPre;tempPre];
end
downFlag = cellPrePost > std(shufData)*2;
upFlag = cellPrePost < -std(shufData)*2;

tChange.up = sum(tFlagPre&upFlag);
tChange.down = sum(tFlagPre&downFlag);
tChange.stay = sum(tFlagPre&(~upFlag & ~downFlag));
tChange.downSwitch = sum(tFlagPre&downFlag&fFlagPost);
prepostTuning.decode.tChange = tChange;

fChange.up = sum(fFlagPre&upFlag);
fChange.down = sum(fFlagPre&downFlag);
fChange.stay = sum(fFlagPre&(~upFlag & ~downFlag));
fChange.upSwitch = sum(fFlagPre&upFlag&tFlagPost);
prepostTuning.decode.fChange = fChange;

nChange.up = sum((~tFlagPre & ~fFlagPre)&upFlag);
nChange.down =  sum((~tFlagPre & ~fFlagPre)&downFlag);
nChange.stay = sum((~tFlagPre & ~fFlagPre)&(~upFlag & ~downFlag));

nChange.upT =  sum((~tFlagPre & ~fFlagPre)&upFlag&tFlagPost);
nChange.upF = sum((~tFlagPre & ~fFlagPre)&upFlag&fFlagPost);
nChange.downT = sum((~tFlagPre & ~fFlagPre)&downFlag&tFlagPost);
nChange.downF = sum((~tFlagPre & ~fFlagPre)&downFlag&fFlagPost);
prepostTuning.decode.nChange = nChange;
prepostTuning.decode.nDecode = sum(decodeFlag);

save([savePath mouse '_prePostTuning.mat'],'prepostTuning');

%% all the functions
function [f1] = plotMapIdx(img, rois, neuronEachPlane,C,colormapIndex,plotIndex, plotFlag)
    f1 = figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.04, 0.6, 0.9]);
    if ~iscell(plotIndex)
        plotIndex = {plotIndex};
    end
    if ~iscell(plotFlag)
        plotFlag = {plotFlag};
    end
    nPlanes = length(neuronEachPlane);
    for i = 1:nPlanes
        neuronPlane = [0 cumsum(neuronEachPlane)];
        for k = 1:length(plotIndex)
            subplot(length(plotIndex),2,(k-1)*2+i);             
            imagesc(img{i});colormap gray;hold on;
            ylim([0 size(img{i},1)]);xlim([0 size(img{1},2)]);
            for j = 1:neuronEachPlane(i)
                cellIndex = neuronPlane(i) + j;
                x = rois{i}{j}(:,1); %freehand rois have the outlines in x-y coordinates
                y = rois{i}{j}(:,2); %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
                if plotFlag{k}(cellIndex)
                    %plot(x,y,'.','color',C(colormapIndex(peakIndex(cellIndex)),:),'MarkerSize',1);
                    patch(x,y,C(colormapIndex(plotIndex{k}(cellIndex)),:),'EdgeColor','none');
                else
                    patch(x,y,[0.8 0.8 0.8],'EdgeColor','none');
                    %plot(x,y,'.','color',[0.8 0.8 0.8],'MarkerSize',1);
                end
            end
            title(['Plane' int2str(i)])
            axis off;
        end
    end
end

function [peakAct,baseAct] = getPeakAct(TC,pretoneFrames,peakFrames)
nNeuron = size(TC,4);
peakAct = zeros(size(TC,2),size(TC,3),size(TC,4));
baseAct = zeros(size(TC,2),size(TC,3),size(TC,4));
for i = 1:nNeuron
    peakAct(:,:,i) = squeeze(TC(pretoneFrames+peakFrames(i),:,:,i));
    baseAct(:,:,1) = squeeze(TC(pretoneFrames,:,:,i));
end
end

function pretoneTC = zscoreTC(pretoneTC)
    nNeuron = size(pretoneTC,4);
    temp = reshape(pretoneTC,8500,nNeuron);
    temp = zscore(temp, 0,1);
    pretoneTC = reshape(temp,50,17,10,nNeuron);
end

function index1 = normalizeIndex(index1, tempMax, tempMin)

tempLim = max([abs(tempMax) abs(tempMin)]);
index1 = index1/tempLim;
index1 = round(index1*127.5 + 128.5);
index1(index1 > 256) = 256;
index1(index1 < 1) = 1;

end

function decoder = runDecoder(data1,data2, percTrain)
    nNeuron = size(data1,2);
    rep = 100;
    nTrial = size(data1,1);
    nTrain = round (percTrain * nTrial);
    nTest = nTrial - nTrain;
    decoder.cellW = zeros(nNeuron,rep);
    decoder.cellAcc = zeros(nNeuron,rep,2);
    decoder.popW = zeros(nNeuron,rep);
    decoder.popAcc = zeros(rep,2);

    decoder.cellWShuf = zeros(nNeuron,rep);
    decoder.cellAccShuf = zeros(nNeuron,rep,2);
    decoder.popWShuf = zeros(nNeuron,rep);
    decoder.popAccShuf = zeros(rep,2);
    for i = 1:rep
        idxT = randperm(nTrial);
        idxF = randperm(nTrial);
        trainT = data1(idxT(1:nTrain),:);
        trainF = data2(idxF(1:nTrain),:);
        testT = data1(idxT(nTrain+1:end),:);
        testF = data2(idxF(nTrain+1:end),:);

        
        [tempCellW, tempCellAcc, tempPopW, tempPopAcc] = getDecoderW(trainT,trainF,testT,testF);
        decoder.cellW(:,i) = tempCellW; decoder.cellAcc(:,i,:) = tempCellAcc;
        decoder.popW(:,i) = tempPopW; decoder.popAcc(i,:) = tempPopAcc;

        allData = [data1;data2];
        idxR = randperm(nTrial*2);
        trainT = allData(idxR(1:nTrain),:);
        trainF = allData(idxR(nTrain+1:nTrain*2),:);
        testT = allData(idxR(nTrain*2+1:nTrain*2+nTest),:);
        testF = allData(idxR(nTrain*2+nTest+1:end),:);
        [tempCellW, tempCellAcc, tempPopW, tempPopAcc] = getDecoderW(trainT,trainF,testT,testF);
        decoder.cellWShuf(:,i) = tempCellW; decoder.cellAccShuf(:,i,:) = tempCellAcc;
        decoder.popWShuf(:,i) = tempPopW; decoder.popAccShuf(i,:) = tempPopAcc;
    end
end


function [cellW, cellAcc, popW, popAcc] = getDecoderW(trainT,trainF,testT,testF)
    nTrain = size(trainT,1);
    nTest = size(testT,1);

    nNeuron = size(trainT,2);
    cellW = zeros(nNeuron,1);
    cellAcc = zeros(nNeuron,2);
    popAcc = zeros(1,2);
    for j = 1:nNeuron    
        tempW = (mean(trainT(:,j)) - mean(trainF(:,j))) / (var(trainT(:,j))/2 + var(trainF(:,j))/2);
        cellW(j) = tempW;
        predictT = tempW * trainT(:,j); predictF = tempW * trainF(:,j); 
        c = mean(predictT)/2 + mean(predictF)/2;
        correctT = predictT>c;correctF = predictF<c;
        cellAcc(j,1) = (sum(correctT) + sum(correctF)) / (2*nTrain);
        predictT = tempW * testT(:,j); correctT = predictT>c;
        predictF = tempW * testF(:,j); correctF = predictF<c;
        cellAcc(j,2) = (sum(correctT) + sum(correctF)) / (2*(nTest));
    end
    warning ('off')
    tempW = (mean(trainT) - mean(trainF)) * inv(cov(trainT)/2 + cov(trainF)/2);
    popW = tempW;
    predictT = trainT * tempW'; predictF = trainF * tempW'; 
    c = mean(predictT)/2 + mean(predictF)/2;
    correctT = predictT>c;correctF = predictF<c;
    popAcc(1) = (sum(correctT) + sum(correctF)) / (2*nTrain);
    predictT = testT * tempW'; predictF = testT * tempW'; 
    correctT = predictT>c;correctF = predictF<c;
    popAcc(2) = (sum(correctT) + sum(correctF)) / (2*(nTest));

end