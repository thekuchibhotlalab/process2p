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
global nPlanes

nPlanes = str2double(p.Results.nPlanes);
filename = p.Results.filename;
fn1 = p.Results.fn1; fn2 = p.Results.fn2;

if iscell(filename) && length(filename) >1
    disp('ERROR - More than one session selected!')
end

%---------CHECK NUMBER OF CHANNELS-----------
%data = load([p.Results.sbxpath sep filename '.mat']); 
%infosbx = data.info;
%if isempty(infosbx.otparam)
%    check_nPlanes = 1;
%else
%    check_nPlanes = infosbx.otparam(3);
%end
%if check_nPlanes ~= nPlanes
%    disp('ERROR - nPlanes in the parameter not consistent with sbx file.')
%    pause;
%end

%---------GET RELEVANT PARAMETERS-----------
[nFuncChannel, functionalChannel, roiType] = func_getFuncChanRoiType(varargin{:});

%---------GET CALCIUM TRACES-----------
[TC, neuronEachPlane, roisCoord] = func_loadTCRoi(varargin{:});

%---------CHECK IF NUMBER OF CHANNELS IS CORRECT-----------
if size(TC,1) ~= nFuncChannel
    disp('ERROR - Number of functional channels not correct')
end
%---------GET TUNING DATA FOR EACH CHANNEL-----------
if nFuncChannel == 1
    %CREATE SAVE FOLDER
    savePath = [p.Results.savepath sep p.Results.filename '_Tuning'];
    if exist(savePath,'dir') ~= 7
        mkdir(savePath);
        mkdir([savePath '/singleNeuron']);
        mkdir([savePath '/population']);
    end
    %PROCESS TUNING DATA
    TC_allPlane = cat(1,TC{:})';
    getTuning_oneChannel(TC_allPlane,neuronEachPlane{1},roisCoord(1,:),savePath,p.Results.suite2ppath,fn1,fn2);

elseif nFuncChannel == 2
    for i = 1:nFuncChannel        
        %CREATE SAVE FOLDER
        savePath = [p.Results.savepath sep p.Results.filename '_Tuning' sep 'chan' int2str(i)];
        if exist(savePath,'dir') ~= 7
            mkdir(savePath);
            mkdir([savePath '/singleNeuron']);
            mkdir([savePath '/population']);
        end
        %PROCESS TUNING DATA
        TC_channel = cat(1,TC{i,:})';
        getTuning_oneChannel(TC_channel,neuronEachPlane{i},roisCoord(i,:),savePath,p.Results.suite2ppath,fn1,fn2);
    end
end
%---------END OF GET TUNING FUNCTION-----------
end

%---------FUNCTION TO PROCESS CALCIUM TRACES-----------
function getTuning_oneChannel(TCraw,neuronEachPlane,roisBound,savePath,suite2ppath,fn1,fn2)

%---------DECLARE SOME PARAMETERS-----------
global nPlanes
sep = '\';
sumTC = sum(TCraw(1:end-1,:),1); 
TCraw(:,sumTC==0) = 5000 + randi(100,size(TCraw,1),sum(sumTC==0));
nNeuron = size(TCraw,2);
%---------DECLARE PROTOCOL PARAMETERS USING CUSTOM INPUT-----------
feval(fn1);
%---------SHIFT THE TC FOR PRETONE PERIOD-----------
TCraw = TCraw(1:nFramesPerTone*nTones*nTrials,:);
TC = circshift(TCraw,1-toneOnset,1); % shift so that tone is played ON FRAME 1
%---------GAUSSIAN FILTER TO SMOOTH THE TRACES-----------
if exist('smoothArg'); TC = smoothdata(TC,1,smoothArg{:}); end
%---------COMPUTE DFF-----------
if exist('dffArg'); TC = fn_getDff(TC,dffArg{:}); end
%---------RESHAPE TC AND PRETONE-----------
TCpretone = circshift(TC,pretoneFrames,1); % the tone is played ON FRAME 11 AFTER THIS SHIFT
TCpretone = reshape(TCpretone,nFramesPerTone,nTones,nTrials,nNeuron);
%---------REORDER TC TO ALIGN WITH TONE FREQUENCY ORDER-----------
TCpretone_reorder = zeros(size(TCpretone));
for x=1:nTones; TCpretone_reorder(:,x,:,:)=TCpretone(:,toneindex(x),:,:); end
%---------COMPUTE MEAN AND MEDIAN TC OF ALL TRIALS-----------
trialMean = squeeze(nanmean(TCpretone_reorder,3));
trialMedian = squeeze(nanmedian(TCpretone_reorder,3));
toneMean = squeeze(nanmean(nanmean(TCpretone_reorder,2),3));
%---------PLOT ALL TONE EVOKED ACTIVITY-----------
[peakValue,peakFrames] = max(toneMean(pretoneFrames+1:pretoneFrames+peakFrameBin,:),[],1);
fn_getTuningPSTH(trialMean,toneMean,pretoneFrames,peakFrames,peakValue,toneLabel,toneindex,savePath); %#ok<*USENS>
%---------COMPUTE TUNING CURVE OF ALL NEURONS-----------
[maxValue,peakFrames] = max(toneMean(pretoneFrames+1:pretoneFrames+peakFrameBin,:),[],1);
switch toneActSel
case 'peak'
    for i = 1:nNeuron
        toneAct(:,:,i) = squeeze(TCpretone_reorder(pretoneFrames+peakFrames(i),:,startTrial:end,i));
    end
case 'sum'
    toneAct = squeeze(nanmean(TCpretone_reorder(pretoneFrames+toneActSumBin,:,startTrial:end,:),1));
end
baseAct = squeeze(nanmean(TCpretone_reorder((pretoneFrames-baselineFrames+1):pretoneFrames,:,startTrial:end,:),1));
baseAct_noAvg = squeeze(TCpretone_reorder((pretoneFrames-baselineFrames+1):pretoneFrames,:,startTrial:end,:));

tuning.toneAct = toneAct; tuning.baseAct = baseAct; tuning.baseAct_noAvg = baseAct_noAvg;
tuning.median = squeeze(median(toneAct,2)); tuning.mean = squeeze(mean(toneAct,2));
magicNum = sqrt(pi/2);tuning.medianSEM = squeeze(nanstd(toneAct,0,2))/sqrt(size(toneAct,2))*magicNum;
tuning.meanSEM = squeeze(nanstd(toneAct,0,2))/sqrt(size(toneAct,2));

%---------SIGNIFICANCE TESTING-----------
tuning.significance = cellfun(@(x,y)(fn_significanceTest(toneAct,baseAct,x,y,toneLabel)),...
    significanceTestList,significanceTestAlpha,'UniformOutput',false);
disp('Significance tests SUCCESSFUL!')
if exist('snrAnalysisList') && contains('roc',snrAnalysisList); tuning.roc = fn_roc(toneAct,baseAct_noAvg); end
disp('SNR analysis SUCCESSFUL!')
if exist('otherAnalysisList') && contains('suppresion',otherAnalysisList); tuning.suppression = fn_suppresion(toneAct,baseAct_noAvg); end
disp('Other analysis SUCCESSFUL!')
%---------SELECT CRITERIA FOR RESPONSIVE CELL-----------
tuning.responsiveCellFlag = tuning.significance{1}.signif;
tuning.responsiveCellToneFlag = tuning.significance{1}.toneH;
%---------COMPUTE POPULATION TUNING DATA-----------
temp = tuning.median; temp(~tuning.responsiveCellToneFlag) = nan; 
popTuning.median = nanmean(temp,2); 
temp = tuning.mean; temp(~tuning.responsiveCellToneFlag) = nan; 
popTuning.mean = nanmean(temp,2); 
popTuning.responsiveCount = sum(tuning.responsiveCellToneFlag,2);
[~,bfIdx] = max(tuning.median,[],1);
popTuning.bfMedian = bfIdx; popTuning.bfCountMedian = histcounts(bfIdx,0.5:1:nTones+0.5);
[~,bfIdx] = max(tuning.mean,[],1);
popTuning.bfMean = bfIdx; popTuning.bfCountMean = histcounts(bfIdx,0.5:1:nTones+0.5);
%---------SAVE RELEVANT DATA-----------
if ~exist('allDataName');  allDataName = {'tuning','popTuning','neuronEachPlane',...
        'peakFrames','TCpretone_reorder','toneLabel'}; end
save([savePath '/population/tuning.mat'],allDataName{:});
disp('Save data SUCCESSFUL!')
%---------LOAD FOV IMAGE FOR TONOTOPY PLOT-----------
neuronPlane = cumsum(neuronEachPlane);
try  neuronPlane = [0;neuronPlane];
catch;  neuronPlane = [0 neuronPlane];disp('check this'); end
for i = 1:nPlanes
    data = load([suite2ppath sep 'plane' num2str(i-1) sep 'Fall.mat']); 
    refImg{i} = data.ops.meanImg;
end
%---------PLOT TONOTOPY AND POPULATION STATISTICS-----------
fn_getTuningPopulation(tuning,popTuning,nNeuron,nPlanes, neuronEachPlane,...
    neuronPlane,toneLabel,toneindex,refImg,roisBound,savePath);

%---------RUN USER DEFINED ANALYTICAL FUNCTION FN2-----------
if ~isempty(fn2); feval(fn2); end 
%---------PLOT SINGLE NEURON-----------
if saveSingleNeuronFlag 
    fn_getTuningSingleNeuron(popTuning,tuning,TCraw,TCpretone_reorder,nFramesPerTrial,...
        nFramesPerTone,nTrials,trialMean,trialMedian,toneAct,baseAct_noAvg,toneLabel,...
        toneindex,nNeuron,nTones,roisBound,refImg,neuronEachPlane,savePath); 
end
end


function testResult = fn_significanceTest(toneAct,baseAct,significanceTest,testAlpha,toneLabel)
    
nNeuron = size(toneAct,3); nTones = size(toneAct,1);
testResult.name = significanceTest; testResult.alpha = testAlpha;
dispFlag = true;

if strcmp(significanceTest,'ttest') || strcmp(significanceTest,'signrank')
    testResult.toneP = zeros(nTones,nNeuron); testResult.toneH = zeros(nTones,nNeuron);
    for i = 1:nNeuron 
        tic;
        for j = 1:nTones
            try
                [h,p] = feval(significanceTest, toneAct(j,:,i),baseAct(j,:,i),'alpha',testAlpha,'tail','right');
                if strcmp(significanceTest,'signrank');temph=h;h=p;p=temph;end
                testResult.toneP(j,i) = p; 
                testResult.toneH(j,i) = h; 
            catch 
                disp(['Error in Cell ' int2str(i) 'Tone ' int2str(j) ' t test']);
                testResult.toneP(j,i) = 1; % set p at 1
                testResult.toneH(j,i) = 0; % set hypothesis wrong
            end          
        end
        t = toc;
        if dispFlag; disp([significanceTest ' -- estimated time = ' num2str(t*nNeuron,'%.2f') ' secs' ]); dispFlag = false; end
    end
    testResult.signif = sum(testResult.toneH,1)>0;
elseif strcmp(significanceTest,'anova')
    nTrials = size(toneAct,2);
    peakANOVA = nan((nTrials)*nTones,nTones+1,nNeuron);
    for i = 1:nNeuron 
        for j = 1:nTones
            % Although most times we are taking 1 frame, but use nanmean here to allow multiple frames
            peakANOVA(((j-1)*(nTrials)+1):(j*(nTrials)),j,i) = squeeze(toneAct(j,:,i));
            peakANOVA(((j-1)*(nTrials)+1):(j*(nTrials)),end,i) = squeeze(baseAct(j,:,i));
        end
    end
    for i= 1:nNeuron
        tic;
        groupNames = cell(1,nTones+1);
        groupNames{end} = 'Baseline';
        for j = 1:nTones
            groupNames{j} = toneLabel{j};
        end
        [p,h,stats] = anova1(peakANOVA(:,:,i),groupNames,'off');
        [results,~,~,~] = multcompare(stats,'Display','off');
        results = results(results(:,2)==(nTones+1),6);
        
        testResult.toneP(:,i) = results;
        testResult.toneH(:,i) = (results<testAlpha);
        testResult.signif(i) = p<0.05 && sum(results<0.05)>0;
        t = toc;
        if dispFlag; disp([significanceTest ' -- estimated time = ' num2str(t*nNeuron,'%.2f') ' secs' ]); dispFlag = false; end
    end
end

end

function results = fn_roc(toneAct,baseAct)
   
nShuffle = 100; nNeuron = size(toneAct,3); nTones = size(toneAct,1);
results.rocAuc = zeros(nTones,nNeuron);results.rocTpr = zeros(nTones,nNeuron);
results.rocAucZscore = zeros(nTones,nNeuron);
dispFlag = true;
for i = 1:nNeuron
    tic;
    for j = 1:nTones
        
        rocBaseAct = squeeze(baseAct(:,j,:,i));
        rocAct = [rocBaseAct(:)' squeeze(toneAct(j,:,i))];
        rocLabel = [zeros(1,numel(rocBaseAct)) ones(1,size(toneAct,2))];
        [tpr, fpr, threshold] = roc(rocLabel, rocAct);
        tempAuc = trapz([0 fpr 1],[0 tpr 1]);
        results.rocAuc(j,i) = tempAuc;
        tempFpr = find(fpr<0.05);
        results.rocTpr(j,i) = tpr(tempFpr(end));
        shuffAuc = zeros(1,nShuffle);
        for k = 1:nShuffle
            shuffLabel = rocLabel(randperm(length(rocLabel)));
            [tprShuff, fprShuff, thresholdShuff] = roc(shuffLabel, rocAct);
            shuffAuc(k) = trapz([0 fprShuff 1],[0 tprShuff 1]);
        end
        results.rocAucZscore(j,i) = (tempAuc - mean(shuffAuc)) / std(shuffAuc);  
    end
    t = toc; 
    if dispFlag; disp(['roc -- estimated time = ' num2str(t*nNeuron,'%.2f') ' secs' ]); dispFlag = false; end
end

end