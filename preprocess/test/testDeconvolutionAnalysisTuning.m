%% run spk session 
clear;
datapath = 'C:\Users\zzhu34\Documents\tempdata\deconv_test_cd036\allSessions\';

load([datapath 'cd036_calman_ar2_foo95_optb_nosmin_s.mat'],'s');

savePath = 'C:\Users\zzhu34\Documents\tempdata\deconv_test_cd036\allSessions\tuning\spk';
sessionSpk = s{2}; sessionSpk = [nan(size(sessionSpk,1),1) sessionSpk];

global nPlanes; nPlanes = 2;
getTuning_oneChannel_noPlot(sessionSpk',savePath,'tuningParam');

%% run dff session 
clear;
%parobj = parpool;
global tempTC ishere nFrames; 
load('C:\Users\zzhu34\Documents\tempdata\excitData\cd036\TC\cd036_TC_plane0.mat','tempTC');
load('C:\Users\zzhu34\Documents\tempdata\excitData\cd036\roi\ishere_plane0.mat','ishere');
ishere = ishere{1};
[nFrames, nFrames_oneplane] = func_readnFrames('cd036','root',...
    'C:\Users\zzhu34\Documents\tempdaawta\excitData\config\mouse');
nFrames = [0 0; nFrames];

for i = 1:size(nFrames,1)-1
    sessionTC = tempTC(:,nFrames(i)+1:nFrames(i+1),1);
    tempSessionDff = (sessionTC - repmat(median(sessionTC,2),1,size(sessionTC,2))) ./ repmat(median(sessionTC,2),1,size(sessionTC,2));
    tempSessionDff = tempSessionDff - smoothdata(tempSessionDff,2,'movmedian',3000);
    sessionDff{i} = tempSessionDff(:,2:end);
end

savePath = 'C:\Users\zzhu34\Documents\tempdata\deconv_test_cd036\allSessions\tuning\dff';
sessionDff = sessionDff{2}; sessionDff = [nan(size(sessionDff,1),1) sessionDff];
global nPlanes; nPlanes = 2;
getTuning_oneChannel_noPlot(sessionDff',savePath,'tuningParam');


%% Data data for comparison
clear;
dff = load('C:\Users\zzhu34\Documents\tempdata\deconv_test_cd036\allSessions\tuning\dff\tuning_smooth3.mat');
spk = load('C:\Users\zzhu34\Documents\tempdata\deconv_test_cd036\allSessions\tuning\spk\tuning_foo90_smooth0_window.mat');

%% Plot significance results for all tones (best tones)
respFlagSpk = spk.anovaPeakCorr(2,:)>0; respFlagDff = dff.anovaPeakCorr(2,:)>0;
plotAllFigures(spk, dff, respFlagSpk,respFlagDff,[]);

%% Plot significance results for target and foil for better comparison

target = 8;%13454;
foil = 11; %22627;

respFlagSpk = spk.anovaSignifToneCorr(target,:)>0; respFlagDff = dff.anovaSignifToneCorr(target,:)>0;
plotAllFigures(spk, dff, respFlagSpk,respFlagDff,target);
respFlagSpk = spk.anovaSignifToneCorr(foil,:)>0; respFlagDff = dff.anovaSignifToneCorr(foil,:)>0;
plotAllFigures(spk, dff, respFlagSpk,respFlagDff,foil);

%% plot the p values for significance test

pTSpk = spk.ttestToneP(target,:); pTDff = dff.ttestToneP(target,:);
figure; histogram(pTDff,0:0.01:0.2)
hold on;histogram(pTSpk,0:0.01:0.2)
%% Plot SI for target and foil 

target = 8;%13454;
foil = 11; %22627;
respFlag = dff.anovaSignifToneCorr(target,:) | dff.anovaSignifToneCorr(foil,:)>0;
plotSI(spk, dff, target, foil, respFlag);

%% Plot single neuron summary for each neuron
savePath = 'C:\Users\zzhu34\Documents\tempdata\deconv_test_cd036\allSessions\tuning\comparison_spk_foo95smin_window_smooth_dff_smooth3\';
plotNeuronSummary(spk, dff,target,foil,savePath);

%% functions
function plotAllFigures(spk, dff, respFlagSpk,respFlagDff,toneIdx)
    if isempty(toneIdx); allToneFlag = true; else; allToneFlag = false; end
    
    bothFlag = (respFlagSpk & respFlagDff);
    dffFlag = (~respFlagSpk & respFlagDff);
    spkFlag = (respFlagSpk & ~respFlagDff);

    fn_plotPieMultiPanel({{bothFlag,dffFlag,spkFlag}}, 'legend', {{'both','dff only','spike only'}},...
            'title',{'responsive cells'});
        
    if ~allToneFlag
        bothPSTH = squeeze(spk.trialMean(:,toneIdx,bothFlag))';
        dffPSTH = squeeze(spk.trialMean(:,toneIdx,dffFlag))';
    else 
        bothPSTH = squeeze(mean(spk.trialMean(:,:,bothFlag),2))';
        dffPSTH = squeeze(mean(spk.trialMean(:,:,dffFlag),2))';
    end
    
    figure; 
    [~, sortIdx ] = sort(max(bothPSTH,[],2),'descend');
    subplot(1,2,1); imagesc(bothPSTH(sortIdx,:)); caxis([prctile(bothPSTH(:),1) prctile(bothPSTH(:),99)]);
    title('Resp cell identified by spikes and dff'); xlabel('frames'); ylabel('neurons')
    [~, sortIdx ] = sort(max(dffPSTH,[],2),'ascend');
    subplot(1,2,2); imagesc(dffPSTH(sortIdx,:)); caxis([prctile(bothPSTH(:),2) prctile(bothPSTH(:),98)]);
    title('Resp cell identified by dff but not spikes'); xlabel('frames'); ylabel('neurons'); 

    if ~allToneFlag
        spkAuc = spk.rocAuc(toneIdx,:);dffAuc = dff.rocAuc(toneIdx,:);
    else
        spkAuc = max(spk.rocAuc,[],1);dffAuc = max(dff.rocAuc,[],1);
    end
    figure;  hold on;
    h = cdfplot(spkAuc(bothFlag)); set( h, 'Color', matlabColors(1), 'LineWidth', 1.5);
    h = cdfplot(spkAuc(dffFlag));set( h, 'Color', matlabColors(1),'LineStyle','--', 'LineWidth', 1.5);
    h = cdfplot(dffAuc(bothFlag)); set( h, 'Color', [0 0 0], 'LineWidth', 1.5);
    h = cdfplot(dffAuc(dffFlag)); set( h, 'Color', [0 0 0], 'LineStyle','--', 'LineWidth', 1.5);
    title('ROC value distribution'); legend('spk-both','spk-dffonly','dff-both','dff-dffonly','Location','Best')
    xlabel('AUC value of ROC (measurement of signal-to-noise)');

    if allToneFlag
        dffTuning = squeeze(mean(dff.peakAct(:,:,respFlagDff),2));dffTuningAll = dffTuning(:); 
        spkTuning = squeeze(mean(spk.peakAct(:,:,respFlagDff),2));spkTuningAll = spkTuning(:);
        tempReg = [spkTuningAll ones(size(spkTuningAll))]; coeff = tempReg \ dffTuningAll; fitSpkTuning  = tempReg*coeff;
        maxLim = max([dffTuningAll spkTuningAll]);minLim = min([dffTuningAll spkTuningAll]);
        
        dffTuningSignifTone = dffTuning; spkTuningSignifTone = spkTuning; respIndex = find(respFlagDff);
        for i = respIndex
            dffTuningSignifTone(~dff.anovaSignifToneCorr(:,i),i) = nan; 
            spkTuningSignifTone(~dff.anovaSignifToneCorr(:,i),i) = nan;
        end
        dffTuningSignifTone = dffTuningSignifTone(:); spkTuningSignifTone = spkTuningSignifTone(:);
        dffTuningSignifTone(isnan(dffTuningSignifTone)) = []; spkTuningSignifTone(isnan(spkTuningSignifTone)) = [];
        tempReg = [spkTuningSignifTone ones(size(spkTuningSignifTone))]; coeff = tempReg \ dffTuningSignifTone; fitSpkTuningSignifTone  = tempReg*coeff;
        figure; subplot(1,3,1);scatter(fitSpkTuningSignifTone,dffTuningSignifTone, 5,'filled');
        hold on; plot([minLim maxLim],[minLim maxLim],...
            'Color',[0.5 0.5 0.5],'LineWidth',2); xlabel('spk tuning (fitted)'); ylabel('dff tuning'); xlim([0 1]);ylim([0 1])
        title(['Tuning (resp tones), spk and dff, corr=' num2str(corr(fitSpkTuningSignifTone,dffTuningSignifTone),'%.2f')]);
        
        subplot(1,3,2);scatter(fitSpkTuning,dffTuningAll, 5,'filled'); hold on; plot([minLim maxLim],[minLim maxLim],...
            'Color',[0.5 0.5 0.5],'LineWidth',2); xlabel('spk tuning (fitted)'); ylabel('dff tuning'); xlim([0 1]);ylim([0 1])
        title(['Tuning, spk and dff, corr=' num2str(corr(fitSpkTuning,dffTuningAll),'%.2f')]);

        dffRoc = dff.rocAuc(:,respFlagDff); spkRoc = spk.rocAuc(:,respFlagDff);
        subplot(1,3,3);scatter(spkRoc(:),dffRoc(:), 10,'filled'); hold on; plot([minLim maxLim],[minLim maxLim],...
            'Color',[0.5 0.5 0.5],'LineWidth',2); xlabel('spk ROC'); ylabel('dff ROC');title('Distribution of spk and dff SI');
        title(['ROC, spk and dff, corr=' num2str(corr(spkRoc(:),dffRoc(:)),'%.2f')]);xlim([0 1]);ylim([0 1])
        
        for i = 1:size(dffTuning,2);neuronCorr(i) = corr(dffTuning(:,i),spkTuning(:,i));end
        figure;histogram(neuronCorr,50); xlabel('Tuning Corr');ylabel('Frequency')
        
        [~,dffBestFreq] = max(dffTuning); [~,spkBestFreq] = max(spkTuning); 
        figure; histogram(spkBestFreq - dffBestFreq,-17:1:17);
        xlabel('spk - dff'); title(['Difference in best freq. Mean Abs Diff = ' ...
            num2str(mean(abs(spkBestFreq - dffBestFreq)),'%.1f')]);
    end

end

function plotSI(spk, dff, targIdx, foilIdx, respFlag)
    
    spkPeakT = squeeze(mean(spk.peakAct(targIdx,:,:),2));
    spkPeakF = squeeze(mean(spk.peakAct(foilIdx,:,:),2));
    dffPeakT = squeeze(mean(dff.peakAct(targIdx,:,:),2));
    dffPeakF = squeeze(mean(dff.peakAct(foilIdx,:,:),2));
    spkSI = abs(spkPeakT - spkPeakF) ./ (abs(spkPeakT) + abs(spkPeakF));
    dffSI = abs(dffPeakT - dffPeakF) ./ (abs(dffPeakT) + abs(dffPeakF));

    figure; hold on
    h = cdfplot(spkSI(respFlag)); set( h, 'Color', matlabColors(1), 'LineWidth', 1.5);
    h = cdfplot(spkSI(~respFlag)); set( h, 'Color', matlabColors(1),'LineStyle','--', 'LineWidth', 1.5);
    h = cdfplot(dffSI(respFlag));set( h, 'Color', [0 0 0],'LineStyle','-', 'LineWidth', 1.5);
    h = cdfplot(dffSI(~respFlag));set( h, 'Color', [0 0 0],'LineStyle','--', 'LineWidth', 1.5);
    xlabel('SI'); ylabel('proportion of cells'); title('Distribution of SI')
    legend('spk-resp','spk-unresp','dff-resp','dff-unresp','Location','Best')

    deltaSI = spkSI - dffSI;
    fn_plotHistMultiPanel({{deltaSI(respFlag),deltaSI(~respFlag)}},'histCountArgIn',{-2:0.1:2},...
        'plotArgIn',{'LineWidth',2});
    ylimm = ylim; plot([0 0], ylimm,'Color',[0.5 0.5 0.5]);
    legend('responsive cells','unresponsive cells'); xlabel('spkSI - dffSI'); ylabel('proportion of cells');
    title('deltaSI of responsive and unresponsive neurons compared')

    figure; scatter(spkSI(respFlag), dffSI(respFlag),10,'filled'); hold on; plot([0 1],[0 1],...
        'Color',[0.5 0.5 0.5],'LineWidth',2); xlabel('spk SI'); ylabel('dff SI');
    title(['SI, spk and dff, corr=' num2str(corr(spkSI(respFlag),dffSI(respFlag)),'%.2f')]);
    lmModel = fitlm(spkSI(respFlag), dffSI(respFlag));


    bin_deltaSI = [-1 -0.5:0.2:0.5 1]; bin_deltaSIAvg = (bin_deltaSI(1:end-1) + bin_deltaSI(2:end))/2;
    [N, edges, nbin] = histcounts(deltaSI(respFlag), bin_deltaSI);
    spkPSTH = squeeze(spk.trialMean(:,targIdx,respFlag))';
    dffPSTH = squeeze(dff.trialMean(:,targIdx,respFlag))';
    fn_figureWholeScreen(); timeBin = 6:30;
    for i = 1:length(N)
        tempIdx = (nbin==i); tempSpk = spkPSTH(tempIdx,timeBin);  tempDff = dffPSTH(tempIdx,timeBin); 
        [~, sortIdx] = sort(max(tempDff(:,6:16),[],2),'descend');
        subplot_tight(2,length(N),i); imagesc(tempSpk(sortIdx,:)); set(gca,'xtick',[]); 
        caxis([prctile(spkPSTH(:),5) prctile(spkPSTH(:),99.9)]); if i==1; ylabel('neuron'); end;
        title(['deltaSI ' num2str(bin_deltaSI(i),'%.1f') ' to ' num2str(bin_deltaSI(i+1),'%.1f') ]);
        subplot_tight(2,length(N),length(N)+i);imagesc(tempDff(sortIdx,:)); 
        caxis([prctile(dffPSTH(:),5) prctile(dffPSTH(:),99.5)]);if i==1; ylabel('neuron'); end;
        xlabel('frames(onset=5)')
    end

    spkPSTH = squeeze(spk.trialMean(:,foilIdx,:))';
    dffPSTH = squeeze(dff.trialMean(:,foilIdx,:))';
    fn_figureWholeScreen(); timeBin = 6:30;
    for i = 1:length(N)
        tempIdx = (nbin==i); tempSpk = spkPSTH(tempIdx,timeBin);  tempDff = dffPSTH(tempIdx,timeBin); 
        [~, sortIdx] = sort(max(tempDff(:,6:16),[],2),'descend');
        subplot_tight(2,length(N),i); imagesc(tempSpk(sortIdx,:)); set(gca,'xtick',[]); 
        caxis([prctile(spkPSTH(:),5) prctile(spkPSTH(:),99.9)]);if i==1; ylabel('neuron'); end
        title(['deltaSI ' num2str(bin_deltaSI(i),'%.1f') ' to ' num2str(bin_deltaSI(i+1),'%.1f') ]);
        subplot_tight(2,length(N),length(N)+i);imagesc(tempDff(sortIdx,:)); 
        caxis([prctile(dffPSTH(:),5) prctile(dffPSTH(:),99.5)]);if i==1; ylabel('neuron'); end
        xlabel('frames (onset=5)')
    end
end

function plotNeuronSummary(spk, dff,targIdx,foilIdx,savePath)
    mkdir(savePath);
    spkPeakT = squeeze(mean(spk.peakAct(targIdx,:,:),2));
    spkPeakF = squeeze(mean(spk.peakAct(foilIdx,:,:),2));
    dffPeakT = squeeze(mean(dff.peakAct(targIdx,:,:),2));
    dffPeakF = squeeze(mean(dff.peakAct(foilIdx,:,:),2));
    spkSI = abs(spkPeakT - spkPeakF) ./ (abs(spkPeakT) + abs(spkPeakF));
    dffSI = abs(dffPeakT - dffPeakF) ./ (abs(dffPeakT) + abs(dffPeakF));

    respFlagSpkT = spk.anovaSignifToneCorr(targIdx,:)>0; 
    respFlagDffT = dff.anovaSignifToneCorr(targIdx,:)>0;
    respFlagSpkF = spk.anovaSignifToneCorr(foilIdx,:)>0; 
    respFlagDffF = dff.anovaSignifToneCorr(foilIdx,:)>0;

    columns = 6;
    for i = 1:size(spk.TC_original,2)
        f = fn_figureWholeScreen('visible','off'); 
        nFramesPerTone = 50; pretoneFrames = 10;
        frameAxis = pretoneFrames:10:nFramesPerTone;
        frameLabel = cellfun(@num2str,num2cell(frameAxis-pretoneFrames),'UniformOutput',false);
        
        % plot 1 - PSTH of each tone
        subplot_tight(2,columns,1); temp = dff.trialMean(:,:,i); imagesc(temp'); 
        caxis([prctile(temp(:),2) prctile(temp(:),99)]); ylabel('Tones(dff)'); title('PSTH, all tones')
        subplot_tight(2,columns,1+columns); temp = spk.trialMean(:,:,i); imagesc(temp'); 
        xticklabels([])
        caxis([prctile(temp(:),2) prctile(temp(:),99.5)]); ylabel('Tones(spk)');xlabel('Frames');
        xticks(frameAxis); xticklabels(frameLabel)
        
        % plot 2 - Peak of each tone and each trial
        subplot_tight(2,columns,2); temp = dff.peakAct(:,:,i); imagesc(temp); 
        caxis([prctile(temp(:),2) prctile(temp(:),99)]); ylabel('Tones(dff)'); title('Peak activity, all trial')
        xticklabels([])
        subplot_tight(2,columns,2+columns); temp = spk.peakAct(:,:,i); imagesc(temp); 
        caxis([prctile(temp(:),2) prctile(temp(:),99.5)]); ylabel('Tones(spk)');xlabel('Trials');

        % plot 3 - Target PSTH, each trial
        subplot_tight(2,columns,3); temp = squeeze(dff.TCpretone_reorder(:,targIdx, :,i)); imagesc(temp'); 
        caxis([prctile(temp(:),2) prctile(temp(:),99)]); ylabel('Trials(dff)'); title('Target PSTH, all trial')
        xticklabels([])
        subplot_tight(2,columns,3+columns); temp = squeeze(spk.TCpretone_reorder(:,targIdx, :,i)); imagesc(temp'); 
        caxis([prctile(temp(:),2) prctile(temp(:),99.5)]); ylabel('Trials(spk)');xlabel('Frames');
        xticks(frameAxis);xticklabels(frameLabel)

        % plot 4 - Foil PSTH, each trial
        subplot_tight(2,columns,4); temp = squeeze(dff.TCpretone_reorder(:,foilIdx, :,i)); imagesc(temp'); 
        caxis([prctile(temp(:),2) prctile(temp(:),99)]); ylabel('Trials(dff)'); title('Foil PSTH, all trial')
        xticklabels([])
        subplot_tight(2,columns,4+columns); temp = squeeze(spk.TCpretone_reorder(:,foilIdx, :,i)); imagesc(temp'); 
        caxis([prctile(temp(:),2) prctile(temp(:),99.5)]); ylabel('Trials(spk)');xlabel('Frames');
        xticks(frameAxis);xticklabels(frameLabel)

        % plot 5 - Target and foil PSTH, averge
        tempT = squeeze(mean(dff.TCpretone_reorder(:,targIdx, :,i),3)); 
        tempF = squeeze(mean(dff.TCpretone_reorder(:,foilIdx, :,i),3)); 
        plotWindow = 6:30;
        subplot_tight(2,columns,5); plot(tempT,'LineWidth',2); hold on; plot(tempF,'LineWidth',2); 
        xlim([1 length(tempF)]); ylimm = ylim; plot([10.5 10.5], ylimm,'Color',[0.5 0.5 0.5])
        ylabel('dff'); title(['Dff PSTH, SI= ' num2str(dffSI(i),'%.2f')]);
        xticklabels([]); legend('Target','Foil'); ylim(ylimm)
        
        tempT = squeeze(mean(spk.TCpretone_reorder(:,targIdx, :,i),3)); 
        tempF = squeeze(mean(spk.TCpretone_reorder(:,foilIdx, :,i),3)); 
        subplot_tight(2,columns,5+columns); plot(tempT,'LineWidth',2); hold on; plot(tempF,'LineWidth',2); 
        xlim([1 length(tempF)]);ylimm = ylim; plot([10.5 10.5], ylimm,'Color',[0.5 0.5 0.5])
        ylabel('spk'); title(['SPK PSTH, SI= ' num2str(spkSI(i),'%.2f')]); xlabel('Frames');
        xticks(frameAxis);xticklabels(frameLabel); legend('Target','Foil'); ylim(ylimm)

        % plot 6 - Tuning curve comparisona nd ROC comparison
        dffTuning = mean(dff.peakAct(:,:,i),2); spkTuning = mean(spk.peakAct(:,:,i),2);
        tempReg = [spkTuning ones(size(spkTuning))]; coeff = tempReg \ dffTuning; 
        subplot_tight(2,columns,6); plot(dffTuning,'LineWidth',2); hold on; plot(tempReg*coeff,'LineWidth',2); 
        ylimm = ylim; plot([targIdx targIdx], ylimm,'Color',[0 0.5 0]); plot([foilIdx foilIdx], ylimm,'Color',[0.5 0 0])
        legend('dff','spk'); title('Tuning Curve'); ylim(ylimm)
        xlim([1 length(dffTuning)]);xticklabels([]);ylabel('dff');

        dffRoc = dff.rocAuc(:,i); spkRoc = spk.rocAuc(:,i);
        subplot_tight(2,columns,6+columns); plot(dffRoc,'LineWidth',2); hold on; plot(spkRoc,'LineWidth',2); 
        ylimm = ylim; plot([targIdx targIdx], ylimm,'Color',[0 0.5 0]); plot([foilIdx foilIdx], ylimm,'Color',[0.5 0 0])
        legend('dff','spk');title('ROC curve'); ylim(ylimm)
        xlim([1 length(dffTuning)]); xlabel('Tones'); ylabel('AUC')


        if respFlagSpkT(i) ~= respFlagDffT(i); colorT = {1 0 0}; else; colorT = {0 0 0}; end
        if respFlagSpkF(i) ~= respFlagDffF(i); colorF = {1 0 0}; else; colorF = {0 0 0}; end
        
        if respFlagDffT(i); txtT = 'T-Resp'; else; txtT = 'T-NotResp'; end
        if respFlagDffF(i); txtF = 'F-Resp'; else; txtF = 'F-NotResp'; end
        
        titleT = sprintf('\\color[rgb]{%f, %f, %f}%s', colorT{:},txtT);
        titleF = sprintf('\\color[rgb]{%f, %f, %f}%s', colorF{:},txtF);
        suptitle(['Neuron ' int2str(i)  ' ' titleT ' ' titleF]);
        saveas(f,[ savePath '/cell' int2str(i) '.png']);
        close(f);
    end
end

