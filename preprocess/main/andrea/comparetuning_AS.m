function tuningidxwt = comparetuning_AS()
close all
clear all
[~,~,table]=xlsread('Z:\Andrea2\2P\CompareTuningDataSets\Tuning_curve6v12mAPP.xlsx');
% [~,~,table] = xlsread('home/jennifer/Dropbox/mfiles/Mice_Behavior/miscinfos/Tuning_curve.xlsx');

% animallist = vertcat(table(2:end,1));
% % animallist = {table{2:end,1}};
% 
% genolist =  vertcat(table(2:end,4));
% applist = find(strcmp(genolist,'APP+'));
% wtlist = find(strcmp(genolist,'WT'));
% datalist = vertcat(table{2:end,2});
% % datalist = table{2:end,2};
% filelist = {table(2:end,3)};
% filelist = filelist{1};
animallist = vertcat(table{2:end,1});
genolist = {table{2:end,4}};
applist = find(strcmp(genolist,'APP+'));
% appKIlist = find(strcmp(genolist,'APPKI'));
wtlist = find(strcmp(genolist,'WT'));
datalist = vertcat(table{2:end,2});
filelist = {table{2:end,3}};
nTones = 17;
nTrials = 10;
nFramesPerTone = 100/2; % 50
nFramesPerTrial = nFramesPerTone * nTones; % 850
startTrial = 2; % the first tone is on frame 0
nFrames = nFramesPerTrial*nTrials; % 4250
frameRate = 31.25/2; % in the future, do not hard code this.
pretoneFrames = 10;
baselineFrames = 5;
smoothWindow = 5;
toneOnset = 20/2; % should be 0 for Celine/Ziyi TDT protocol, and 20 for Aaron and Fangchen
tonePeakLim = 0.66;
RootAddress1 = ['Z:\Andrea2\2P\APP_PS1\AS037\preprocessData\'];
tuning=load([RootAddress1 'AS037_000_000_tuning\population\tuning.mat']);
responsiveCellFlag = tuning.anovaPeakCorr(2,:);

trialMean1=tuning.trialMean;
trialMedian1=tuning.trialMedian;

for i = 1:17
  [maxValue,tempPeakIdx] = max(squeeze(trialMean1(pretoneFrames+1:pretoneFrames+ceil(frameRate*tonePeakLim),i,:)),[],1); 
  peakActAvg(:,i) = maxValue; peakFrames(:,i) =  tempPeakIdx; 
end

for ani = [1:length(animallist)]
    clear tuningMean tuningMedian responsiveCellFlag tuningData
        disp(ani)
%         datalist(ani)
        if datalist(ani) == 2
%             RootAddress = ['X:\LabData', num2str(datalist(ani)), '\Kelly\Behavior_Imaging\', animallist{ani}, '\preprocessData\'];
            RootAddress = ['X:\LabData' num2str(datalist(ani)) '\Kelly\Behavior_Imaging\' animallist(ani,1:end) '\preprocessData\'];
%             file= filelist{ani,1};
            %file2= file(ani);
%             fullpath_ = fullfile(RootAddress, 'TuningCurve_', file, '/population/tuningMean.mat');
%             load([RootAddress 'TuningCurve_' file '/population/tuningMean.mat']);
            load([RootAddress 'TuningCurve_' filelist{ani} '/population/tuningMean.mat']);
            load([RootAddress 'TuningCurve_' filelist{ani} '/population/tuningMedian.mat']);
            load([RootAddress 'TuningCurve_' filelist{ani} '/population/responsiveCellFlag.mat']);
%             load([RootAddress 'TuningCurve_' filelist{ani} '/population/neuronRespTable.mat'])
            % load([RootAddress 'TuningCurve_' filelist{ani} '/population/tuningData.mat'])
        elseif datalist(ani) == 3
            RootAddress = ['W:\LabData' num2str(datalist(ani)) '\Kelly\' animallist(ani,1:end) '\preprocessData\'];
            load([RootAddress 'TuningCurve_' filelist{ani} '/population/tuningMean.mat']);
            load([RootAddress 'TuningCurve_' filelist{ani} '/population/tuningMedian.mat']);
            load([RootAddress 'TuningCurve_' filelist{ani} '/population/responsiveCellFlag.mat']);
        elseif datalist(ani) == 1
%         RootAddress = ['home/jennifer/Dropbox/mfiles/Mice_Behavior/' animallist(ani,1:end) '/'];
%             RootAddress = ['Z:\Andrea2\2P\APP_NLGF\' animallist(ani,1:end) '\preprocessData\'];
            RootAddress = ['Z:\Andrea2\2P\APP_PS1\' animallist(ani,1:end) '\preprocessData\'];
            tuning=load([RootAddress '' filelist{ani} '_tuning\population\tuning.mat']);
            responsiveCellFlag = tuning.anovaPeakCorr(2,:);
%             tuningMean=tuning.poptuningMean;
%             tuningMedian=tuning.poptuningMedian;
            trialMean=tuning.trialMean;
            nNeuron=size(trialMean,3);
            trialMedian=tuning.trialMedian;
            for i = 1:nNeuron
                tuningMean(:,i)= squeeze(trialMean(pretoneFrames+peakFrames(i),:,i));
                tuningMedian(:,i) = squeeze(trialMedian(pretoneFrames+peakFrames(i),:,i));
            end
         else 
%         RootAddress = ['home/jennifer/Dropbox/mfiles/Mice_Behavior/' animallist(ani,1:end) '/'];
            RootAddress = ['Z:\Andrea2\2P\APP_NLGF\' animallist(ani,1:end) '\preprocessData\'];
%             RootAddress = ['Z:\Andrea2\2P\APP_PS1\' animallist(ani,1:end) '\preprocessData\'];
            tuning=load([RootAddress '' filelist{ani} '_tuning\population\tuning.mat']);
            responsiveCellFlag = tuning.anovaPeakCorr(2,:);
%             tuningMean=tuning.poptuningMean;
%             tuningMedian=tuning.poptuningMedian;
            trialMean=tuning.trialMean;
            nNeuron=size(trialMean,3);
            trialMedian=tuning.trialMedian;
            for i = 1:nNeuron
                tuningMean(:,i)= squeeze(trialMean(pretoneFrames+peakFrames(i),:,i));
                tuningMedian(:,i) = squeeze(trialMedian(pretoneFrames+peakFrames(i),:,i));
            end

        end
%     RootAddress = ['Z:\Andrea2\2P\APP_PS1\' animallist(ani,1:end) '\preprocessData\'];
% %     tuning=load([RootAddress '' filelist{ani} '_tuning\population\tuning.mat']);
%     tuning=load([RootAddress '' filelist{ani} '_tuning\population\tuning.mat']);
%     responsiveCellFlag = tuning.anovaPeakCorr(2,:);
% %             tuningMean=tuning.poptuningMean;
% %             tuningMedian=tuning.poptuningMedian;
%     trialMean=tuning.trialMean;
%     nNeuron=size(trialMean,3);
%     trialMedian=tuning.trialMedian;
%     for i = 1:nNeuron
%        tuningMean(:,i)= squeeze(trialMean(pretoneFrames+peakFrames(i),:,i));
%        tuningMedian(:,i) = squeeze(trialMedian(pretoneFrames+peakFrames(i),:,i));
%     end
    Meantuning{ani}=tuningMean;
    Mediantuning{ani}=tuningMedian;
    TunedNeurons{ani}=logical(responsiveCellFlag);
%     ToneEvokedMat{ani} = tuningData;
%     Significantones{ani} = neuronRespTable(1:17,TunedNeurons{ani})
%     Bandwidth_Large{ani} = find(Meantuning{ani}(:,TunedNeurons{ani}));%;cellfun(@max,Significantones(TunedNeurons))]%sum(neuronRespTable(1:end-1,TunedNeurons+1));% refine this, how about nb of octive it's responsive to ?
%     Qfactor = BF'./Bandwidth_Large; % BF/Bandwidth according to Yue 2014 %inf is no spread since log2(1) = 0

end

PopTuningMean = cell2mat(cellfun(@(x,y) mean(x(:,y),2),Meantuning,TunedNeurons,'UniformOutput',0));
PopTuningMedian = cell2mat(cellfun(@(x,y) mean(x(:,y),2),Mediantuning,TunedNeurons,'UniformOutput',0));
NormPopTuningMean = (PopTuningMean-min(PopTuningMean))./max((PopTuningMean-min(PopTuningMean)));
NormPopTuningMedian = (PopTuningMedian-min(PopTuningMedian))./max((PopTuningMedian-min(PopTuningMedian)));
NormTuningNeuronstmp = cellfun(@(x,y) (x(:,y)-min(x(:,y)))./max(x(:,y)-min(x(:,y))),Meantuning,TunedNeurons,'UniformOutput',0);
NormTuningNeurons = cell2mat(NormTuningNeuronstmp)';
indexneuronmouse = cell2mat(cellfun(@(x,y) ones(1,x)*y,mat2cell(cellfun(@(x) size(x,2),NormTuningNeuronstmp),1,ones(1,length(animallist))),num2cell([1:length(animallist)]),'UniformOutput',false));
[neuronnbbf,indbfneurons] = find(NormTuningNeurons == 1);
nboftunedneurons = cellfun(@(x) numel(find(x==1)),TunedNeurons);
nbofroi = cellfun(@length,TunedNeurons);
prctevoked = nboftunedneurons./nbofroi;
for i = 1:length(NormTuningNeurons)
    Bandwidthneurons(i) = numel(find(NormTuningNeurons(i,:)>0.3));
    Bandwidthind{i} = (find(NormTuningNeurons(i,:)>0.3));
end

[Bandwithvalue,Sortedcellsind] = sort(Bandwidthneurons);
%center neurons around bf
CenterNormTuningNeurons = nan(size(NormTuningNeurons,1),(size(NormTuningNeurons,2)*2)-1);
for i = 1:size(NormTuningNeurons,1)
     CenterNormTuningNeurons(i,(max(indbfneurons)-(find(NormTuningNeurons(i,:)==1)):max(indbfneurons)+(16-find(NormTuningNeurons(i,:)==1)))+1) = NormTuningNeurons(i,:);
end
% for i = 1:size(NormTuningNeurons,1)
%     CenterNormTuningNeurons(i,:) = s;
% end 

toneorder = [45254.834 8000 13454.34264 4756.82846 5656.854249,...
        22627.417 64000 53817.37058 4000 9513.65692,...
        16000 6727.171322 19027.31384 26908.6852 32000,...
        11313.7085 38054.62768];
freqAxis = log2(sort(toneorder));

% c = [0.8 0.8 0.8; 1 1 1; 1 0.8 0.8 ; 1 0.7 0.7; 1 0.6 0.6; 1 0.5 0.5; 1 0.4 0.4; 1 0.3 0.3; 1 0.2 0.2; 1 0.1 0.1; 1 0 0; 0.8 0 0; 0.6 0 0];

c = [0.8 0.8 0.8]
logvalue = 1.22.^(1:20)/1.22^20;%log([1:1.7/100:2.7]);
% rmlog = flip((logvalue));
for cc = 1:length(logvalue)
   c = [c; [1 1-logvalue(cc) 1-logvalue(cc)]];
end
% c = [c; 0 0 0];

appneurons = find(ismember(indexneuronmouse,applist));
wtneurons = find(ismember(indexneuronmouse,wtlist));

tuningxvalues = [17-4:17+4];
tuningxvalues(find(tuningxvalues == 17)) = [];

h = figure(2); hold on; 
subplot(1,2,1); hold on; title('Tuning Curve APP+ n = 6')
% [nr,nc] = size(CenterNormTuningNeurons(ismember(indexneuronmouse,applist),:));
% pcolor([CenterNormTuningNeurons(ismember(indexneuronmouse,applist),:) nan(nr,1); nan(1,nc+1)]);
% shading flat;
xvalues = [-1:0.25:1];
xvalues(5)= [];
imagesc(xvalues,1:numel(appneurons),CenterNormTuningNeurons(Sortedcellsind(find(ismember(Sortedcellsind,appneurons))),tuningxvalues));
set(gca, 'ydir', 'reverse');
set(gca, 'Color')
xlabel('Tone distance from bf');
ylabel('Cell number');
ylim([1 sum(ismember(indexneuronmouse,applist))]);
colormap(c);
xlim([xvalues(1) xvalues(end)])
xticks([-1 -0.5 0 0.5 1])
xticklabels({'-1','-1/2','0','1/2','1'})
subplot(1,2,2); hold on; title('Tuning Curve APP- n = 3')
% [nr,nc] = size(CenterNormTuningNeurons(ismember(indexneuronmouse,wtlist),:));
% pcolor([CenterNormTuningNeurons(ismember(indexneuronmouse,wtlist),:) nan(nr,1); nan(1,nc+1)]);
% shading flat;
imagesc(xvalues,1:numel(wtneurons),CenterNormTuningNeurons(Sortedcellsind(find(ismember(Sortedcellsind,wtneurons))),tuningxvalues));
set(gca, 'ydir', 'reverse');
set(gca, 'Color')                                                                                                                                            
xlabel('Tone distance from bf [oct.]');
ylabel('Cell number');
ylim([1 sum(ismember(indexneuronmouse,wtlist))]);
colormap(c);
xlim([xvalues(1) xvalues(end)])
xticks([-1 -0.5 0 0.5 1])
xticklabels({'-1','-1/2','0','1/2','1'})
% Tuningsorted = ;
% BandwidthSorted = ;


% figure(1); hold on; 
% for i = 1:length(applist)
% subplot(1,2,1); hold on;
% title('APP mean tuning activity'); plot(freqAxis,PopTuningMean(:,applist(i))-1);
% end
% set(gca, 'XTick', log2([4000 8000 16000 32000 64000]));
% set(gca, 'XTickLabel', [4 8 16 32 64]);
% xlabel('Frequency (kHz)');
% ylabel('mean F/F0');
% 
% for i = 1:length(wtlist)
% subplot(1,2,2); hold on;
% title('WT mean tuning activity'); plot(freqAxis,PopTuningMean(:,wtlist(i))-1);
% end
% set(gca, 'XTick', log2([4000 8000 16000 32000 64000]));
% set(gca, 'XTickLabel', [4 8 16 32 64]);
% xlabel('Frequency (kHz)');
% ylabel('mean F/F0');
% suptitle('Mean tuning activity APP+ and WT')
% % saveas(gcf,['X:\LabData2\kelly\comparetuningmean'],'pdf');
% 
% figure(2); hold on; 
% for i = 1:length(applist)
% subplot(1,2,1); hold on;
% title('APP median tuning activity'); plot(freqAxis,PopTuningMedian(:,applist(i))-1);
% end
% set(gca, 'XTick', log2([4000 8000 16000 32000 64000]));
% set(gca, 'XTickLabel', [4 8 16 32 64]);
% xlabel('Frequency (kHz)');
% ylabel('mean F/F0');
% 
% for i = 1:length(wtlist)
% subplot(1,2,2); hold on;
% title('WT median tuning activity'); plot(freqAxis,PopTuningMedian(:,wtlist(i))-1);
% end
% set(gca, 'XTick', log2([4000 8000 16000 32000 64000]));
% set(gca, 'XTickLabel', [4 8 16 32 64]);
% xlabel('Frequency (kHz)');
% ylabel('mean F/F0');
% suptitle('Median tuning activity APP+ and WT')
% % saveas(gcf,['X:\LabData2\kelly\comparetuningmedian'],'pdf');

% plot center of bf
[~,bfmeanind] = max(PopTuningMean);
[~,bfmedianind] = max(PopTuningMedian);

%reshape mat
meannewmattuning = nan(length(animallist),50);%max(bfmeanind)-min(bfmeanind)+17);
mediannewmattuning = nan(length(animallist),50);%max(bfmeanind)-min(bfmeanind)+17);

for i = 1:length(animallist)
 meannewmattuning(i,(max(bfmeanind)-(bfmeanind(i)):max(bfmeanind)+(16-bfmeanind(i)))+1) = NormPopTuningMean(:,i);
 mediannewmattuning(i,(max(bfmedianind)-(bfmedianind(i)):max(bfmedianind)+(16-bfmedianind(i)))+1) = NormPopTuningMedian(:,i);
end
MeanNormTuningnewmat = meannewmattuning(:,1:27);
MedianNormTuningnewmat = mediannewmattuning(:,1:27);
% figure(3); hold on; 
% for i = 1:length(applist)
% subplot(1,2,1); hold on;
% xvalues = [1:17]-bfmeanind(applist(i));
% title('APP median tuning activity'); plot(xvalues,PopTuningMedian(:,applist(i))-1);
% end
% % set(gca, 'XTick', log2([4000 8000 16000 32000 64000]));
% % set(gca, 'XTickLabel', [4 8 16 32 64]);
% xlabel('Frequency (kHz)');
% ylabel('mean F/F0');
% 
% for i = 1:length(wtlist)
% subplot(1,2,2); hold on;
% xvalues = [1:17]-bfmeanind(wtlist(i));
% title('WT median tuning activity'); plot(xvalues,PopTuningMedian(:,wtlist(i))-1);
% end
% % set(gca, 'XTick', log2([4000 8000 16000 32000 64000]));
% % set(gca, 'XTickLabel', [4 8 16 32 64]);
% xlabel('Frequency (kHz)');
% ylabel('mean F/F0');
% suptitle('Median tuning activity APP+ and WT')
% saveas(gcf,['X:\LabData2\kelly\centercomparetuningmedian'],'pdf');

% figure(4); hold on; 
% for i = 1:length(applist)
% subplot(1,2,1); hold on;
% xvalues = [1:17]-bfmeanind(applist(i));
% title('APP mean tuning activity'); plot(xvalues,PopTuningMean(:,applist(i))-1);
% end
% % set(gca, 'XTick', log2([4000 8000 16000 32000 64000]));
% % set(gca, 'XTickLabel', [4 8 16 32 64]);
% xlabel('Frequency (kHz)');
% ylabel('mean F/F0');
% 
% for i = 1:length(wtlist)
% subplot(1,2,2); hold on;
% xvalues = [1:17]-bfmeanind(wtlist(i));
% title('WT mean tuning activity'); plot(xvalues,PopTuningMean(:,wtlist(i))-1);
% end
% % set(gca, 'XTick', log2([4000 8000 16000 32000 64000]));
% % set(gca, 'XTickLabel', [4 8 16 32 64]);-
% xlabel('Frequency (kHz)');
% ylabel('mean F/F0');
% suptitle('Mean tuning activity APP+ and WT')
% saveas(gcf,['X:\LabData2\kelly\centercomparetuningmean'],'pdf');
h = figure(4); hold on; 
plot(-13:13,nanmean(MedianNormTuningnewmat(wtlist,:)), 'LineStyle', '-','Color','r','LineWidth',2);
plot(-13:13,nanmean(MedianNormTuningnewmat(applist,:)), 'LineStyle', '--','Color','r','LineWidth',2);
xlim([-13 13]);
xlabel('Number of tone relative to best frequency');
ylabel('median F/F0 [norm.]');
legend('APP/PS1+  (6m n=12)', 'APP/PS1+ (12m n=3)'); legend boxoff

figure(5); hold on; 
for i = 1:length(applist)
subplot(1,2,1); hold on;
xvalues = [1:17]-bfmeanind(applist(i));
title('APP mean tuning activity'); plot(xvalues,NormPopTuningMean(:,applist(i)));
end
plot(-13:13,nanmean(MeanNormTuningnewmat(applist,:)),'Color',[0,0,0],'LineWidth',2);
set(gca, 'XTick', log2([4000 8000 16000 32000 64000]));
set(gca, 'XTickLabel', [4 8 16 32 64]);
xlabel('Frequency (kHz)');
ylabel('mean F/F0');

for i = 1:length(wtlist)
subplot(1,2,2); hold on;
xvalues = [1:17]-bfmeanind(wtlist(i));
title('WT mean tuning activity'); plot(xvalues,NormPopTuningMean(:,wtlist(i)));
end
plot(-13:13,nanmean(MeanNormTuningnewmat(wtlist,:)),'Color',[0,0,0],'LineWidth',2);

% set(gca, 'XTick', log2([4000 8000 16000 32000 64000]));
% set(gca, 'XTickLabel', [4 8 16 32 64]);
xlabel('Frequency (kHz)');
ylabel('mean F/F0');
%suptitle('Mean tuning activity APP+ and WT')
% saveas(gcf,['X:\LabData2\kelly\centercomparetuningmean'],'pdf');

figure(6); hold on; 
for i = 1:length(applist)
subplot(1,2,1); hold on;
xvalues = [1:17]-bfmeanind(applist(i));
title('APP median tuning activity'); plot(xvalues,NormPopTuningMedian(:,applist(i)));
end
plot(-13:13,nanmean(MedianNormTuningnewmat(applist,:)),'Color',[0,0,0],'LineWidth',2);
% set(gca, 'XTick', log2([4000 8000 16000 32000 64000]));
% set(gca, 'XTickLabel', [4 8 16 32 64]);
xlabel('Frequency (kHz)');
ylabel('mean F/F0');

for i = 1:length(wtlist)
subplot(1,2,2); hold on;
xvalues = [1:17]-bfmeanind(wtlist(i));
title('WT median tuning activity'); plot(xvalues,NormPopTuningMedian(:,wtlist(i)));
end
plot(-13:13,nanmean(MedianNormTuningnewmat(wtlist,:)),'Color',[0,0,0],'LineWidth',2);

% set(gca, 'XTick', log2([4000 8000 16000 32000 64000]));
% set(gca, 'XTickLabel', [4 8 16 32 64]);
xlabel('Frequency (kHz)');
ylabel('mean F/F0');
%suptitle('Median tuning activity APP+ and WT')
% saveas(gcf,['X:\LabData2\kelly\centercomparetuningmedian'],'pdf');

% plot dcf for bandwith as normalized tuning curve fro wt and app
for ani = 1:length(animallist)
 animap{ani} = CenterNormTuningNeurons(find(indexneuronmouse==ani),:);
end
appmap=CenterNormTuningNeurons(appneurons,:);
wtmap=CenterNormTuningNeurons(wtneurons,:);
tuningidxapp = sort(nansum(appmap,2));
tuningidxwt = sort(nansum(wtmap,2));
tuningidani = (cellfun(@(x) sort(nansum(x,2)),animap,'UniformOutput',false));
xwt=1/length(tuningidxwt):1/length(tuningidxwt):1;
xapp=1/length(tuningidxapp):1/length(tuningidxapp):1;

figure;
plot(xapp,tuningidxapp,'r');
hold on;
plot(xwt,tuningidxwt,'b');


figure(); hold on; 
for ani = 1:length(animallist)
   [xvalues2{ani},yvalues{ani},h(1,2+ani)] = cdfplot2(tuningidani{ani}) 
   if ismember(ani,applist)
       set( h(:,2+ani), 'LineStyle', '-', 'Color', [189,193,193]./255);% [0.6,0.6,0.9]
   else
       set( h(:,2+ani), 'LineStyle', '-', 'Color', [0.9,0.6,0.6]);
   end
end

maxxvalues = max(cellfun(@(x) max(x(2:end-1)),xvalues2));
minxvalues = min(cellfun(@(x) min(x(2:end-1)),xvalues2));
minbin = min(cellfun(@length,xvalues2));
maxbin = max(cellfun(@length,xvalues2));

for ani = 1:length(animallist)
    clear x indx
[x,indx] = unique((xvalues2{ani}));
interpx = minxvalues:((maxxvalues-minxvalues)/maxbin):maxxvalues-((maxxvalues-minxvalues)/maxbin);
interpcurve(ani,:) = interp1(x(2:end-1),yvalues{ani}(indx(2:end-1)),interpx);
end

rgb=[65 65 65]/256;
h(1,1) = plot(interpx,nanmedian(interpcurve(wtlist,:)),	'color',rgb,'LineWidth',2)% was 'b' for wt
h(1,2) = plot(interpx,nanmedian(interpcurve(applist,:)),'r','LineWidth',2)

set( h(:,1), 'LineStyle', '-', 'Color', 'k');
set( h(:,2), 'LineStyle', '-', 'Color', 'r');
xlabel('Area under the tuning curve', 'FontSize', 14);
ylabel('% of cells', 'FontSize', 14);
% xticks([0 0.2 0.4 0.6])
% xticklabels({'0','0.2','0.4','0.6'})
yticks([0 0.2 0.4 0.6 0.8 1])
yticklabels({'0','20','40','60','80','100'})

% legend([h(:,1),h(:,2),h(:,3),h(:,6)],{'Median wild type','Median AAP+','APP+','Wild type'},'Location','SouthEast'); legend boxoff
legend([h(:,1),h(:,2)],{'WT (6mo n=12) ','APP PS1+ (6mo n=12)'},'Location','SouthEast'); legend boxoff

% Check bandwith per bf tuning
toneseq = mat2cell(round(sort(toneorder)),1,ones(1,17));
colormapIndex = round(linspace(1,64,17));
C = colormap('jet');
figure(); hold on
for iii = 1:17
    [a(iii,:),b(iii,:)]=hist(Bandwidthneurons(neuronnbbf(find(indbfneurons==iii))),[1:17]);
    plot(b(iii,:),a(iii,:),'Color',C(colormapIndex(iii),:));
    meanbandwidth(iii) = mean(Bandwidthneurons(neuronnbbf(find(indbfneurons==iii))));
    stdbandwidth(iii) = std(Bandwidthneurons(neuronnbbf(find(indbfneurons==iii))));
end

xlabel('Bandwidth'); ylabel('Number of cells');
legend(cellfun(@num2str,toneseq,'UniformOutput',false)); legend boxoff;
title('Bandwidth distribution per frequency across animals')

figure(11); hold on; 
errorbar(cell2mat(toneseq),meanbandwidth,stdbandwidth,'.','MarkerSize',20);
xlabel('Tone frequency [Hz]'); ylabel('Bandwidth');

% Test bandwidth per frequency
[p,~,stats] = anova1(Bandwidthneurons(neuronnbbf),indbfneurons);
[a,b] = multcompare(stats);


figure(); hold on; 
for ani = 1:9
    [a,b]=hist(indbfneurons(intersect(neuronnbbf,find(indexneuronmouse==ani))),[1:17]);
    if ismember(ani,applist)
        plot(b,a,'Color',[1,0,0]);
    else
        plot(b,a,'Color',[0,0,1]);
    end
end
% xticks(cell2mat(toneseq));
xlabel('Tone Number');
ylabel('Number of cells');
legend({animallist}); legend boxoff;
title('Frequency distribution per animal')


% plot prct evoked response app vs wt
figure(); hold on; bar(1,mean(prctevoked(wtlist)),'FaceColor',[0.2,0.4,0.8]); plot(ones(1,length(wtlist)),prctevoked(wtlist),'.','MarkerSize',20,'Color',[0,0,0.8]);
bar(2,mean(prctevoked(applist)),'FaceColor',[1,0.2,0.2]); plot(ones(1,length(applist))*2,prctevoked(applist),'.','MarkerSize',20,'Color',[0.8,0,0]);
xticks([1,2])
xticklabels({'APP-','APP+'});
ylabel(' % of tone evoked neurons');

end
