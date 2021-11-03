function Plot_tone_evoked()
clear all

load('w:\LabData3\jenni\cd014\preprocessData\TuningCurve_cd014_001_000\population\tuningMean.mat')
load('w:\LabData3\jenni\cd014\preprocessData\TuningCurve_cd014_001_000\population\tuningMedian.mat')
load('w:\LabData3\jenni\cd014\preprocessData\TuningCurve_cd014_001_000\population\TCtoneMean.mat')

NbFrame = 100; 
Neuron = 46;% 46
framelim = [0:8]*NbFrame;
MeanResponseTone = [nan(1,50) TCtoneMean(1:NbFrame,1,Neuron)'];
for tone = 2:17
MeanResponseTone = [MeanResponseTone nan(1,50) TCtoneMean(NbFrame-3:NbFrame,tone-1,Neuron)' TCtoneMean(1:NbFrame,tone,Neuron)' nan(1,50)];
end
% end
% MeanResponseTone = cell2mat(cellfun(@(x) mean(x,1),ResponseTone,'UniformOutput',false)');
% hFig = figure(); hold on; 
% for tone = 1:17
%    plot(TCtoneMean(1:NbFrame,tone,306),'b') 
% end
 colormapIndex = round(linspace(40,250,17));
 C = colormap('jet');
xvalues(1,:) = [1,150];
for tone = 2:17
xvalues(tone,:) = [xvalues(tone-1,end)+1,xvalues(tone-1,end)+200];
end



toneorder = [45254.834 8000 13454.34264 4756.82846 5656.854249,...
        22627.417 64000 53817.37058 4000 9513.65692,...
        16000 6727.171322 19027.31384 26908.6852 32000,...
        11313.7085 38054.62768];

hFig = figure(); hold on;
for tone = 1:17
plot([xvalues(tone,1):xvalues(tone,2)],MeanResponseTone(1,xvalues(tone,1):xvalues(tone,2)),'Color',C(colormapIndex(tone),:));
end
color = get(hFig,'Color');
set(gca,'XColor',color,'YColor',color,'TickDir','out')
title(['CD014 - ROI:' num2str(Neuron)])

legend({'4000','4757','5657','6727','8000','9516','11314','16000','19027','22627','26909','32000','38055','45255','53817','64000'},'FontSize',6)
figure(); hold on; 
for tone = 1:17
plot((1:100)*1/30,smooth(TCtoneMean(:,tone,Neuron)-1,3),'Color',C(colormapIndex(tone),:));
end
ylabel('Fluorescence \DeltaF/F_0')
xlabel('Time [s]')
xlim([1*1/30 100*1/30])
ylim([0 0.3])
title(['CD014 - ROI:' num2str(Neuron)])

legend({'4000','4757','5657','6727','8000','9516','11314','16000','19027','22627','26909','32000','38055','45255','53817','64000'},'FontSize',6)
legend boxoff 
