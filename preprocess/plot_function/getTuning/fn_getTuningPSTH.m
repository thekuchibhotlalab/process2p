function fn_getTuningPSTH(trialMean,toneMean,pretoneFrames,peakFrames,peakValue,toneLabel,toneindex,savePath)
%fn_getTuningPSTH - Description
%
% Syntax: fn_getTuningPSTH(input)
%
% Long description

nTones = size(trialMean,2); nFramesPerTone = size(trialMean,1);
%---------PLOT PSTH-----------
[nRow, nColumn] = fn_sqrtInt(ceil((nTones+1)/2));nColumn = nColumn*2;
psthFig = fn_figureSmartDim('widthHeightRatio',2.3,'hSize',0.7);
for i = 1:nTones
    subplot(nRow,nColumn,i);
    imagesc(squeeze(trialMean(:,i,:))');
    caxis([prctile(trialMean(:),5) prctile(trialMean(:),95)])
    title([toneLabel{toneindex(i)} 'HZ'])
    if i>(nRow-1)*nColumn
        xlabel('frames')
        frameAxis = pretoneFrames:10:nFramesPerTone;
        frameLabel = cellfun(@num2str,num2cell(frameAxis-pretoneFrames),'UniformOutput',false);
        xticks(frameAxis)
        xticklabels(frameLabel)
    else
        xticklabels([])
    end
    if mod(i,nColumn) == 1; ylabel('neurons');
    else; yticklabels([]); end
end
subplot(nRow,nColumn,nTones+1);imagesc(toneMean')
xticks(frameAxis);xticklabels(frameLabel);xlabel('frames')
caxis([prctile(trialMean(:),5) prctile(trialMean(:),95)])
yticklabels([]);title('Mean Tone Evoked')
saveas(gcf,[ savePath ...
        '/population/populationTonePSTH.png']);
saveas(gcf,[ savePath ...
    '/population/populationTonePSTH.fig']);

%---------PLOT PEAK FRAME-----------
latFig = figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.04, 0.45, 0.56]);% Enlarge figure to full screen.

subplot(2,2,1);histogram(peakFrames)
xlabel('peak frame');ylabel('frequency');title('peak frame of all neurons')

subplot(2,2,2); hold on; scatter(peakFrames,peakValue,15,'filled');  xlimm = xlim; 
plot(xlimm, [prctile(peakValue,50) prctile(peakValue,50)],'--','Color',[0.8 0.8 0.8]);
xlabel('peak frame'); ylabel('peak amplitude'); title('peak frame & amplitude')

subplot(2,2,3); [~,sortIndex] = sort(peakFrames); imagesc(toneMean(:,sortIndex)')
xticks(frameAxis); xticklabels(frameLabel); ylabel('roi'); title('average dff of all tones')
caxis([prctile(toneMean(:),5) prctile(toneMean(:),95)])

subplot(2,2,4); plot(toneMean,'Color',[0.8 0.8 0.8]);hold on; plot(mean(toneMean,2))
xticks(frameAxis); xticklabels(frameLabel); xlim([0 nFramesPerTone])
ylabel('dff');title('population average dff'); 
saveas(gcf,[ savePath ...
        '/population/populationTonePeak.png']); 
saveas(gcf,[ savePath ...
    '/population/populationTonePeak.fig']); 
end