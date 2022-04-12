function fn_getTuningSingleNeuron(popTuning,tuning,TC_original,TCpretone_reorder,...
    nFramesPerTrial,nFramesPerTone,nTrials,trialMean,trialMedian,toneAct,baseAct_noAvg,...
    toneLabel,toneindex, nNeuron,nTones,roisBound,refImg,neuronEachPlane,savePath)
%fn_getTuningSingleNeuron - Description
%
% Syntax: fn_getTuningSingleNeuron()
%
% Long description
global nPlanes
trialMeanTrialTC = reshape(trialMean,[nFramesPerTrial,nNeuron]);
trialMedianTrialTC = reshape(trialMedian,[nFramesPerTrial,nNeuron]);
neuronPlane = cumsum(neuronEachPlane);
try  neuronPlane = [0;neuronPlane];
catch;  neuronPlane = [0 neuronPlane];disp('check this'); end
C = colormap('jet');colormapIndex = round(linspace(1,size(C,1),nTones));

for j = 1:nPlanes
    for i = 1:neuronEachPlane(j)
        tic;
        tuningFig = figure('visible','off');
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.04, 0.6, 0.9]);% Enlarge figure to full screen.
        cellIndex = neuronPlane(j) + i;
        % REMEMBER TO CHANGE THIS LINE
        
        subplot(2,2,1)
        imagesc(refImg{j});colormap gray;hold on;
        if ~isempty(roisBound{j}{i})
            x=roisBound{j}{i}(:,1); %freehand rois have the outlines in x-y coordinates
            y=roisBound{j}{i}(:,2); %matlab matrices are inverted so x values are the 2nd column of mn coordinates, and y is the 1st columna
            patch(x,y,C(colormapIndex(popTuning.bfMedian(cellIndex)),:),'EdgeColor','none');
            title(['Cell #' int2str(cellIndex) ' Plane #' int2str(j) ' FoV'])
        end
        subplot(4,2,2)
        plot(TC_original(:,cellIndex),'LineWidth',0.5);
        xlim([1 size(TC_original,1)]);
        xlabel('Frames')
        ylabel('Fluorescence')
        title('Raw Fluorescence')

        subplot(4,2,4)
        tempSampleRate = 1;
        timeAxis = (0:tempSampleRate:nFramesPerTrial-1) * nPlanes / 30;
        TC_trial = reshape(TCpretone_reorder,[nFramesPerTrial,nTrials,nNeuron]);
        for k = 1:nTrials     
            plot(timeAxis,TC_trial(1:tempSampleRate:nFramesPerTrial,k,cellIndex),'color',[0.9 0.9 0.9],'LineWidth',0.5); hold on
        end

        timeAxis = (0:nFramesPerTrial-1) * nPlanes / 30;
        
        tempNTones = 8; colormapIndex = round(linspace(1,64,17)); 
        jjet = jet; toneColor = jjet(colormapIndex,:);
        for k = 1:tempNTones
            tempIdx =((k-1)*50+1):(k*50);
            p1 = plot(timeAxis(tempIdx), trialMedianTrialTC(tempIdx,cellIndex),'LineWidth',3,'color',toneColor(k*2,:));
        end
        %p1 = plot(timeAxis, trialMedianTrialTC(:,cellIndex),'LineWidth',1.2,'color',[0.0000 0.4470 0.7410]);
        %p2 = plot(timeAxis, trialMeanTrialTC(:,cellIndex),'LineWidth',1.2,'color',[0.8500 0.3250 0.0980]);
        yMax = 0.7 * max(trialMedianTrialTC(:,cellIndex)) + 0.3 * max(max(TC_trial(1:tempSampleRate:nFramesPerTrial,:,cellIndex)));
        yMin = 0.7 * min(trialMeanTrialTC(:,cellIndex)) + 0.3 * min(min(TC_trial(1:tempSampleRate:nFramesPerTrial,:,cellIndex)));
        onsetTime = (1:nFramesPerTone:nFramesPerTrial) * nPlanes / 30;
        for k = 1:nTones
            plot([onsetTime(k) onsetTime(k)],[yMin yMax],'LineWidth',0.5,'color',[0.4 0.4 0.4]);
        end
        nTicks = 5; tickLocation = round(linspace(1,nTones,nTicks));temp = onsetTime+(onsetTime(2)-onsetTime(1))/2;
        xticks(temp(tickLocation));xticklabels(toneLabel(toneindex(tickLocation)));xlim([0 timeAxis(end)]);xlabel('Frequency(kHz)')
        ylim([yMin yMax]);ylabel('dF/F');title('Mean Activity Across Trials')

        subplot(2,2,3)
        [rocMax,rocMaxIdx] = max(tuning.roc.rocAuc(:,cellIndex));
        %toneAct = squeeze(TCpretone_reorder(pretoneFrames+peakFrames(cellIndex),:,:,cellIndex));
        %baseAct_noAvg = squeeze(TCpretone_reorder(pretoneFrames-4:pretoneFrames,:,:,cellIndex));
        tempBaseAct = squeeze(baseAct_noAvg(:,rocMaxIdx,:,cellIndex));
        rocAct = [tempBaseAct(:)' squeeze(toneAct(rocMaxIdx,:,cellIndex))];
        rocLabel = [zeros(1,numel(tempBaseAct)) ones(1,size(toneAct,2))];
        [tpr, fpr, threshold] = roc(rocLabel, rocAct);
        h_auc = plot([0 fpr 1],[0 tpr 1],'Color', [0 0 0], 'LineWidth', 1.2);hold on; 
        plot([0 1],[0 1],'Color',[0.8,0.8,0.8],'LineWidth',0.5)
        xlabel('false positive')
        ylabel('true positive')
        xlim([0 1])
        ylim([0 1])
        title(['ROC ' toneLabel{toneindex(rocMaxIdx)} 'HZ'])
        legend(h_auc,['AUC ' num2str(rocMax,'%0.2f')])

        signifTonePoint = tuning.responsiveCellToneFlag(:,cellIndex);
        signifTuningMedian = tuning.median(logical(signifTonePoint),cellIndex);
        signifTonePoint = find(logical(signifTonePoint));
        %neuronRespTable(1:end-1,[1 cellIndex+1]);
        %signifTuningMedian = trialMedian(:,cellIndex);
        %signifTuningMedian = signifTuningMedian(signifTonePoint(:,2)==1);
        %signifTonePoint = signifTonePoint(signifTonePoint(:,2)==1,1);
        
        nTicks = 5; tickLocation = round(linspace(1,nTones,nTicks));
        
        subplot(2,2,4)
        h_med = errorbar(1:nTones, tuning.median(:,cellIndex),...
            tuning.medianSEM(:,cellIndex),...
            'LineWidth',2,'color',[0.0000 0.4470 0.7410]); hold on;
        h_mean = errorbar(1:nTones, tuning.mean(:,cellIndex),...
            tuning.meanSEM(:,cellIndex),...
            'LineWidth',2,'color',[0.8500 0.3250 0.0980]);
        scatter(signifTonePoint,signifTuningMedian,40,[0.2 0.2 0.2],'*');
        
        set(gca, 'XTick', tickLocation); set(gca, 'XTickLabel', toneLabel(toneindex(tickLocation)));
        xlim([1 nTones]); xlabel('Frequency (kHz)');ylabel('dF/F');

        yyaxis right
        h_roc = plot(1:nTones,tuning.roc.rocAuc(:,cellIndex),'Color',[0.8 0.8 0.8],'LineWidth',2);
        ylabel('AUC')        
        title('Tuning Curve');

        if tuning.responsiveCellFlag(cellIndex)
            ylimit = get(gca,'YLim');
            scatter(2,ylimit(2)*0.9+ylimit(1)*0.1,40,[0 0 0],'*')
        end
        legend([h_med, h_mean, h_roc],'Median','Mean','ROC');

        saveas(tuningFig,[savePath ...
            '/singleNeuron/Neuron' num2str(cellIndex,'%03d') '.png']);
        saveas(tuningFig,[savePath ...
            '/singleNeuron/Neuron' num2str(cellIndex,'%03d') '.fig']);
        
        close(tuningFig);
        timeElapsed = toc;

        disp(['Cell #' int2str(cellIndex) ' Plane #' int2str(j) ' Time=' num2str(timeElapsed,'%03f') ' secs'])

    end
end

end