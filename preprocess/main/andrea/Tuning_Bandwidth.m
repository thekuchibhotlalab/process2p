function Tuning_Bandwidth(User)
% Input example
% User = 'kelly';
% Animal = 'kf024';
% FileNb = '000_002';
% Planenb = '1';
% a = genpath('Y:\Shared\Matlab\preprocessing\');
% addpath(a)

% careful the planes are not taken into account !
[~,~,table] = xlsread(['T:\LabData2\' User '\Behavior_Imaging\Tuning_curve.xlsx']);
animallist = vertcat(table{2:end,1});
genolist = {table{2:end,4}};
applist = find(strcmp(genolist,'APP+'));
wtlist = find(strcmp(genolist,'WT'));
datalist = vertcat(table{2:end,2});
filelist = {table{2:end,3}};

for ani = 1:length(animallist)
    Animal = animallist(ani,:);
    
    clearvars -except RootAddress table applist wtlist filelist datalist Animal ani animallist User numroiperbin indcorrtun indcorrtunstd indcorrtunshuffletuning indcorrtunshuffledist indcorrtunshuffletuning xvalues2 alldist allcorrtun
    % Load necessary data
%     if datalist(ani) == 2
%         RootAddress = ['T:\LabData' num2str(datalist(ani)) '\' User '\Behavior_Imaging\' animallist(ani,1:end) '\preprocessData\'];
%     elseif datalist(ani) == 3
%         RootAddress = ['Z:\LabData' num2str(datalist(ani)) '\' User '\' animallist(ani,1:end) '\preprocessData\'];
%     end
%         RootAddress = ['X:\LabData' num2str(datalist(ani)) '\' User '\' Animal '\preprocessData\'];
%     
%     if ani == 4
%         load('X:\LabData2\kelly\kf026\preprocessData\TuningCurve_kf026_000_001\population\tuningMean.mat')
%         load('X:\LabData2\kelly\kf026\preprocessData\TuningCurve_kf026_000_001\population\tuningMedian.mat')
%         load('X:\LabData2\kelly\kf026\preprocessData\TuningCurve_kf026_000_001\population\responsiveCellFlag.mat')
%         
%     elseif ani == 5
%         load('X:\LabData2\kelly\kf028\preprocessData\TuningCurve_kf028_000_000\population\tuningMean.mat')
%         load('X:\LabData2\kelly\kf028\preprocessData\TuningCurve_kf028_000_000\population\tuningMedian.mat')
%         load('X:\LabData2\kelly\kf028\preprocessData\TuningCurve_kf028_000_000\population\responsiveCellFlag.mat')
%         
%     elseif ani == 6
%         load('X:\LabData2\kelly\kf007\preprocessData\TuningCurve_kf007_003_002\population\tuningMean.mat')
%         load('X:\LabData2\kelly\kf007\preprocessData\TuningCurve_kf007_003_002\population\tuningMedian.mat')
%         load('X:\LabData2\kelly\kf007\preprocessData\TuningCurve_kf007_003_002\population\responsiveCellFlag.mat')
%     else
%         load([RootAddress 'TuningCurve_'  filelist{ani} '\population\tuningMean.mat'])
%         load([RootAddress 'TuningCurve_'  filelist{ani} '\population\responsiveCellFlag.mat'])
%         load([RootAddress 'TuningCurve_'  filelist{ani} '\population\neuronRespTable.mat'])
         tuning=load([RootAddress 'TuningCurve_'  filelist{ani} '\population\tuning.mat']);
         responsiveCellFlag = tuning.anovaPeakCorr(2,:);
         tuningMean=tuning.popTuningMean;
%          tuningMedian=tuning.poptuningMedian;
%     end
    load([RootAddress(1:end-15)  filelist{ani} '.mat']);
    
    % Get recording parameters such as magnification ect
    MagnificationValue = info.config.magnification_list(info.config.magnification); % apparently magnification of 5 in info config corresponds to 1
    Size_Image = 1000/str2num(MagnificationValue);
    Pixel_Micron = Size_Image/513; % from the size of the roi mean image
    % Get Tuned Neurons
    TunedNeuronstmp = find(responsiveCellFlag);%[freqind,TunedNeuronstmp] = find((neuronRespTable(1:end-1,2:end))>0); % see Ziyi's ttest% more conservative mutlcompare unlike ResponsiveCellFlag
    TunedNeurons = unique(TunedNeuronstmp);
    
    if isempty(info.otparam) % JL 04/30/19
        nPlanes = 1;%nPlanes = 4;
        NeuronsPlanesInd = [];
        roiName = [RootAddress 'ROI\' filelist{ani} '_roi' num2str(nPlanes) '_green_aligned_on_green.zip']; %JL 26/04/19 change name green on green
        rois = ReadImageJROI(roiName);
        NeuronsPlanesInd = [NeuronsPlanesInd ones(1,length(ReadImageJROI(roiName)))*1];
    else
        nPlanes = length(info.otwave_um);%nPlanes = 4;
        NeuronsPlanesInd = [];
        rois = [];
        for i = 1:nPlanes
            roiName = [RootAddress 'ROI\' filelist{ani} '_roi' num2str(i) '_green_aligned_on_green.zip']; %JL 26/04/19 change name green on green
            rois = [rois ReadImageJROI(roiName)];
            NeuronsPlanesInd = [NeuronsPlanesInd ones(1,length(ReadImageJROI(roiName)))*i];
        end
    end
    
    
    
    for i = 1:size(rois,2)
        % Significantones{TunedNeurons(i)} = neuronRespTable(freqind(TunedNeuronstmp==TunedNeurons(i)),1);
        xy_median{i} = [median(rois{i}.mnCoordinates(:,1)),median(rois{i}.mnCoordinates(:,2))];
    end
    xvalues = cellfun(@(x) x(1),xy_median);
    yvalues = cellfun(@(x) x(2),xy_median);
    % for kf028 remove ROI with x > 450
    % if strcmp(Animal,'kf028')
    % TunedNeurons(ismember(TunedNeurons,find(xvalues-450>0))) = [];
    % end
    [~,bffreqind] = max(tuningMean(:,TunedNeurons));
    BF = neuronRespTable(bffreqind,1);
    % Bandwidth_Large = [cellfun(@(x) log2(max(x)/min(x)),Significantones(TunedNeurons))];%;cellfun(@max,Significantones(TunedNeurons))]%sum(neuronRespTable(1:end-1,TunedNeurons+1));% refine this, how about nb of octive it's responsive to ?
    % Qfactor = BF'./Bandwidth_Large; % BF/Bandwidth according to Yue 2014 %inf is no spread since log2(1) = 0
    cellbf = mat2cell(BF,ones(1,size(BF,1)),1);
    
    for p = 1:nPlanes
        clear Distance_pairwise initialmat Tuning_corr shuffletmp shufflemat Shuffle_corr  alldist{ani} allcorrtun{ani} allcorrshuffle startbin indcorrtunshuffletuningtest randomperm indcorrtunshuffledisttest
        listneuronplane{p} = intersect(find(NeuronsPlanesInd == p),TunedNeurons);
        for i = 1:size(listneuronplane{p},2) % this only takes into acount 1 plane
            Distance_pairwise(i,:) = cellfun(@(x,y) sqrt((x(1)-y(1))^2+(x(2)-y(2))^2),xy_median(listneuronplane{p}),repmat(xy_median(listneuronplane{p}(i)),1,size(listneuronplane{p},2)));
            [RelativeDist{p}{i},NeuronInd{i}] = sort(Distance_pairwise(i,1:find(Distance_pairwise(i,:)==0)));
        end
        
        initialmat = tuningMean(:,listneuronplane{p})-1;
        Tuning_corr = corrcoef(initialmat);%./max(tuningMean(:,TunedNeurons)));
        shuffletmp = cellfun(@(x) datasample([1:17],17,'Replace',false)',mat2cell([1:size(listneuronplane{p},2)],1,ones(1,size(listneuronplane{p},2))),'UniformOutput',false);
        shufflemat = initialmat(horzcat(shuffletmp{:}));
        Shuffle_corr = corrcoef(shufflemat);
        % Tuning_corr2 = corrcoef(tuningMean(:,TunedNeurons)-mean(tuningMean(:,TunedNeurons)));
        
        for i = 1:size(listneuronplane{p},2)
            Tuning_corr_sorted{p}{i} = Tuning_corr(i,NeuronInd{i});
            Tuning_corr_sorted_suffle{p}{i} = Shuffle_corr(i,NeuronInd{i});
            %     Tuning_corr_sorted_2{i} = Tuning_corr2(i,NeuronInd{i});
        end
    end
        alldist{ani} = cell2mat([RelativeDist{:}]);
        allcorrtun{ani} = cell2mat([Tuning_corr_sorted{:}]);%horzcat(Tuning_corr_sorted{:});
        % allcorrtun{ani}2 = horzcat(Tuning_corr_sorted{:});
        allcorrshuffle = cell2mat([Tuning_corr_sorted_suffle{:}]);%horzcat(Tuning_corr_sorted_suffle{:});
        %     figure();plot(alldist{ani},allcorrtun{ani},'.');
    
        distancebin = 20;
        [startbin] = prctile(alldist{ani},[0:100/distancebin:100]);%[bindist,startbin]=discretize(alldist{ani},50,'IncludedEdge','right');% prctile(X,[25 50 75],1)
        cellstart = mat2cell(startbin,1,ones(1,size(startbin,2)));
        % [binvalue,groupbin] = sort(bindist);
        % figure();plot(binvalue,allcorrtun{ani}(groupbin),'.');
        numroiperbin{ani} = cellfun(@(x,y,d) numel(find(x>(y) & x<=d)),repmat({alldist{ani}},1,distancebin),cellstart(1:end-1),cellstart(2:end));
        indcorrtun{ani} = cellfun(@(x,y,d,z) median(z(find(x>(y) & x<=d))),repmat({alldist{ani}},1,distancebin),cellstart(1:end-1),cellstart(2:end),repmat({allcorrtun{ani}},1,distancebin));
        indcorrtunstd{ani} = cellfun(@(x,y,d,z) std(z(find(x>(y) & x<=d))),repmat({alldist{ani}},1,distancebin),cellstart(1:end-1),cellstart(2:end),repmat({allcorrtun{ani}},1,distancebin));
        
        % indcorrtun2 = cellfun(@(x,y,d,z) median(z(find(x>(y) & x<=d))),repmat({alldist{ani}},1,40),cellstart(1:end-1),cellstart(2:end),repmat({allcorrtun{ani}2},1,40));
        for nbperm = 1:100
            randomperm(nbperm,:) = (randperm(size(alldist{ani},2)));
            indcorrtunshuffledist{ani}(nbperm,:) = cellfun(@(x,y,d,z) median(z(find(x>(y) & x<=d))),repmat({alldist{ani}(randomperm(nbperm,:))},1,distancebin),cellstart(1:end-1),cellstart(2:end),repmat({allcorrtun{ani}},1,distancebin));
            indcorrtunshuffledisttest(nbperm,:) = cellfun(@(x,y,d,z) (z(find(x>(y) & x<=d))),repmat({alldist{ani}(randomperm(nbperm,:))},1,distancebin),cellstart(1:end-1),cellstart(2:end),repmat({allcorrtun{ani}},1,distancebin),'UniformOutput',0);
        end
        indcorrtunshuffletuning{ani} = cellfun(@(x,y,d,z) median(z(find(x>(y) & x<=d))),repmat({alldist{ani}},1,distancebin),cellstart(1:end-1),cellstart(2:end),repmat({allcorrshuffle},1,distancebin));
        
        xvalues2{ani} = cumsum(diff(startbin))*Pixel_Micron;
        % statistical test, test distribution of distqnce bins agains random
        % version
        % indcorrshuffletuningtest = cellfun(@(x,y,d,z) (z(find(x>(y) & x<=d))),repmat({alldist{ani}},1,distancebin),cellstart(1:end-1),cellstart(2:end),repmat({allcorrshuffle},1,distancebin),'UniformOutput',0);
        % indcorrtuntest = cellfun(@(x,y,d,z) (z(find(x>(y) & x<=d))),repmat({alldist{ani}},1,distancebin),cellstart(1:end-1),cellstart(2:end),repmat({allcorrtun{ani}},1,distancebin),'UniformOutput',0);
        % for i = 1:distancebin
        %     pvaldist(i) = ranksum(cell2mat(indcorrtunshuffledisttest(:,i)'),indcorrtuntest{i});
        %     pvaltuning(i) = ranksum(indcorrshuffletuningtest{i},indcorrtuntest{i});
        % end
        
        %     figure(2); hold on;
        %     if p == 1 && nPlanes<4
        %         h =openfig([RootAddress 'TuningCurve_' Animal '_' FileNb '\population\tuningMap.fig'])
        %         x0=10;
        %         y0=10;
        %         width=900;
        %         height=1000;
        %         set(gcf,'position',[x0,y0,width,height])
        %         hold on;
        %     elseif p == 1 && nPlanes ==4
        %         h = figure();
        %         x0=10;
        %         y0=10;
        %         width=900;
        %         height=1000;
        %         set(gcf,'position',[x0,y0,width,height])
        %         hold on;
        %     end
        %
        %     if nPlanes>1 && nPlanes<4
        %
        %         subplot(2,nPlanes,p+nPlanes); hold on;
        %         A = shadedErrorBar(cumsum(diff(startbin))*Pixel_Micron,indcorrtun,indcorrtunstd,{'Color',[0,0,0]},1); % change the xvalue to center of bins plot(cumsum(diff(startbin)),indcorrtun2,'r');
        %         B = plot(cumsum(diff(startbin))*Pixel_Micron,mean(indcorrtunshuffledist),'Color',[0.2,0.2,0.8]);%shadedErrorBar(cumsum(diff(startbin)),mean(indcorrtunshuffledist),std(indcorrtunshuffledist),{'Color',[0.1,0.1,0.1]});
        %         C = plot(cumsum(diff(startbin))*Pixel_Micron,indcorrtunshuffletuning,'Color',[0.8,0.2,0.2]);
        %         axis tight; xlabel(['Intracellular distance [µm] cells=' num2str(length(listneuronplane{p}))]); ylabel('Median tuning Correlation');
        %         [hh, icons] = legend([A.mainLine B C],{'Median tuning corr.','Shuffle distance','Shuffle tuning'}); legend boxoff
        %         %         p1 = icons(1).Position;
        %         %         icons(1).Position = [0.3 p1(2) 0];
        %         %         icons(2).XData = [0.05 0.1];
        %     elseif nPlanes == 4
        %         subplot(2,2,p); hold on; title(['plane ' num2str(p)]);
        %         A = shadedErrorBar(cumsum(diff(startbin))*Pixel_Micron,indcorrtun,indcorrtunstd,{'Color',[0,0,0]},1); % change the xvalue to center of bins plot(cumsum(diff(startbin)),indcorrtun2,'r');
        %         B = plot(cumsum(diff(startbin))*Pixel_Micron,mean(indcorrtunshuffledist),'Color',[0.2,0.2,0.8]);%shadedErrorBar(cumsum(diff(startbin)),mean(indcorrtunshuffledist),std(indcorrtunshuffledist),{'Color',[0.1,0.1,0.1]});
        %         C = plot(cumsum(diff(startbin))*Pixel_Micron,indcorrtunshuffletuning,'Color',[0.8,0.2,0.2]);
        %         axis tight; xlabel(['Intracellular distance [µm] cells=' num2str(length(listneuronplane{p}))]); ylabel('Median tuning Correlation');
        %         [hh, icons] = legend([A.mainLine B C],{'Median tuning corr.','Shuffle distance','Shuffle tuning'}); legend boxoff
        %         %         p1 = icons(1).Position;
        %         %         icons(1).Position = [0.3 p1(2) 0];
        %         %         icons(2).XData = [0.05 0.1];
        %     elseif nPlanes == 1
        %         subplot(2,2,3); hold on;
        %         A = shadedErrorBar(cumsum(diff(startbin))*Pixel_Micron,indcorrtun,indcorrtunstd,{'Color',[0,0,0]},1); % change the xvalue to center of bins plot(cumsum(diff(startbin)),indcorrtun2,'r');
        %         B = plot(cumsum(diff(startbin))*Pixel_Micron,mean(indcorrtunshuffledist),'Color',[0.2,0.2,0.8]);%shadedErrorBar(cumsum(diff(startbin)),mean(indcorrtunshuffledist),std(indcorrtunshuffledist),{'Color',[0.1,0.1,0.1]});
        %         C = plot(cumsum(diff(startbin))*Pixel_Micron,indcorrtunshuffletuning,'Color',[0.8,0.2,0.2])
        %         axis tight; xlabel(['Intracellular distance [µm] cells=' num2str(length(listneuronplane{p}))]); ylabel('Median tuning Correlation');
        %         legend([A.mainLine B C],{'Median tuning corr.','Shuffle distance','Shuffle tuning'}); legend boxoff
        %     end
    
    
end
animalwithenoughcells = find(cellfun(@(x) x(end),numroiperbin)>100);

% plot correlation
% for ani = animalwithenoughcells
% figure(1); hold on; plot(xvalues2{ani},indcorrtun{ani}-mean(indcorrtunshuffledist{ani}));
% figure(2); hold on; plot(xvalues2{ani},indcorrtun{ani}-indcorrtunshuffletuning{ani});
% %  plot(xvalues2{ani},mean(indcorrtunshuffledist{ani}));
% end
appxvalues = [xvalues2{1} xvalues2{3}];
appindcorr = [indcorrtun{1}-(indcorrtunshuffletuning{1}) indcorrtun{3}-(indcorrtunshuffletuning{3})]
wtxvalues = [xvalues2{4} xvalues2{5}];
wtindcorr = [indcorrtun{4}-(indcorrtunshuffletuning{4}) indcorrtun{5}-(indcorrtunshuffletuning{5})]
[a,lim] = hist([xvalues2{animalwithenoughcells}],30);

for bin = 2:length(lim)
appval{bin}= appindcorr(appxvalues>lim(bin-1) & appxvalues<=lim(bin));
wtval{bin}= appindcorr(wtxvalues>lim(bin-1) & wtxvalues<=lim(bin));
end
appval{1} = appindcorr(appxvalues<lim(1));
wtval{1}= wtindcorr(wtxvalues<lim(1));


figure(1); hold on; 
A = shadedErrorBar(lim(1,1:9),cellfun(@median,appval(1:9)),cellfun(@std,appval(1:9)),{'Color',[1,0,0]},1);
B = shadedErrorBar(lim(1,1:9),cellfun(@median,wtval(1:9)),cellfun(@std,wtval(1:9)),{'Color',[0,0,1]},1);
xticks(round(lim(1,1:9)))
xlim([lim(1) lim(9)])
xlabel(['Intracellular distance [µm] (bins)']);ylabel('Correlation');
legend([A.mainLine B.mainLine],{'APP+','WT'}); legend boxoff

%individual curves
figure(2); hold on; 
A = plot(xvalues2{1},indcorrtun{1}-indcorrtunshuffletuning{1},'Color',[1,0,0]);
A1 = plot(xvalues2{3},indcorrtun{3}-indcorrtunshuffletuning{3},'Color',[1,0,0]);
B = plot(xvalues2{4},indcorrtun{4}-indcorrtunshuffletuning{4},'Color',[0,0,1]);
B1 = plot(xvalues2{5},indcorrtun{5}-indcorrtunshuffletuning{5},'Color',[0,0,1]);


C = plot(xvalues2{1},mean(indcorrtunshuffledist{1}),'-k')
plot(xvalues2{1},mean(indcorrtunshuffledist{3}),'-k')
plot(xvalues2{1},mean(indcorrtunshuffledist{4}),'-k')
plot(xvalues2{1},mean(indcorrtunshuffledist{5}),'-k')
legend([A A1 B B1],{'APP+: kf022','APP+: kf024','WT: kf026','WT: kf028'}); legend boxoff
xlabel(['Intracellular distance [µm]']);ylabel('Correlation corrected for tuning');
end