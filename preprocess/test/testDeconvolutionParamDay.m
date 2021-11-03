mice = {'cd017','cd036','cd037','cd041','cd042','cd044'};
multidayComp = {'ar1_allday_fixBadSession', 'ar1_eachday'};
%%
figure;
for i = 1:length(mice)
    subplot(2,3,i); hold on; 
    datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mice{i} '\'];
    comp = {'ar1_allday','ar1_allday_smooth3'};
    plotPars = {};
    for j = 1:length(comp)
        filename = [datapath comp{j} '\' mice{i} '_calman_ar1_foo90_pars_allday_p.mat'];
        try 
            load(filename,'pars');
        catch
            filename = [datapath comp{j} '\' mice{i} '_calman_ar2_foo90_pars_allday_p.mat'];
            load(filename,'pars');
        end
        plotPars{j} = pars;
        histogram(plotPars{j},0:0.025:1); 
    end
    legend('no smooth','smooth','Location','Best');
    title(mice{i});xlabel('decay param'); ylabel('count')
end
suptitle('effect of smoothing on all day fitting')

figure;
for i = 1:length(mice)
    subplot(2,3,i); hold on; 
    datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mice{i} '\'];
    comp = {'ar1_eachday','ar1_eachday_smooth3'};
    plotPars = {};
    for j = 1:length(comp)
        filename = [datapath comp{j} '\' mice{i} '_calman_ar1_foo90_pars_allday_p.mat'];
        try 
            load(filename,'pars'); plotPars{j} = pars;
        catch
            filename = [datapath comp{j} '\' mice{i} '_calman_ar2_foo90_pars_day_p.mat'];
            load(filename,'pars');
            plotPars{j} = mean(cell2mat(pars),2);
        end
        histogram(plotPars{j},0:0.025:1); 
    end
    legend('no smooth','smooth','Location','Best');
    title(mice{i});xlabel('decay param'); ylabel('count')
end
suptitle('effect of smoothing on each day fitting')
%% compare all day and each day
figure;
for i = 1:length(mice)
    subplot(2,3,i); hold on; 
    datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mice{i} '\'];
    comp = {'ar1_allday','ar1_allday_fixBadSession'};
    plotPars = {};
    tempMean = [];
    for j = 1:length(comp)
        filename = [datapath comp{j} '\' mice{i} '_calman_ar1_foo90_pars_allday_p_plane0.mat'];
        try 
            load(filename,'pars'); plotPars{j} = pars;
        catch
            filename = [datapath comp{j} '\' mice{i} '_calman_ar1_foo90_pars_allday_p.mat'];
            load(filename,'pars');plotPars{j} = pars;
            %plotPars{j} = mean(cell2mat(pars),2);
        end
        histogram(plotPars{j},0:0.025:1); 
        temp = plotPars{j}; temp(temp==0) = [];
        tempMean(j) = mean(temp,1);
    end
    legend('fit all day','fix bad session','Location','Best');
    title([mice{i} ' mean=' num2str(tempMean(1),'%.2f') ' ' num2str(tempMean(2),'%.2f')]); xlabel('decay param'); ylabel('count')
end
suptitle('all day vs. each day fitting')

figure;
for i = 1:length(mice)
    subplot(2,3,i); hold on; 
    datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mice{i} '\'];
    comp = {'ar1_allday','ar1_allday_fixBadSession'};
    plotPars = {};
    tempMean = [];
    for j = 1:length(comp)
        filename = [datapath comp{j} '\' mice{i} '_calman_ar1_foo90_pars_allday_p_plane0.mat'];
        try 
            load(filename,'pars'); plotPars{j} = pars;
        catch
            filename = [datapath comp{j} '\' mice{i} '_calman_ar1_foo90_pars_allday_p.mat'];
            load(filename,'pars');plotPars{j} = pars;
            %plotPars{j} = mean(cell2mat(pars),2);
        end
        histogram(log(0.5)./log(plotPars{j})* (1/15.63),0:0.025:2.0); 
        temp = logical(plotPars{j}~=0);
        tempMean(j) = mean(log(0.5)./log(plotPars{j}(temp))* (1/15.63),1);
    end
    legend('fit all day','fix bad session','Location','Best');
    title([mice{i} ' mean=' num2str(tempMean(1),'%.2f') ' ' num2str(tempMean(2),'%.2f')]); xlabel('decay time (sec)'); ylabel('count')
end
suptitle('all day vs. each day fitting')
%% compare TC
methodCorr = {};
for i = 1:length(mice)
   
    datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mice{i} '\'];
    comp = {'ar1_allday','ar1_eachday'};
    allS = {};
    for j = 1:length(comp)
        filename = [datapath comp{j} '\' mice{i} '_calman_ar1_foo90_pars_allday_s.mat'];
        try 
            load(filename,'s'); allS{j} = s;
        catch
            filename = [datapath comp{j} '\' mice{i} '_calman_ar2_foo90_pars_day_s.mat'];
            load(filename,'s');
            allS{j} = s;
        end
        
        disp(1);
    end
    tempCorr = [];
    for k = 1:length(allS{1})
        allDayS = allS{1}{k};eachDayS = allS{2}{k};
        for l = 1:size(allDayS,1)
            corrValue = corr(allDayS(l,:)',eachDayS(l,:)');
            tempCorr(k,l) = corrValue;
        end
    end
    methodCorr{i} = tempCorr;
end
%% Plot correlation of TC
figure;
for i = 1:length(mice)
    subplot(2,3,i); hold on; 
    histogram(mean(methodCorr{i},1),0.3:0.01:1)
    
    title(mice{i}); xlabel('spike corr allday-eachday'); ylabel('cell count')
end
suptitle('all day vs. each day fitting')

%% stability 
figure;
for i = 1:length(mice)
    subplot_tight(1,6,i);
    datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mice{i} '\'];
    comp = 'ar1_eachday';
    plotPars = {};
    
    filename = [datapath comp '\' mice{i} '_calman_ar1_foo90_pars_allday_p.mat'];
    try 
        load(filename,'pars');
    catch
        filename = [datapath comp '\' mice{i} '_calman_ar2_foo90_pars_day_p.mat'];
        load(filename,'pars');
    end
    pars = cell2mat(pars);
    imagesc(pars); caxis([0.6 1]);
    title(mice{i}); xlabel('days')
end
colorbar;
suptitle('decay parameter across days')




