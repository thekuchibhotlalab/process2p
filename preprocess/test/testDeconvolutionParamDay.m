mice = {'cd017','cd036','cd037','cd041','cd042','cd044'};
multidayComp = {'ar1_allday', 'ar1_eachday'};

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

figure;
for i = 1:length(mice)
    subplot(2,3,i); hold on; 
    datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mice{i} '\'];
    comp = {'ar1_allday','ar1_eachday'};
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
    legend('fit all day','fit each day','Location','Best');
    title(mice{i}); xlabel('decay param'); ylabel('count')
end
suptitle('all day vs. each day fitting')
%%
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

