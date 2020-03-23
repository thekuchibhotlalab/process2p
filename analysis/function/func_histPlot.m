function f1 = func_histPlot (data, legends, titles, xlimm)
    %myFun - Description
    %
    % Syntax: f1 = func_piePlot (data, legends, title)
    %
    % Long description
        
    nPlot = length(data);
    f1 = figure; 
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.025, 0.3, 0.95, 0.4]);
    
    for i = 1:nPlot
        subplot(1,nPlot,i)
        temp = data{i};
        for j = 1:length (temp)
            hold on;
            histogram(temp{j},'FaceAlpha',0.3);
        end 
        legend(legends{i},'Location','Best'); title (titles{i})
        if ~isempty (xlimm)
           xlim(xlimm{i});
        end
    end
    
    end