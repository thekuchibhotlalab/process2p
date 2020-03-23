function f1 = func_piePlot (data, legends, titles)
%myFun - Description
%
% Syntax: f1 = func_piePlot (data, legends, title)
%
% Long description
    
nPlot = length(data);
f= figure; 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.04, 0.6, 0.9]);

for i = 1:nPlot
    subplot(1,nPlot,i)
    temp = data{i};
    newData = [];
    for j = 1:length(temp)
        newData = [newData sum(temp{j})];
    end
    newData(newData==0) = 0.1;
    pie(newData / sum(newData),ones(1,length(newData)));
    legend(legends{i},'Location','Best'); title (titles{i})

end

end