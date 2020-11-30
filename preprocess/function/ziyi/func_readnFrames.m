function [nFrames, nFrames_oneplane] = func_readnFrames(mouse, varargin)

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('root', pwd)
p.parse(varargin{:});

filename = [p.Results.root '\' mouse '_config.csv'];
configTable = readtable(filename);

nFrames_oneplaneStr = configTable.nFrames_oneplane;
if ~iscell(nFrames_oneplaneStr); nFrames_oneplaneStr = {nFrames_oneplaneStr}; end
for i = 1:length(nFrames_oneplaneStr)
    nFrames_oneplaneSplit = strsplit(nFrames_oneplaneStr{i});
    for j = 1:2
        nFrames_oneplane(i,j) = str2double(nFrames_oneplaneSplit{j});
    end
end
nFrames = cumsum(nFrames_oneplane);

end