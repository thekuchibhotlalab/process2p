function [TC,neuronEachPlane] = func_loadTCmanual(varargin)
%func_loadTCmanual - Description
%
% Syntax: [TC,neuronEachPlane] = func_loadTCmanual(varargin)
%
% Long description
    
p = func_createInputParser();
p.addParameter('file', [])
p.parse(varargin{:});
sep = '\';
%---------GET NUMBER OF CHANNELS-----------
nPlanes = str2double(p.Results.nPlanes);
if iscell(p.Results.functionalChannel)
    nFuncChannel = length(p.Results.functionalChannel);
    functionalChannel = p.Results.functionalChannel;
else
    nFuncChannel = 1;
    functionalChannel = {p.Results.functionalChannel};
end
nFrames_oneplane = p.Results.nFrames_oneplane;



tcFileSplit = strsplit(p.Results.tcFile);
tcFile = reshape(tcFileSplit,nFuncChannel,nPlanes);
TC = cell(1,nFuncChannel);
for i = 1:nFuncChannel
    for j = 1:nPlanes
        tcName = [p.Results.datapath sep tcFile{i,j}];
        tempTC = load(tcName);
        tempTC = tempTC.tempTC;
        % TEMP solution: file input specifies which file in concat TC to load
        if ~isempty(p.Results.file)
            tempTC = tempTC(:,1:nFrames_oneplane(p.Results.file+1,j));
        end
        
        if size(TC{i},1)>0
            TC{i} = [TC{i} [tempTC';nan(1,size(tempTC,1))]];
        else
            TC{i} = [TC{i} tempTC']; %TC: time by neuron
        end
        neuronEachPlane{i}(j) = size(tempTC,1);
    end
end

end