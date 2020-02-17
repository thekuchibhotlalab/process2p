function [TC,neuronEachPlane] = func_loadTCmanual(varargin)
%func_loadTCmanual - Description
%
% Syntax: [TC,neuronEachPlane] = func_loadTCmanual(varargin)
%
% Long description
    
p = func_createInputParser();
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
nFrames_oneplane_select = nFrames_oneplane(logical(p.Results.filenameTCFlag),:);
nFrames_oneplane_all = nFrames_oneplane;

nFrames_oneplane_all = cumsum(nFrames_oneplane_all);
nFrames_oneplane_all = [zeros(1,nPlanes);nFrames_oneplane_all];

tcFileSplit = strsplit(p.Results.tcFile);
tcFile = reshape(tcFileSplit,nFuncChannel,nPlanes);
TC = cell(nFuncChannel,nPlanes);
for i = 1:nFuncChannel
    TC_thisPlane = cell(1,nPlanes);
    for j = 1:nPlanes
        tcName = [p.Results.datapath sep tcFile{i,j}];
        tempTC = load(tcName);
        try 
            tempTC = tempTC.TC;
        catch
            tempTC = tempTC.tempTC;
            disp('TC naming follows old convention')
        end
        
        % if the target file is part of the TC, only extract that part
        if ~all(p.Results.filenameTCFlag)
            frameIndex_thisPlane = [];
            fileIndex = find(p.Results.filenameTCFlag==1);
            for k = fileIndex
                frameIndex_thisPlane = [frameIndex_thisPlane ...
                    (nFrames_oneplane_all(k,j)+1):nFrames_oneplane_all(k+1,j)]; 
            end
            tempTC = tempTC(:,frameIndex_thisPlane); 
        end
        TC_thisPlane{j} = tempTC;
        neuronEachPlane{i}(j) = size(tempTC,1);
    end
    TC_thisPlane = func_attachNanFrames(TC_thisPlane, nFrames_oneplane_select,varargin{:});
    TC(i,:) = TC_thisPlane;
end

end