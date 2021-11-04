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
nFrames_oneplane = p.Results.nFrames_oneplane_TC;
nFrames_oneplane_cumsum = [zeros(1,nPlanes); cumsum(nFrames_oneplane,1)];
%nFrames_oneplane_select = nFrames_oneplane(logical(p.Results.filenameTCFlag),:);


tcFileSplit = strsplit(p.Results.tcFile);
tcFile = reshape(tcFileSplit,nFuncChannel,nPlanes);
TC = cell(nFuncChannel,1); neuronEachPlane = cell(nFuncChannel,1);
for i = 1:nFuncChannel
    for j = 1:nPlanes
        tcName = [p.Results.datapath sep tcFile{i,j}];
        tempTC = load(tcName);
        try 
            tempTC = tempTC.TC;
        catch
            tempTC = tempTC.tempTC;
            disp('TC naming follows old convention')
        end
        
        
        
        fileIndex = p.Results.filenameTCIdx;
        for k = 1:p.Results.nFiles
            frameIndex_thisPlane = (nFrames_oneplane_cumsum(fileIndex(k),j)+1):nFrames_oneplane_cumsum(fileIndex(k)+1,j); 
            TC_thisPlane{k,j} = tempTC(:,frameIndex_thisPlane);
        end

        neuronEachPlane{i}(j) = size(tempTC,1);
        disp(['Loading part of the TCfile. From ' int2str(min(frameIndex_thisPlane)) ' to '...
            int2str(max(frameIndex_thisPlane)) ' frames. Total ' int2str(length(frameIndex_thisPlane)) ' frames.']);
    end
    TC_thisPlane = func_attachNanFrames(TC_thisPlane);
    TC{i} = TC_thisPlane;
end
if nFuncChannel == 1; TC = TC{1}; neuronEachPlane = neuronEachPlane{1}; end

end