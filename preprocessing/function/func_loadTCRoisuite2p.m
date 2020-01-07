function [TC, neuronEachPlane, roisCoord] = func_loadTCRoisuite2p(varargin)

p = func_createInputParser();
p.addParameter('file', [])
p.parse(varargin{:});

%---------GET NUMBER OF CHANNELS-----------
nPlanes = str2double(p.Results.nPlanes);
if iscell(p.Results.functionalChannel)
    nFuncChannel = length(p.Results.functionalChannel);
    functionalChannel = p.Results.functionalChannel;
    roiType = p.Results.roiType{1};
else
    nFuncChannel = 1;
    functionalChannel = {p.Results.functionalChannel};
    roiType = p.Results.roiType;
end
nFrames_oneplane = p.Results.nFrames_oneplane;
%[nFrames_oneplane,nFrames,nFrames_add, nPlanes] = func_getSbxFrames('mouse',p.Results.mouse,'root',p.Results.root,'filename',{filename},'sbxpath',p.Results.sbxpath);

cd(p.Results.suite2ppath);
TC = []; neuronEachPlane = nan(nPlanes,1);
roisCoord = cell(1,nPlanes);

for i=1:nPlanes
    cd([p.Results.suite2ppath '\plane' num2str(i-1)]);
    data = load('Fall.mat'); 
    tc = data.F;
    
    switch roiType
    case 'axon'
        iscellFlag = ~data.iscell(:,1);
    case 'cell'
        iscellFlag = data.iscell(:,1);
    end
    
    % if the target file is part of the TC, only extract that part
    if ~all(p.Results.filenameTCFlag)
        frameIndex_thisPlane = [];
        fileIndex = find(p.Results.filenameTCFlag==1);
        for k = fileIndex
            frameIndex_thisPlane = [frameIndex_thisPlane ...
                (nFrames_oneplane(k,j)+1):nFrames_oneplane(k+1,j)]; 
        end
        tempTC = tempTC(logical(iscellFlag),frameIndex_thisPlane);
        
        % TEMP solution: 2nd plane generally have 1 less frame. 
        if size(TC,1)>0
            TC = [TC [tc(logical(iscellFlag),:)';nan(1,sum(iscellFlag))]];
        else
            TC = [TC tc(logical(iscellFlag),:)']; % time by neuron
        end     
    end
    
    neuronEachPlane(i) = sum(iscellFlag);
    roisCoord{i} = data.stat(logical(iscellFlag))';
end
