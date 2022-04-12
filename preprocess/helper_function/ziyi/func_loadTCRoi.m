function [TC, neuronEachPlane, roisCoord] = func_loadTCRoi(varargin)

p = func_createInputParser();
p.parse(varargin{:});
sep = '\';

nPlanes = str2double(p.Results.nPlanes);

[nFuncChannel, functionalChannel, roiType] = func_getFuncChanRoiType(varargin{:});

disp('----Loading Timecourse and ROI----')
tic;
switch p.Results.roiMethod
    case 'suite2p' % suite2p 
        [suite2pTC, neuronEachPlane, tempRoi] = func_loadTCRoisuite2p(varargin{:});
        roisCoord = cell(1,nPlanes); % only 1 functional channel for suite2p
        for j = 1:nPlanes % number of planes
            roisCoord{1,j} = cell(1,length(tempRoi{j}));
            for k = 1:length(tempRoi{j}) % neuron in each plane 
                bound = boundary(double(tempRoi{j}{k}.xpix)', double(tempRoi{j}{k}.ypix)',1); % restricted bound
                tempCoord = [tempRoi{j}{k}.xpix(bound)' tempRoi{j}{k}.ypix(bound)'];
                roisCoord{1,j}{k} = tempCoord;
            end
        end
    case 'manual'
        roiFileSplit = strsplit(p.Results.roiFile);
        roiFile = reshape(roiFileSplit,nFuncChannel,nPlanes);
        roisCoord = cell(nFuncChannel,nPlanes);
        for i = 1:nFuncChannel
            for j = 1:nPlanes
                roiName = [p.Results.datapath sep roiFile{i,j}];
                tempRoi = ReadImageJROI(roiName);
                roisCoord{i,j} = cell(1,length(tempRoi));
                for k = 1:length(tempRoi) % number of neurons in this plane
                    % flip coord 1 and 2 here since roi mask flips x and y
                    roisCoord{i,j}{k} = tempRoi{k}.mnCoordinates;
                end
            end
            
        end
end
switch p.Results.tcMethod
    case 'suite2p'
        TC = suite2pTC;
        neuronEachPlane = {neuronEachPlane};
    case 'manual'
        [TC, neuronEachPlane] = func_loadTCmanual(varargin{:}); 
        neuronEachPlane = {neuronEachPlane};
end
totalTime = toc;
disp(['Loading Complete. Time = ' num2str(totalTime,'%.2f')]);

end