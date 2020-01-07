function [roisMask, roisCoord, neuronEachPlane] = func_loadRoi(varargin)

p = func_createInputParser();
p.addParameter('xlen', 697)
p.addParameter('ylen', 403)
p.parse(varargin{:});
sep = '\';

nPlanes = str2double(p.Results.nPlanes);

[nFuncChannel, functionalChannel, roiType] = func_getFuncChanRoiType(varargin{:});
xlen = p.Results.xlen;
ylen = p.Results.ylen;
switch p.Results.roiMethod
    case 'suite2p' % suite2p 
        [~, neuronEachPlane, tempRoi] = func_loadTCRoisuite2p(varargin{:});
        roisCoord = cell(1,nPlanes); % only 1 functional channel for suite2p
        roisMask = cell(1,nPlanes);
        for j = 1:nPlanes % number of planes
            roisCoord{1,j} = cell(1,length(tempRoi{j}));
            roisMask{1,j} = cell(1,length(tempRoi{j}));
            for k = 1:length(tempRoi{j}) % neuron in each plane 
                bound = boundary(double(tempRoi{j}{k}.xpix)', double(tempRoi{j}{k}.ypix)',1); % restricted bound
                tempCoord = [tempRoi{j}{k}.xpix(bound) tempRoi{j}{k}.ypix(bound)];
                roisCoord{1,j}{k} = tempCoord;
                roisMask{1,j}{k} = [double(tempRoi{j}{k}.xpix), double(tempRoi{j}{k}.ypix)];
            end
        end
    case 'manual'
        roiFileSplit = strsplit(p.Results.roiFile);
        roiFile = reshape(roiFileSplit,nFuncChannel,nPlanes);
        roisCoord = cell(nFuncChannel,nPlanes);
        roisMask = cell(nFuncChannel,nPlanes);
        for i = 1:nFuncChannel
            for j = 1:nPlanes
                roiName = [p.Results.datapath sep roiFile{i,j}];
                tempRoi = ReadImageJROI(roiName);
                roisCoord{i,j} = cell(1,length(tempRoi));
                for k = 1:length(tempRoi) % number of neurons in this plane
                    % flip coord 1 and 2 here since roi mask flips x and y
                    %xlen = 800; ylen = 600;
                    tempMask = roipoly(zeros(xlen,ylen),tempRoi{k}.mnCoordinates(:,2),...
                        tempRoi{k}.mnCoordinates(:,1));   
                    roisMask{i,j}{k} = logical(tempMask);
                    roisCoord{i,j}{k} = tempRoi{k}.mnCoordinates;
                end
            end
            neuronEachPlane{i}(j) = length(tempRoi);
        end
end

end