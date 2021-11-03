function [roisMask, roisCoord, neuronEachPlane] = func_loadRoi(varargin)

p = func_createInputParser();

p.parse(varargin{:});
sep = '\';

% FILENAME OF THE FINAL REVISION OF ROIS
if p.Results.redrawFile
    redrawFile = 'roi\roi_redrawn_plane0_final.mat roi\roi_redrawn_plane1_final.mat';
end
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
        if ~p.Results.redrawFile
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
        else
            roiFileSplit = strsplit(redrawFile);
            roiFile = reshape(roiFileSplit,nFuncChannel,nPlanes);
            roisCoord = cell(nFuncChannel,nPlanes);
            roisMask = cell(nFuncChannel,nPlanes);
            for i = 1:nFuncChannel
                for j = 1:nPlanes
                    roiName = [p.Results.datapath sep roiFile{i,j}];
                    load(roiName,'roi_redrawn'); % here tempRoi has size of cell*day*3
                    roisCoord{i,j} = cell(size(roi_redrawn,2),size(roi_redrawn,1)); % size of day * cell
                    roisMask{i,j} = cell(size(roi_redrawn,2),size(roi_redrawn,1)); % size of day * cell
                    for k = 1:size(roi_redrawn,1) % number of neurons in this plane
                        % flip coord 1 and 2 here since roi mask flips x and y
                        % xlen = 697; ylen = 403;
                        tempidx = squeeze(roi_redrawn(k,:,2));% tempimdx has size of 1*day
                        tempCoord = squeeze(roi_redrawn(k,:,3));
                        roisCoord{i,j}(:,k) = tempCoord;
                        for l = 1:size(roi_redrawn,2)
                            tempMask = zeros(xlen, ylen);
                            for m = 1:size(tempidx{l},1)
                                tempMask(tempidx{l}(m,1),tempidx{l}(m,2)) = 1;
                            end
                            roisMask{i,j}{l,k} = logical(tempMask);  
                        end
                    end
                end
                neuronEachPlane{i}(j) = size(roi_redrawn,1);
            end
        end
end

end