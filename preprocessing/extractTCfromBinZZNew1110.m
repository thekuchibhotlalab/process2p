function [TC, neuroPil] = extractTCfromBinZZNew1110(mouse, varargin)%,disk)

global info;

infoPath = pwd;

% default parameters if none given
p = inputParser;
p.addParameter('infoPath', infoPath)
p.addParameter('Separator', sep)
%p.addParameter('roiType', 'imageJ')
%p.addParameter('dataType', 'suite2p')
%p.addParameter('neuroPil','mean')
%p.addParameter('neuroPil','mean')
p.parse(varargin{:});

% load the metaDataFile
load([mouse '_ImagingMetaData.mat']);
%nFrames_oneplane = func_loadSbxFrames(mouse);
%h5path = func_loadMasterConfig('Mouse',mouse, 'ReturnData', 'h5path');

cd(suite2ppath);
TC = cell(nPlanes,1);
%-------------------------- neuropil code
neuroPil = cell(nPlanes,1);
%-------------------------- neuropil code

for i=1:nPlanes
    %planePath = [mouseInfo.suite2ppath sep 'plane' num2str(i-1) ];
    cd([suite2ppath 'plane' num2str(i-1)]);    
    %[tempTC, tempNeuroPil] = func_extractTC('roiType',p.Results.roiType, 'dataType', p.Results.dataType);

    data = load(['Fall.mat']);
    
    % check if nb of frames per plane computed is true in suite2p files
    nFrameThisPlane = size(data.Fneu,2);
    if nFrameThisPlane==sum(nFrames_oneplane(:,i))
        disp(['Correct nb of frame for plane ' num2str(i) '. Good to go!']);
    else
        error(['Bad nb of frame computed for plane ' num2str(i) ', check it out!']);
    end
    
    ly = data.ops.Ly;
    lx = data.ops.Lx;

    fileID = fopen('data.bin','r'); % open binary file    
    %roiName = [suite2ppath 'plane' num2str(i-1) '\' mouse '_roi' int2str(i-1) '.zip'];
    
    switch imgConfig.ROIMethod
    case 'imageJ'
        roiName = [mouse '_roi' int2str(i-1) '.zip'];
        rois = ReadImageJROI(roiName); %read imagej rois
    case 'suite2p'
        disp('Work in progress!')    
    end
    
    nCells = size(rois,2);
    nFramesPlane = size(data.F,2);
    tempTC = nan(nCells, nFramesPlane);
    kall = 0;
    
    %-------------------------- neuropil code
    switch imgConfig.neuropil
    case 'mean'
        tempNeuroPil = nan(nCells, nFramesPlane);
        weight= getNeuropilWeight([lx,ly], rois, 'minRadii', 2);
    case 'suite2p'
        disp('Work in progress!')
    end
    %-------------------------- neuropil code
    
    % need to add multiple channels
    for j=1:nFiles
        k=0;
        nimg = nFrames_oneplane(j+1,i);
        blksize = 7000;%2000; % nb of frames loaded at a time (depend on RAM)
        to_read = min(blksize,nimg-k);  
        while to_read>0
            A = fread(fileID,ly*lx*to_read,'*int16');
            A = reshape(A,lx,ly,[]);
            %------- added solely for the sake of neuropil
            reshapeA = reshape(A,lx*ly,[]);
            %--------
            avgA = mean(A,3); % figure;subplot(2,2,1);imagesc(avgA);colormap gray;
            for c=1:nCells
                x = rois{1,c}.mnCoordinates(:,1); %freehand rois have the outlines in x-y coordinates
                y = rois{1,c}.mnCoordinates(:,2);
%                subplot(2,2,2);imagesc(avgA);colormap gray;hold on;patch(y,x,'g','EdgeColor','none');
                bw = roipoly(avgA,y,x); % y=col; x=rows
                [xt,yt]= find(bw>0);
                tempROIAct = nan(to_read,length(xt));
                
                %-------------------------- neuropil code
                %neuroPil = double(A) .* repmat(weight{c},[1 1 to_read]);
                %tempNeuroPil(c,kall+1:kall+to_read) = 
                % note: this code only works for taking the mean
                switch p.Results.neuroPil
                case 'mean'
                    weightFlag = (weight{i}>0);
                    neu = reshapeA(reshape(weightFlag,lx*ly,[]),:);
                    tempNeuroPil(c,kall+1:kall+to_read)= sum(neu,1)./sum(weightFlag(:)) * 0.7;
                case 'suite2p'
                    disp('Work in progress!')
                end
                
                %-------------------------- neuropil code
                
                for fr = 1:length(xt)
                    tempROIAct(1:to_read,fr) = squeeze(A(xt(fr),yt(fr),:)); 
                end
                tempTC(c,kall+1:kall+to_read)= mean(tempROIAct,2)'; 
                
            end
            k = k+to_read;
            kall = kall+to_read; % x frames in the whole file
            to_read = min(blksize,nimg-k);
        end  
        disp(['File ' num2str(j) '/' num2str(nFiles) ' DONE! ' num2str(kall) ' frames so far.']);
    end
    %-------------------------- neuropil code
    if ~strcmp(p.Results.neuroPil,'off')
        neuroPil{i} = tempNeuroPil;
        tempTC = tempTC - tempNeuroPil;
    end
    %-------------------------- neuropil code
    TC{i} = tempTC;
    try
        save([mouse '_TC_plane' num2str(i-1) '.mat'],'tempTC','-v7.3');
    catch
        keyboard
    end
    %-------------------------- neuropil code
    
    %-------------------------- neuropil code
    %fclose all;
    
    
%     figure;hold on;
%     imshow(data.ops.meanImgE);colormap gray;
%     for j=1:nCells
%         y = rois{1,j}.mnCoordinates(:,1); %freehand rois have the outlines in x-y coordinates
%         x = rois{1,j}.mnCoordinates(:,2);
%         scatter(y,x,'r.');
%         tempIndex = (yt-1) .* size(avgA,1) + xt;
%                 tempTC(i,frame) = mean(img(tempIndex));
%                 
%     end
end




        