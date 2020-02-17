function [TC, neuroPil] = extractTCfromBinZZNew(mouse, varargin)%,disk)

global info;

root = pwd;

% default parameters if none given
p = inputParser;
p.addParameter('infoRoot', root)
p.addParameter('Separator', sep)
p.parse(varargin{:});

% get the info file for the mouse
allMouseInfo = importData([root 'mouseInfo.txt']);
if contains(fieldnames(allMouseInfo),mouse)
    mouseInfo = allMouseInfo.(mouse);
else
    disp(['ERROR:' mouse ' not in the mouseInfo.txt file'])
    pause;
end


%if strcmp(mouse,'cd017')
%    suite2ppath = 'H:\celine\cd017\suite2p\';
%    h5path = 'W:\LabData4\celine\cd017\h5\'; % h5path = 'H:\celine\cd017\h5\';
%    sbxpath = 'W:\LabData4\celine\cd017\';
%end

%if strcmp(mouse,'cd036')
%    suite2ppath = 'H:\celine\cd036\suite2p\';
%    h5path = 'H:\celine\cd036\'; 
%    sbxpath = 'T:\LabData5\cd036\';
%end

%if strcmp(mouse,'cd040')
%    suite2ppath = 'H:\celine\cd040\suite2p\';
%    h5path = 'H:\celine\cd040\'; 
%    sbxpath = 'T:\LabData5\cd040\';
%end

cd(h5path);
files = dir('*.h5');
nFiles = length(files);
names = cell(nFiles,1);
for i=1:nFiles, names{i} = files(i).name(1:end-3); end

cd(sbxpath);
% sbxread(names{2},1,1);
% nPlanes = info.otparam(3);
% nFrames_oneplane = nan(nFiles,nPlanes);

nFrames = nan(nFiles,1);%nFrames_add = nan(nFiles,1);
for i=1:nFiles
    sbxread(names{i},1,1);
    nFrames(i) = info.max_idx;    
%     if i>1
%         nFrames_add(i) = nFrames(i)+nFrames_add(i-1); 
%     else
%         nFrames_add(i) = nFrames(i); 
%     end
%     
%     if mod(nFrames(i),2)
%         if i>1
%             nFrames_oneplane(i,:) = [round(nFrames(i)/nPlanes)+nFrames_oneplane(i-1,1) round(nFrames(i)/nPlanes)-1+nFrames_oneplane(i-1,2)];
%         else
%             nFrames_oneplane(i,:) = [round(nFrames(i)/nPlanes) round(nFrames(i)/nPlanes)-1];
%         end
%     else
%         nFrames_oneplane(i,:) = [nFrames(i)/nPlanes+nFrames_oneplane(i-1,1) nFrames(i)/nPlanes+nFrames_oneplane(i-1,2)];
%     end        
end
% nFrames_oneplane = [[0 0];nFrames_oneplane];

% nFrames_oneplane = nan(nFiles,nPlanes);
% nFrames_oneplane(logical(mod(nFrames_add,2)),:) = [round(nFrames_add(logical(mod(nFrames_add,2)))/nPlanes) round(nFrames_add(logical(mod(nFrames_add,2)))/nPlanes)-1];
% nFrames_oneplane(~mod(nFrames_add,2),:) = [nFrames_add(~mod(nFrames_add,2))/nPlanes nFrames_add(~mod(nFrames_add,2))/nPlanes];
% nFrames_oneplane = [[0 0];nFrames_oneplane];

nPlanes = info.otparam(3);
% Here, the nb of frames / plane is NOT cumulative
nFrames_oneplane = nan(nFiles,nPlanes);
nFrames_oneplane(logical(mod(nFrames,2)),:) = [round(nFrames(logical(mod(nFrames,2)))/nPlanes) round(nFrames(logical(mod(nFrames,2)))/nPlanes)-1];
nFrames_oneplane(~mod(nFrames,2),:) = [nFrames(~mod(nFrames,2))/nPlanes nFrames(~mod(nFrames,2))/nPlanes];
nFrames_oneplane = [zeros(1,nPlanes);nFrames_oneplane];

cd(suite2ppath);
TC = cell(nPlanes,1);
%-------------------------- neuropil code
neuroPil = cell(nPlanes,1);
%-------------------------- neuropil code

for i=1:nPlanes
    cd([suite2ppath 'plane' num2str(i-1)]);    
    data = load('Fall.mat');
    
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
    roiName = [suite2ppath 'plane' num2str(i-1) '\' mouse '_roi' int2str(i-1) '.zip'];
    rois = ReadImageJROI(roiName); %read imagej rois
    nCells = size(rois,2);
    nFramesPlane = size(data.F,2);
    tempTC = nan(nCells, nFramesPlane);
    kall = 0;
    
    %-------------------------- neuropil code
    tempNeuroPil = nan(nCells, nFramesPlane);
    weight= getNeuropilWeight([lx,ly], rois, 'minRadii', 2);
    %-------------------------- neuropil code
    
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
%                 subplot(2,2,2);imagesc(avgA);colormap gray;hold on;patch(y,x,'g','EdgeColor','none');
                bw = roipoly(avgA,y,x); % y=col; x=rows
                [xt,yt]= find(bw>0);
                tempROIAct = nan(to_read,length(xt));
                
                %-------------------------- neuropil code
                %neuroPil = double(A) .* repmat(weight{c},[1 1 to_read]);
                %tempNeuroPil(c,kall+1:kall+to_read) = 
                % note: this code only works for taking the mean
                weightFlag = (weight{i}>0);
                neu = reshapeA(reshape(weightFlag,lx*ly,[]),:);
                tempNeuroPil(c,kall+1:kall+to_read)= sum(neu,1)./sum(weightFlag(:)) * 0.7;
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
    tempTC_removed = tempTC - tempNeuroPil;
    %-------------------------- neuropil code
    TC{i} = tempTC;
    try
        save([mouse '_TC_removed_neuropil_exclude2_plane' num2str(i-1) '.mat'],'tempTC_removed','-v7.3');
    catch
        keyboard
    end
    %-------------------------- neuropil code
    neuroPil{i} = tempNeuroPil;
    %-------------------------- neuropil code
    fclose all;
    
    
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




        