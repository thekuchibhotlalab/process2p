function TC = extractTCfromBin(mouse,disk)

global info;

if disk==1
    path1 = 'X:\LabData1\celine';
elseif disk==4
    path1 = 'U:\LabData4\celine';
elseif disk==5
    path1 = 'T:\LabData5';
elseif disk==0
    path1 = 'H:\celine';
else
    error('Sure about the disk? Check the code.');
end

if disk==0
    path = [path1 '\' mouse '\h5\'];
else
    path = [path1 '\' mouse '\h52\'];
end
cd(path);
files = dir('*.h5');
nFiles = length(files);
names = cell(nFiles,1);
for i=1:nFiles, names{i} = files(i).name(1:end-3); end
if disk~=0
    cd ..
else
    path = ['W:\LabData4\celine\' mouse];cd(path);
end
nFrames = nan(nFiles,1);nFrames_add = nan(nFiles,1);
for i=1:nFiles
    sbxread(names{i},1,1);
    nFrames(i) = info.max_idx;    
    if i>1, nFrames_add(i) = nFrames(i)+nFrames_add(i-1); else, nFrames_add(i) = nFrames(i); end
end
nPlanes = info.otparam(3);
nFrames_oneplane = nan(nFiles,nPlanes);
% nFrames_oneplane(logical(mod(nFrames_add,2)),:) = [round(nFrames_add(logical(mod(nFrames_add,2)))/nPlanes) round(nFrames_add(logical(mod(nFrames_add,2)))/nPlanes)-1];
% nFrames_oneplane(~mod(nFrames_add,2),:) = [nFrames_add(~mod(nFrames_add,2))/nPlanes nFrames_add(~mod(nFrames_add,2))/nPlanes];
% nFrames_oneplane = [[0 0];nFrames_oneplane];
nFrames_oneplane(logical(mod(nFrames,2)),:) = [round(nFrames(logical(mod(nFrames,2)))/nPlanes) round(nFrames(logical(mod(nFrames,2)))/nPlanes)-1];
nFrames_oneplane(~mod(nFrames,2),:) = [nFrames(~mod(nFrames,2))/nPlanes nFrames(~mod(nFrames,2))/nPlanes];
nFrames_oneplane = [[0 0];nFrames_oneplane];

if disk==0
    path = [path1 '\' mouse '\suite2p\'];
else
    path = [path1 '\' mouse '\h52\suite2p\'];
end
cd(path);
TC = cell(nPlanes,1);
%-------------------------- neuropil code
neuroPil = cell(nPlanes,1);
%-------------------------- neuropil code

for i=1:nPlanes
    cd([path 'plane' num2str(i-1)]);    
    data = load('Fall.mat');
    ly = data.ops.Ly;
    lx = data.ops.Lx;
%     cd(['H:\plane' num2str(i-1)]); % for speed purpose, access from local disk
    fileID = fopen('data.bin','r'); % open binary file    
    roiName = [ path 'plane' num2str(i-1) '\' mouse '_roi' int2str(i-1) '.zip'];
    rois = ReadImageJROI(roiName); %read imagej rois
    nCells = size(rois,2);
    nFramesPlane = size(data.F,2);
    tempTC = nan(nCells, nFramesPlane);
    
    %-------------------------- neuropil code
    tempNeuroPil = nan(nCells, nFramesPlane);
    weight = getNeuropilWeight([lx,ly], rois);
    %-------------------------- neuropil code
    kall = 0;    
    for j=1:nFiles
        k=0;
        nimg = nFrames_oneplane(j+1,i);
        blksize = 7000;%2000; % nb of frames loaded at a time (depend on RAM)
        to_read = min(blksize,nimg-k);  
        while to_read>0
            A = fread(fileID,ly*lx*to_read,'*int16');
            A = reshape(A,lx,ly,[]); 
            avgA = mean(A,3); % figure;subplot(2,2,1);imagesc(avgA);colormap gray;
            for c=1:nCells
                x = rois{1,c}.mnCoordinates(:,1); %freehand rois have the outlines in x-y coordinates
                y = rois{1,c}.mnCoordinates(:,2);
%                 subplot(2,2,2);imagesc(avgA);colormap gray;hold on;patch(y,x,'g','EdgeColor','none');
                bw = roipoly(avgA,y,x); % y=col; x=rows
                [xt,yt]= find(bw>0);
%                 tempTC(c,kall+1:kall+to_read) = squeeze(mean(mean(A(xt,yt,:),2)));
                tempROIAct = nan(to_read,length(xt));
                
                %-------------------------- neuropil code
                neuroPil = double(A) .* repmat(weight{c},[1 1 to_read]);
                tempNeuroPil(c,kall+1:kall+to_read) = squeeze(sum(sum(neuroPil,1),2));
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
    TC{i} = tempTC;
    %-------------------------- neuropil code
    neuroPil{i} = tempNeuroPil;
    %-------------------------- neuropil code
    
    %save([mouse '_TC_plane' num2str(i-1) '.mat'],'tempTC');
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


        