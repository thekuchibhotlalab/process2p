function getMeanImg(varargin)

p = func_createInputParser();
p.parse(varargin{:});
sep = '\';

%---------GET RELEVANT PARAMETERS-----------
global info
[nFuncChannel, functionalChannel, roiType] = func_getFuncChanRoiType(varargin{:});
filenames = strsplit(p.Results.filename);
nFiles = length(filenames);
nPlanes = str2double(p.Results.nPlanes);
mouse = p.Results.mouse;
sbxpath = p.Results.sbxpath;
suite2ppath = p.Results.suite2ppath;
h5path = p.Results.h5path;
datapath = p.Results.datapath;
roiFile = p.Results.roiFile;
nFrames_oneplane = p.Results.nFrames_oneplane;
%---------CHECK NUMBER OF CHANNELS-----------
data = load([p.Results.sbxpath sep filenames{1} '.mat']); 
infosbx = data.info;
if isempty(infosbx.otparam)
    check_nPlanes = 1;
else
    check_nPlanes = infosbx.otparam(3);
end
if check_nPlanes ~= nPlanes
    disp('ERROR - nPlanes in the parameter not consistent with sbx file.')
    pause;
end
%---------GET RELEVANT PARAMETERS-----------
[nFuncChannel, functionalChannel, roiType] = func_getFuncChanRoiType(varargin{:});

%---------SAVE MEAN IMAGE OF ALL SESSIONS-----------
imgPath = [savepath sep 'meanImg' sep]; % 2 channel does not work, suite2p only saves meanImg for 2nd channel
for i = 1:nPlanes
    for j = 1:nFuncChannel
        load([suite2ppath sep 'plane' int2str(j-1) sep 'Fall.mat'],'ops');
        mkdir(imgPath)
        maxProj = zeros(size(ops.meanImgE));
        yStart = ops.xrange(1);
        xStart = ops.yrange(1);
        maxProj(xStart+1:xStart+size(ops.max_proj,1), yStart+1:yStart+size(ops.max_proj,2)) = double(ops.max_proj);
        
        if nFuncChannel>1
            imwrite(double(ops.meanImgE),[savepath sep 'meanImgE_plane' int2str(i-1) '_' functionalChannel{chan} '.tiff'])
            imwrite(uint16(maxProj),[savepath sep 'maxProj_plane' int2str(i-1) '_' functionalChannel{chan} '.tiff'])
        else
            imwrite(double(ops.meanImgE),[savepath sep 'meanImgE_plane' int2str(i-1) '.tiff'])
            imwrite(uint16(maxProj),[savepath sep 'maxProj_plane' int2str(i-1) '.tiff'])
        end
    end
end

imshowFlag = true;
avg_sessions = cell(nFiles,2);
for i=1:nPlanes
    for j = 1:nFuncChannel
        data = load([suite2ppath sep 'plane' num2str(i-1) sep 'Fall.mat']);
        ly = data.ops.Ly;
        lx = data.ops.Lx;
        fileID = fopen([suite2ppath sep 'plane' num2str(i-1) sep 'data.bin'],'r'); % open binary file
        for j=1:nFiles
            k=0; a=1;
            nimg = nFrames_oneplane(j+1,i);
            blksize = 7000;%2000; % nb of frames loaded at a time (depend on RAM)
            to_read = min(blksize,nimg-k);  
            avgA = []; avgA = nan(lx,ly);
            while to_read>0
                A = fread(fileID,ly*lx*to_read,'*int16');
                A = reshape(A,lx,ly,[]);
                avgA(:,:,a) = mean(A,3);
                a=a+1;
                k = k+to_read;
                to_read = min(blksize,nimg-k);
            end
            % to check the bloc averages:   figure;for l=1:9,subplot(3,3,l);hold on;imagesc(avgA(:,:,l));colormap gray;title(num2str(l));end
            avg_sessions{j,i} = mean(avgA,3)';
        end   
        % Save the mean img of this plane   
        sessionMeanImg = avg_sessions(:,i);
        if nFuncChannel>1
            save([savepath sep mouse '_MeanImgPerSessions_Plane' num2str(i) '_' functionalChannel{chan} '.mat'],'sessionMeanImg');
        else
            save([savepath sep mouse '_MeanImgPerSessions_Plane' num2str(i) '.mat'],'sessionMeanImg');
        end
        fclose all;  
    end  
end

end