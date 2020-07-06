function getMeanImg(varargin)

p = func_createInputParser();
p.parse(varargin{:});
sep = p.Results.sep;

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
savepath = p.Results.savepath;
nFrames_oneplane = p.Results.nFrames_oneplane;
alignMethod = p.Results.alignMethod;
if ~isempty(alignMethod); alignMethod = strsplit(alignMethod); end 
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
mkdir(imgPath)
for i = 1:nPlanes
    for chan = 1:nFuncChannel
        load([suite2ppath sep 'plane' int2str(i-1) sep 'Fall.mat'],'ops');
        
        maxProj = zeros(size(ops.meanImgE));
        yStart = ops.xrange(1);
        xStart = ops.yrange(1);
        maxProj(xStart+1:xStart+size(ops.max_proj,1), yStart+1:yStart+size(ops.max_proj,2)) = double(ops.max_proj);
        
        if nFuncChannel>1
            imwrite(double(ops.meanImgE),[imgPath sep 'meanImgE_plane' int2str(i-1) '_' functionalChannel{chan} '.tiff'])
            imwrite(uint16(maxProj),[imgPath sep 'maxProj_plane' int2str(i-1) '_' functionalChannel{chan} '.tiff'])
        else
            imwrite(double(ops.meanImgE),[imgPath sep 'meanImgE_plane' int2str(i-1) '.tiff'])
            imwrite(uint16(maxProj),[imgPath sep 'maxProj_plane' int2str(i-1) '.tiff'])
        end
    end
end

avg_sessions = cell(nFiles,nFuncChannel,nPlanes);
for i=1:nPlanes
    for chan = 1:nFuncChannel
        data = load([suite2ppath sep 'plane' num2str(i-1) sep 'Fall.mat']);
        ly = data.ops.Ly;
        lx = data.ops.Lx;
        if chan == 1
            fileID = fopen([suite2ppath sep 'plane' num2str(i-1) sep 'data.bin'],'r'); % open binary file
        else
            fileID = fopen([suite2ppath sep 'plane' num2str(i-1) sep 'data_chan' int2str(chan) '.bin'],'r'); % open binary file
        end
        for j=1:nFiles
            if  ~isempty(alignMethod) && ~strcmp(alignMethod{j},'suite2p')
                loadFilename = [alignMethod{j} sep filenames{j} '_plane' int2str(i-1) '_stackreg.h5'];
                tempA = h5read(loadFilename,'/data');
                avg_sessions{j,i,chan} = mean(tempA,3);
            else
                k=0; a=1;
                nimg = nFrames_oneplane(j,i);
                blksize = 5000;%2000; % nb of frames loaded at a time (depend on RAM)
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
                % to check the bloc averages: 
                avg_sessions{j,i,chan} = mean(avgA,3)';
            end
        end   
        % Save the mean img of this plane   
        sessionMeanImg = avg_sessions(:,i);
        if nFuncChannel>1
            save([imgPath sep mouse '_MeanImgPerSessions_Plane' num2str(i-1) '_' functionalChannel{chan} '.mat'],'sessionMeanImg');
        else
            save([imgPath sep mouse '_MeanImgPerSessions_Plane' num2str(i-1) '.mat'],'sessionMeanImg');
        end
        fclose all;  
    end  
end

end