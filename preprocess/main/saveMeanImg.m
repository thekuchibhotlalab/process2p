function saveMeanImg(varargin)

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
nFrames_oneplane = p.Results.nFrames_oneplane_bin;
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

if contains(suite2ppath,'None')
    suite2ppath_split = strsplit(suite2ppath);
    suite2ppath_removeNone = suite2ppath_split(~strcmp(suite2ppath_split,'None'));
    suite2ppath = suite2ppath_removeNone{1};
end

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
            imwrite(fn_convert2uint16(ops.meanImgE),[imgPath sep 'meanImgE_plane' int2str(i-1) '_' functionalChannel{chan} '.tiff'])
            imwrite(uint16(maxProj),[imgPath sep 'maxProj_plane' int2str(i-1) '_' functionalChannel{chan} '.tiff'])
        else
            imwrite(fn_convert2uint16(ops.meanImgE),[imgPath sep 'meanImgE_plane' int2str(i-1) '.tiff'])
            imwrite(uint16(maxProj),[imgPath sep 'maxProj_plane' int2str(i-1) '.tiff'])
        end
    end
end

avg_sessions = cell(nFiles,nPlanes,nFuncChannel);
top_sessions1 = cell(nFiles,nPlanes,nFuncChannel);
top_sessions2 = cell(nFiles,nPlanes,nFuncChannel);
for i=1:nPlanes
    for chan = 1:nFuncChannel
        bin_file_count = 0;
        data = load([suite2ppath sep 'plane' num2str(i-1) sep 'Fall.mat']);
        ly = data.ops.Ly;
        lx = data.ops.Lx;
        if chan == 1
            fileID = fopen([suite2ppath sep 'plane' num2str(i-1) sep 'data.bin'],'r'); % open binary file
        else
            fileID = fopen([suite2ppath sep 'plane' num2str(i-1) sep 'data_chan' int2str(chan) '.bin'],'r'); % open binary file
        end
        for j=1:nFiles
            if j == 1; tic; end
            if  ~isempty(alignMethod) && ~strcmp(alignMethod{j},'suite2p')
                loadFilename = [alignMethod{j} sep filenames{j} '_plane' int2str(i-1) '_stackreg.h5'];
                tempA = h5read(loadFilename,'/data');
                avg_sessions{j,i,chan} = mean(tempA,3);
                top_sessions1{j,i,chan} = prctile(tempA,70,3);
                top_sessions2{j,i,chan} = prctile(tempA,85,3);
            else
                k=0; a=1;
                bin_file_count = bin_file_count + 1;
                nimg = nFrames_oneplane(bin_file_count,i);
                blksize = 5000;%2000; % nb of frames loaded at a time (depend on RAM)
                to_read = min(blksize,nimg-k);  
                avgA = nan(lx,ly);
                topA1 = nan(lx,ly);
                topA2 = nan(lx,ly);
                while to_read>0
                    A = fread(fileID,ly*lx*to_read,'*int16');
                    A = reshape(A,lx,ly,[]);
                    avgA(:,:,a) = mean(A,3);
                    topA1(:,:,a) = prctile(A,70,3);
                    topA2(:,:,a) = prctile(A,85,3);
                    a=a+1;
                    k = k+to_read;
                    to_read = min(blksize,nimg-k);
                end
                % to check the bloc averages: 
                avg_sessions{j,i,chan} = mean(avgA,3)';
                top_sessions1{j,i,chan} = prctile(topA1,75,3)';
                top_sessions2{j,i,chan} = prctile(topA2,75,3)';
            end
            if j ==1; tempTime = toc; disp(['Estim. time = ' num2str(tempTime*nFiles,'%.3f') ' secs.']); end 
        end   
        % Save the mean img of this plane   
        sessionMeanImg = avg_sessions(:,i,chan);
        if nFuncChannel>1
            save([imgPath sep mouse '_MeanImgPerSessions_Plane' num2str(i-1) '_' functionalChannel{chan} '.mat'],'sessionMeanImg');
        else
            save([imgPath sep mouse '_MeanImgPerSessions_Plane' num2str(i-1) '.mat'],'sessionMeanImg');
        end
        
        sessionMeanImg = top_sessions1(:,i,chan);
        if nFuncChannel>1
            save([imgPath sep mouse '_Top1ImgPerSessions_Plane' num2str(i-1) '_' functionalChannel{chan} '.mat'],'sessionMeanImg');
        else
            save([imgPath sep mouse '_Top1ImgPerSessions_Plane' num2str(i-1) '.mat'],'sessionMeanImg');
        end
        sessionMeanImg = top_sessions2(:,i,chan);
        if nFuncChannel>1
            save([imgPath sep mouse '_Top2ImgPerSessions_Plane' num2str(i-1) '_' functionalChannel{chan} '.mat'],'sessionMeanImg');
        else
            save([imgPath sep mouse '_Top2ImgPerSessions_Plane' num2str(i-1) '.mat'],'sessionMeanImg');
        end
        fclose all;  
    end  
end

end