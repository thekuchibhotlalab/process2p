function tcFilename = extractRedraw(varargin)

p = func_createInputParser();
p.parse(varargin{:});

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
day = p.Results.day;

TCpath = [datapath filesep 'TC' filesep];
if exist(TCpath,'dir')==0; mkdir(TCpath); end

%---------CHECK NUMBER OF CHANNELS-----------
data = load([p.Results.sbxpath filesep filenames{1} '.mat']); 
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
%---------LOAD ROI BASED ON ROI METHOD-----------
switch p.Results.roiMethod
case 'manual'
    data = load([suite2ppath filesep 'plane0' filesep 'Fall.mat'] );
    ly = data.ops.Ly;
    lx = data.ops.Lx;
    [roisMask, roisCoord, neuronEachPlane] = func_loadRoi(varargin{:},'xlen',lx,'ylen',ly,'redrawFile',true);
case 'suite2p'
    disp('Work in progress!')
end
%---------EXTRACT BASED ON DATA TYPE-----------

tcFilename = cell(nFuncChannel,nPlanes);
switch p.Results.dataType
case 'suite2p' % suite2p only support one functional channel!!!
    nFrames = nan(nFiles,1);%nFrames_add = nan(nFiles,1);
    for i=1:nFiles
        sbxread([sbxpath filesep filenames{i}],1,1);
        nFrames(i) = info.max_idx; 
    end

    for i=1:nPlanes 
        tic;
        nFramesPlane = sum(nFrames_oneplane(:,i));
            
        %check if nb of frames per plane computed is true in suite2p files
        data = load([suite2ppath filesep 'plane' num2str(i-1) filesep 'Fall.mat'] );
        nFramesPlaneCheck = size(data.F,2);
        if nFramesPlaneCheck==nFramesPlane
            disp(['Correct nb of frame for plane ' num2str(i) '. Good to go!']);
        else
            error(['Bad nb of frame computed for plane ' num2str(i) ', check it out!']);
        end
        for chan = 1:nFuncChannel
            if chan == 1
                fileID = fopen([suite2ppath filesep 'plane' num2str(i-1) filesep 'data.bin'],'r'); % open binary file 
            else 
                fileID = fopen([suite2ppath filesep 'plane' num2str(i-1) filesep 'data_chan' int2str(chan) '.bin'],'r'); % open binary file 
            end
            % get roi and preallocate resource
            rois = roisMask{chan,i};
            nCells = size(rois,2);
            weight= getNeuropilWeight([lx,ly], roisCoord{chan,i}, 'minRadii', 2);
            %preallocate size for TC and neuropil
            TC = nan(nCells, nFramesPlane);
            neuroPil = nan(nCells, nFramesPlane);
            kall = 0;
            for j=1:nFiles
                k=0;
                nimg = nFrames_oneplane(j,i);
                blksize = 7000;%2000; % nb of frames loaded at a time (depend on RAM)
                to_read = min(blksize,nimg-k);
                while to_read>0
                    A = fread(fileID,lx*ly*to_read,'*int16');
                    A = reshape(A,lx*ly,[]);
                    for c=1:nCells
                        % get roi and neuropil mask
                        xyflag = logical(reshape(rois{day(j)+1,c},lx*ly,1));
                        weightFlag = reshape((weight{i}>0),lx*ly,1);
                        % extract timecourse
                        neuroPil(c,kall+1:kall+to_read)= mean(A(weightFlag,:),1);
                        TC(c,kall+1:kall+to_read)= mean(A(xyflag,:),1); 
                    end
                    
                    k = k+to_read;
                    kall = kall+to_read; % x frames in the whole file
                    to_read = min(blksize,nimg-k);
                end  
                disp(['File ' num2str(j) '/' num2str(nFiles) ' DONE! ' num2str(kall) ' frames so far.']);
            end
            % save data. file name are different depending on # of channels
            if nFuncChannel>1
                
                save([TCpath filesep mouse '_TC_plane' num2str(i-1) '_' functionalChannel{chan} '_final.mat'],'TC','-v7.3');
                save([TCpath filesep mouse '_neuroPil_plane' num2str(i-1) '_' functionalChannel{chan} '_final.mat'],'neuroPil','-v7.3');
                tcFilename{chan,i} = ['TC' filesep mouse '_TC_plane' num2str(i-1) '_' functionalChannel{chan} '_final.mat'];

            else
                save([TCpath filesep mouse '_TC_plane' num2str(i-1) '_final.mat'],'TC','-v7.3');
                save([TCpath filesep mouse '_neuroPil_plane' num2str(i-1) '_final.mat'],'neuroPil','-v7.3');
                tcFilename{1,i} = ['TC' filesep mouse '_TC_plane' num2str(i-1) '_final.mat'];
            end
            fclose all;
        end
        disp(['Plane ' num2str(i) ' is done, time elapsed: ' num2str(toc) ' sec']);
    end
case 'sbx' % more than 1 functional channel goes in sbx
    disp('Work in progress!')
end   
        
end
    
    