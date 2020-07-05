function AlignAllMovies(mouse,disk)

if disk==1
    path1 = 'X:\LabData1\celine';
elseif disk==4
    path1 = 'U:\LabData4\celine';
else
    error('Sure about the disk? Check the code.');
end

path = [path1 '\' mouse];
cd(path);

% Defaults values
nFrameInBloc = 100;
transformationType = '[Rigid Body]';
nPlanes = 2;
channel = 1;
x_crop_start = 110; %x-value to start strips
x_crop_end = 0;
y_crop_start = 100; %y-value to start strips
y_crop_end = 0;
% unalignMov = 'U:\LabData4\celine\cd017\preprocessData\movie\cd017_000_000_green_unaligned_plane1.avi';
% fileName = 'cd017_000_000';
% path = 'U:\LabData4\celine\cd017\';
% cd(path);

% Load mean img movie of all sessions in Fiji
pathMeanImgs = [path '\preprocessData\ROI\RefImg\'];
cd(pathMeanImgs);
javaaddpath ('E:\KishoreLab\Shared\Matlab\preprocessing\MIJI\mij.jar')
javaaddpath ('E:\KishoreLab\Shared\Matlab\preprocessing\MIJI\ij-1.52i.jar')
javaaddpath ('E:\KishoreLab\Shared\Matlab\preprocessing\MIJI\StackReg_.jar')
javaaddpath ('E:\KishoreLab\Shared\Matlab\preprocessing\MIJI\MultiStackReg1.45_.jar')
MIJ.start('E:\KishoreLab\Shared\Matlab\Fiji.app')
meanImgPlane1 = [pathMeanImgs mouse '_AllMeanImg_Align_Plane1.tif'];
MIJ.run('Open...',['path=' meanImgPlane1]); % open in Fiji
% Choose which image to be used as a template
tempStr = sprintf('Choose the template frame:');
prompt = tempStr;
UItitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,UItitle,dims);
templateIdx(1) = str2double(answer{1});
MIJ.run("Make Substack...",['slices=' num2str(templateIdx(1))]);  % extract template frame
img = MIJ.getCurrentImage; 
z_img_plane1 = uint16(repmat(img,1,1,nFrameInBloc)); % do a z-stack for subsequent alignment

% Same for plane 2
meanImgPlane2 = [pathMeanImgs mouse '_AllMeanImg_Align_Plane2.tif'];
MIJ.run('Open...',['path=' meanImgPlane2]); % open in Fiji
answer = inputdlg(prompt,UItitle,dims);
templateIdx(2) = str2double(answer{1}); 
MIJ.run("Make Substack...",['slices=' num2str(templateIdx(2))]);  % extract template frame
img = MIJ.getCurrentImage; % do a z-stack for subsequent alignment
z_img_plane2 = uint16(repmat(img,1,1,nFrameInBloc));

% Load the z-stacks in Fiji
MIJ.run("Close All");
MIJ.createImage('Template_plane1',z_img_plane1,1); 
MIJ.createImage('Template_plane2',z_img_plane2,1); 

% Align all movies (all sessions)
cd(path);
sbxfiles = dir('*.sbx');
nFiles = length(cellfun(@numel,{sbxfiles.name}));

for i=1:nFiles
    % Load movie
    d = sbxreadmmap(sbxfiles(i).name(1:13));
    img = d.Data.img;
    img = squeeze(img);
    img = intmax('uint16')-img;
    img = permute(img,[2 1 3]);    
    max_idx = size(img,3);
    
    planeMaxIndex = zeros(1,nPlanes);
    imgStack = cell(1,nPlanes);
    for j=1:nPlanes
        planeMaxIndex(j) = floor(max_idx/nPlanes)+ double( j <= mod(max_idx,nPlanes));
        planeIndex = j:nPlanes:max_idx;
        imgStack{:,j} = img(:,:,planeIndex); % split
        imgStack{:,j} = imgStack{:,j}(x_crop_start:(end-x_crop_end),y_crop_start:(end-y_crop_end),:); % crop
        
        nBlocs = floor(max_idx/nPlanes/nFrameInBloc); % because of space limitation in Fiji
        for k=1:nBlocs
            if k==nBlocs
                bloc = imgStack{:,j}(:,:,(k-1)*nFrameInBloc+1:max_idx/nPlanes);
                nFrames = size(bloc,3);   
            else
                bloc = imgStack{:,j}(:,:,(k-1)*nFrameInBloc+1:k*nFrameInBloc);
            end
                     
%             comb = [];comb(:,:,1) = template;comb(:,:,2:nFrames) = bloc; % combine template and movie
            MIJ.createImage('ToAlign2',bloc,1); % load this movie
            fileNameAlign = [path '\preprocessData\alignment\' sbxfiles(i).name(1:13) '_plane' num2str(j) '_test_stackreg_alignment.txt'];    
%             MIJ.run('MultiStackReg',['stack_1=ToAlign action_1=Align file_1=[] transformation=' transformationType ' save=' fileNameAlign]);
%             MIJ.run('MultiStackReg',['stack_1=Template_plane' num2str(j) ' action_1=Align file_1=[] transformation=[Rigid Body] save=' fileNameAlign]);
            MIJ.run('MultiStackReg',['stack_1=Template_plane' num2str(j) ' action_1=[Use as Reference] file_1=[] stack_2=ToAlign action_2=[Align to First Stack] file_2=[] transformation=[Rigid Body] save=' fileNameAlign]);
        end
    end
end





    
