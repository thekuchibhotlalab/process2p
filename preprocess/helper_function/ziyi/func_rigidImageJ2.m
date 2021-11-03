function func_rigidImageJ2(mouse, sessionName,varargin)

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('imageJpath', 'E:\KishoreLab\Shared\Matlab\Fiji.app')
p.addParameter('MIJIpath', 'E:\KishoreLab\Shared\Matlab\preprocessing\MIJI')
p.addParameter('method', 'Rigid Body')
p.addParameter('alignRep', 1)
p.parse(varargin{:});
sep = '\';
MIJIpath = p.Results.MIJIpath;
imageJpath = p.Results.imageJpath;
alignMethod = p.Results.method;

filename = [mouse '_config.csv'];
configTable = readtable(filename);
h5path = configTable.h5path{1};
datapath = configTable.datapath{1};

% Defaults values
nFrameInBloc = 100;
transformationType = '[Rigid Body]';
nPlanes = 2;
channel = 1;

% Load reference img movie of all sessions in Fiji
refImgPath = [datapath sep 'meanImg'];
load([refImgPath sep mouse '_ops_plane0.mat'],'ops');
img_plane{1} = ops.refImg;
load([refImgPath sep mouse '_ops_plane1.mat'],'ops');
img_plane{2} = ops.refImg;

javaaddpath ([ MIJIpath sep 'mij.jar'])
javaaddpath ([ MIJIpath sep 'ij-1.52i.jar'])
javaaddpath ([ MIJIpath sep 'StackReg_.jar'])
javaaddpath ([ MIJIpath sep 'MultiStackRegistration_-1.46.2.jar'])
MIJ.start(imageJpath)

% Load the z-stacks in Fiji

% Align all movies (all sessions)

% Load movie
img = h5read([h5path sep sessionName '.h5'],'/data');
img = permute(img,[2 1 3]);    
max_idx = size(img,3);

planeMaxIndex = zeros(1,nPlanes);
imgStack = cell(1,nPlanes);
savepath = [h5path sep sessionName];
mkdir(savepath);

for j=1:nPlanes
    tic;
    h5name = [h5path sep sessionName '_plane' int2str(j-1) '_stackReg.h5'];
    planeMaxIndex(j) = floor(max_idx/nPlanes)+ double( j <= mod(max_idx,nPlanes));
    planeIndex = j:nPlanes:max_idx;
    imgStack{j} = img(:,:,planeIndex); % split
    img_save = uint16(zeros(size(imgStack{j})));
    MIJ.run("Close All");
    
    
    nBlocs = ceil(max_idx/nPlanes/nFrameInBloc); % because of space limitation in Fiji
    for k=1:nBlocs
        start = (k-1)*nFrameInBloc+1;
        if k==nBlocs            
            bloc = imgStack{j}(:,:,start:end);
            nFrames = size(bloc,3);
            MIJ.selectWindow(['Template_plane' int2str(j)]); MIJ.run('Close')
            sessionRefImg =  uint16(repmat( mean(imgStack{j}(:,:,1:nFrameInBloc),3),1,1,nFrames));
            MIJ.createImage(['Template_plane' int2str(j)],sessionRefImg,1); 
        else    
            nFrames = nFrameInBloc;
            bloc = imgStack{j}(:,:,start:start+nFrames-1);
            if k==1
                sessionRefImg =  uint16(repmat( imgStack{j}(:,:,end),1,1,nFrameInBloc));
                MIJ.createImage(['Template_plane' int2str(j)],sessionRefImg,1); 
            elseif k==2
                MIJ.selectWindow(['Template_plane' int2str(j)]); MIJ.run('Close')
                sessionRefImg =  uint16(repmat( mean(imgStack{j}(:,:,1:nFrameInBloc),3),1,1,nFrameInBloc));
                MIJ.createImage(['Template_plane' int2str(j)],sessionRefImg,1); 
            end
        end
       
        MIJ.createImage('ToAlign',bloc,1); % load this movie
        fileNameAlign = [savepath sep sessionName '_plane' int2str(j-1) '_bloc' num2str(k,'%03d') '.txt'];  
        MIJ.run('MultiStackReg',['stack_1=Template_plane' int2str(j) ' action_1=[Use as Reference] file_1=[]' ...
            ' stack_2=ToAlign action_2=[Align to First Stack] file_2=' fileNameAlign ' transformation=[Translation] save']);
        
        MIJ.selectWindow('ToAlign');
        img_aligned = uint16(MIJ.getCurrentImage);
        imgStack{j}(:,:,start:start+nFrames-1) = img_aligned;
        
        MIJ.selectWindow('ToAlign');MIJ.run('Close');

    end
    MIJ.selectWindow(['Template_plane' int2str(j)]); MIJ.run('Close')
    
    z_img_plane =  uint16(repmat(img_plane{j},1,1,nFrameInBloc));
    MIJ.createImage(['Template_plane' int2str(j)],z_img_plane,1); 
    meanImg = uint16(mean(imgStack{j},3));
    MIJ.createImage(['Mean_plane' int2str(j)],repmat(meanImg,1,1,nFrameInBloc),1);
    
    if p.Results.alignRep > 1
        fileNameAlign = [];
        for i = 1:p.Results.alignRep
            fileNameAlign{i} = [savepath sep sessionName '_plane' int2str(j-1) '_rotation' int2str(i) '.txt'];  
            MIJ.run('MultiStackReg',['stack_1=Template_plane' int2str(j) ...
                    ' action_1=[Use as Reference] file_1=[] stack_2=Mean_plane' int2str(j) ' action_2=[Align to First Stack] file_2='...
                    fileNameAlign{i} ' transformation=[' alignMethod '] save']);
        end
    else
        fileNameAlign = [savepath sep sessionName '_plane' int2str(j-1) '_rotation.txt'];  
        MIJ.run('MultiStackReg',['stack_1=Template_plane' int2str(j) ...
                ' action_1=[Use as Reference] file_1=[] stack_2=Mean_plane' int2str(j) ' action_2=[Align to First Stack] file_2='...
                fileNameAlign ' transformation=[' alignMethod '] save']);
        
    end
        
    MIJ.selectWindow(['Template_plane' int2str(j)]);MIJ.run('Close');
    MIJ.selectWindow(['Mean_plane' int2str(j)]);MIJ.run('Close');
    for k=1:nBlocs
        start = (k-1)*nFrameInBloc+1;
        if k==nBlocs            
            bloc = imgStack{j}(:,:,start:end);
            nFrames = size(bloc,3);
            z_img_plane =  uint16(repmat(img_plane{j},1,1,nFrames));
            MIJ.createImage(['Template_plane' int2str(j)],z_img_plane,1); 
            MIJ.createImage(['Mean_plane' int2str(j)],repmat(meanImg,1,1,nFrames),1);
            if p.Results.alignRep > 1
                fileNameAlign = [];
                for i = 1:p.Results.alignRep
                    fileNameAlign{i} = [savepath sep sessionName '_plane' int2str(j-1) '_rotation_temp' int2str(i) '.txt'];  
                    MIJ.run('MultiStackReg',['stack_1=Template_plane' int2str(j) ...
                        ' action_1=[Use as Reference] file_1=[] stack_2=Mean_plane' int2str(j) ' action_2=[Align to First Stack] file_2='...
                        fileNameAlign{i} ' transformation=[' alignMethod '] save']);
                end
            else
                fileNameAlign = [savepath sep sessionName '_plane' int2str(j-1) '_rotation_temp.txt'];  
                MIJ.run('MultiStackReg',['stack_1=Template_plane' int2str(j) ...
                    ' action_1=[Use as Reference] file_1=[] stack_2=Mean_plane' int2str(j) ' action_2=[Align to First Stack] file_2='...
                    fileNameAlign ' transformation=[' alignMethod '] save']);
            end
            MIJ.selectWindow(['Template_plane' int2str(j)]);MIJ.run('Close');
            MIJ.selectWindow(['Mean_plane' int2str(j)]);MIJ.run('Close');

        else
            nFrames = nFrameInBloc;
            bloc = imgStack{j}(:,:,start:start+nFrames-1);
        end
        MIJ.createImage('ToAlign',bloc,1);
        
        if p.Results.alignRep > 1
            for i = 1:p.Results.alignRep
                MIJ.run('MultiStackReg',['stack_1=ToAlign'...
                    ' action_1=[Load Transformation File] file_1=' fileNameAlign{i} ' stack_2=None action_2=Ignore file_2=[] transformation=[' alignMethod ']']);
            end
        else
            MIJ.run('MultiStackReg',['stack_1=ToAlign'...
                ' action_1=[Load Transformation File] file_1=' fileNameAlign ' stack_2=None action_2=Ignore file_2=[] transformation=[' alignMethod ']']);   
        end
        MIJ.selectWindow('ToAlign');
        img_aligned = uint16(MIJ.getCurrentImage);
        img_save(:,:,start:start+nFrames-1) = img_aligned;
        MIJ.run('Close');
    end
    h5create(h5name,'/data',size(img_save),'DataType','uint16');
    h5write(h5name,'/data',img_save);
    toc;
end
MIJ.exit;
end
