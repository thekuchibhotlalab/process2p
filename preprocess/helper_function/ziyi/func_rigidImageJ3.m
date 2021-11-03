function func_rigidImageJ3(mouse, sessionName,varargin)

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('imageJpath', 'E:\KishoreLab\Shared\Matlab\Fiji.app')
p.addParameter('MIJIpath', 'E:\KishoreLab\Shared\Matlab\preprocessing\MIJI')
p.parse(varargin{:});
sep = '\';
MIJIpath = p.Results.MIJIpath;
imageJpath = p.Results.imageJpath;

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
javaaddpath ([ MIJIpath sep 'MultiStackReg1.45_.jar'])
MIJ.start(imageJpath)

% Load the z-stacks in Fiji

% Align all movies (all sessions)

% Load movie
img = h5read([h5path sep sessionName '.h5'],'/data');
img = permute(img,[2 1 3]);    
max_idx = size(img,3);

planeMaxIndex = zeros(1,nPlanes);
imgStack = cell(1,nPlanes);

for j=1:nPlanes
    
    tic;
    h5name = [h5path sep sessionName '_plane' int2str(j-1) '_stackReg.h5'];
    
     
    planeMaxIndex(j) = floor(max_idx/nPlanes)+ double( j <= mod(max_idx,nPlanes));
    planeIndex = j:nPlanes:max_idx;
    imgStack{j} = img(:,:,planeIndex); % split
    img_save = uint16(zeros(size(imgStack{j})));
    MIJ.run("Close All");
    z_img_plane =  uint16(repmat(img_plane{j},1,1,nFrameInBloc));
    MIJ.createImage(['Template_plane' int2str(j)],z_img_plane,1); 
    
    nBlocs = ceil(max_idx/nPlanes/nFrameInBloc); % because of space limitation in Fiji
    for k=1:nBlocs
        start = (k-1)*nFrameInBloc+1;
        if k==nBlocs            
            bloc = imgStack{:,j}(:,:,start:end);
            nFrames = size(bloc,3);
            MIJ.selectWindow(['Template_plane' int2str(j)]); MIJ.run('Close')
            z_img_plane =  uint16(repmat(img_plane{j},1,1,nFrames));        
            MIJ.createImage(['Template_plane' int2str(j)],z_img_plane,1); 
        else
            nFrames = nFrameInBloc;
            bloc = imgStack{:,j}(:,:,start:start+nFrames-1);
        end
                    
        MIJ.createImage('ToAlign',bloc,1); % load this movie
        fileNameAlign = [h5path sep sessionName '_plane' int2str(j-1) '_stackreg.txt'];  
        MIJ.run('MultiStackReg',['stack_1=Template_plane' int2str(j) ...
            ' action_1=[Use as Reference] file_1=[] stack_2=ToAlign action_2=[Align to First Stack] file_2=[] transformation=[Rigid Body]']);
        
        MIJ.selectWindow('ToAlign');
        img_aligned = uint16(MIJ.getCurrentImage);
        
        img_save(:,:,start:start+nFrames-1)=img_aligned;
        
        MIJ.selectWindow('ToAlign');
        MIJ.run('Close');

    end
   
    h5create(h5name,'/data',size(img_save),'DataType','uint16');
    h5write(h5name,'/data',img_save);
    toc;
end

end
