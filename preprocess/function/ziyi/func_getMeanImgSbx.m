function meanImg = func_getMeanImgSbx(varargin)

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('sbxpath', pwd)
p.parse(varargin{:});

a = dir([p.Results.sbxpath '/*.sbx']);
cd(p.Results.sbxpath);

filenames= {a.name}; 
avgFrames = 1000;
nPlanes = 2;
meanImg = [];

for i = 1:length(filenames); temp = filenames{i} ;
    tic; disp(['Loading File ' temp]);
    frames = sbxread(temp(1:end-4),1,avgFrames * nPlanes); 
    frames = squeeze(frames(:,:,:,1:nPlanes:end)); 
    
    [m,T] = sbxalignxMat(frames,1:avgFrames);
    
    for j=1:avgFrames
        frames_aligned(:,:,j) = circshift(frames(:,:,j),T(j,:));
    end

    meanImg(:,:,i) = uint16(mean(frames_aligned,3)); 
    t = toc; disp(['Time: ' num2str(t,'%.1f') ' secs']);
end

save('meanImg_example_plane0.mat', 'meanImg');

[m,T] = sbxalignxMat(meanImg,1:length(filenames));
    
for j=1:length(filenames)
    meanImg_aligned(:,:,j) = circshift(meanImg(:,:,j),T(j,:));
end
save('meanImg_example_aligned_plane0.mat', 'meanImg_aligned');