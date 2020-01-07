% directory = ('H:\celine\cd017\suite2p\plane1\');

% directory = ('T:\LabData5\cd037\suite2p\plane1\');
% directory = ('H:\celine\cd036\suite2p\plane0\');
% directory = ('H:\celine\cd016\suite2p\plane1\');
directory = ('H:\celine\cd040\suite2p\plane1\');


load([directory 'Fall.mat'],'ops');
mkdir([directory 'meanImgTiff\'])

imwrite(double(ops.meanImgE),[directory 'meanImgTiff\meanImgE.tiff'])

maxProj = zeros(size(ops.meanImgE));
% xStart = floor((size(ops.meanImg,1) - size(ops.max_proj,1))/2);
% yStart = floor((size(ops.meanImg,2) - size(ops.max_proj,2))/2);
% yStart = min(sum(ops.meanImgE(:,1:round(size(ops.meanImgE,2)/2))==0,2));
% xStart = min(sum(ops.meanImgE(1:round(size(ops.meanImgE,2)/2),:)==0,2));
yStart = ops.xrange(1);
xStart = ops.yrange(1);
maxProj(xStart+1:xStart+size(ops.max_proj,1), yStart+1:yStart+size(ops.max_proj,2)) = double(ops.max_proj);

imwrite(uint16(maxProj),[directory 'meanImgTiff\max_proj_fixed.tiff'])

%% 
% 
% meanImg = ops.meanImgE;
% meanImgShift = circshift(meanImg,[-1,5]);
% imwrite(double(meanImgShift),[directory 'meanImgTiff\meanImgE_fixed.tiff'])


