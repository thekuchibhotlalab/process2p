clear;
dataPath = 'J:\automatic_tracking\';
mouse = 'cd016';
nPlane = 1;
matPath = [dataPath filesep mouse '_plane' int2str(nPlane) '.mat'];
xlen = 697; ylen = 403;

roiData = load(matPath); temp = fields(roiData);
roiData = roiData.(temp{1});


for i = 1:length(roiData.roi)
    stat = regionprops(roiData.roi{i});
    croi(i,:) = round(stat.Centroid);
    
    
    
end 
roiData.croi = croi;


hold on;

extraMargin = 20;
if (ylimm(1)-extraMargin) < 1  || (xlimm(1)-extraMargin) < 1 ||...
        (ylimm(2)+extraMargin)> size(avg_day{k},1) || (xlimm(2)+extraMargin) > size(avg_day{k},2)
    extraMargin = min([ylimm(1)-1, xlimm(1)-1 size(avg_day{k},1)-ylimm(2)   size(avg_day{k},2)-xlimm(2)]);
    marginFlag = 1;
end
% code for original image enhancement
%localImg = uint16(avg_day{k}(ylimm(1)-extraMargin:ylimm(2)+extraMargin, xlimm(1)-extraMargin:xlimm(2)+extraMargin));
%localImgList{k} = imadjust(adapthisteq(localImg, 'NBins', 256)); 

%localImgList{k} = avg_day{k}(ylimm(1)-extraMargin:ylimm(2)+extraMargin, xlimm(1)-extraMargin:xlimm(2)+extraMargin);
%imagesc(imadjust(adapthisteq(localImg, 'NBins', 256),[0 1],[0 1],0.2)); 
localImgList{k} = avg_day{k}(ylimm(1)-extraMargin:ylimm(2)+extraMargin, xlimm(1)-extraMargin:xlimm(2)+extraMargin);
imagesc(localImgList{k}); 
xlim([extraMargin+1 xlimm(2)-xlimm(1)+extraMargin+1]);ylim([extraMargin+1 ylimm(2)-ylimm(1)+extraMargin+1]);
axis off;
title(['D' num2str(k)]);
colormap gray;
plot(extraMargin+1+croi(2)-xlimm(1),extraMargin+1+croi(1)-ylimm(1),'r+'); % plot(1-697, 1-403)
pat_day{k} = patch(yroi-round(xlimm(1)-extraMargin)+1, xroi-round(ylimm(1)-extraMargin)+1,'g', 'FaceColor','None');
        