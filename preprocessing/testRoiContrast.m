i = 1;
k = 7;
js = 20:30;
sbpMargin = 0.01;
extraMargin = 20;
figure; 
for l = 1:length(js)
    j = js(l);
    yroi = rois{1,j}.mnCoordinates(:,1); %freehand rois have the outlines in x-y coordinates
    xroi = rois{1,j}.mnCoordinates(:,2);
    croi = [mean(minmax(xroi)) mean(minmax(yroi))];
    ylimm = [croi(1)-20 croi(1)+20];
    xlimm = [croi(2)-20 croi(2)+20];
    
    subplot_tight(5,length(js),l,sbpMargin)
    imagesc(avg_day{k,i})
    xlim(xlimm);ylim(ylimm);
    axis off;

    subplot_tight(5,length(js),length(js)+l,sbpMargin)
    hold on; imagesc(data.ops.meanImgE)
    xlim(xlimm);ylim(ylimm);
    colormap gray;
    imgg = imadjust(adapthisteq(uint16(avg_day{k,i})));
    axis off;
    
    subplot_tight(5,length(js),length(js)*2+l,sbpMargin)
    hold on; imagesc(imgg);
    xlim(xlimm);ylim(ylimm);
    axis off;

    subplot_tight(5,length(js),length(js)*3+l,sbpMargin); hold on
    localImg = uint16(avg_day{k,i}(ylimm(1):ylimm(2), xlimm(1):xlimm(2)));
    imagesc(imadjust(adapthisteq(localImg, 'NBins', 256),[0 1],[0 1])); xlim([1 41]); ylim([1 41])
    axis off;
    
    subplot_tight(5,length(js),length(js)*4+l,sbpMargin); hold on
    localImg = uint16(avg_day{k,i}(ylimm(1)-extraMargin:ylimm(2)+extraMargin, xlimm(1)-extraMargin:xlimm(2)+extraMargin));
    imagesc(imadjust(adapthisteq(localImg, 'NBins', 256),[0 1],[0 1],0.2)); xlim([1 41]); ylim([1 41])
    xlim([extraMargin xlimm(2)-xlimm(1)+extraMargin]);ylim([extraMargin ylimm(2)-ylimm(1)+extraMargin]);
    axis off;
end