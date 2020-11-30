function cellnum = func_getRedrawNum(redrawPath,checkPlane)
% e.g. path: E:\excitData\cd042\roi\checkCell\revision
% input plane as 0 or 1
allFileStruct = dir([redrawPath '\*.png']);
allFile = {allFileStruct.name};
count = 0;
for i = 1:length(allFile)
    temp = strsplit(allFile{i},'.');
    tempFile = temp{1};
    temp = strsplit(tempFile,'cell');
    if strcmp(temp{2}(end),'_')
        plane =  str2double(temp{2}(end-1));
    else 
        plane =  str2double(temp{2}(end));
    end
    if plane == checkPlane
        count = count+1;
        cellnum(count) = str2double(temp{3});
    end
end
disp(['plane' int2str(plane)])
cellnum = sort(cellnum);
end

