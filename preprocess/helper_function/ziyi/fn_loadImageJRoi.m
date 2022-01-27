roiName = 'cd042_roi0.zip';
xlen = 697; ylen = 403;
tempRoi = ReadImageJROI(roiName);
roisCoord = cell(1,length(tempRoi));

for k = 1:length(tempRoi) % number of neurons in this plane
    % flip coord 1 and 2 here since roi mask flips x and y
    %xlen = 800; ylen = 600;
    tempMask = roipoly(zeros(xlen,ylen),tempRoi{k}.mnCoordinates(:,2),...
        tempRoi{k}.mnCoordinates(:,1));   
    roisMask{k} = logical(tempMask);
    roisCoord{k} = tempRoi{k}.mnCoordinates;
end