function func_extractTC(fileName,varargin)


p = inputParser;
p.addParameter('roiType', 'imageJ')
p.addParameter('dataType', 'suite2p')
p.addParameter('neuroPil','mean')
p.addParameter('neuroPil','mean')
p.addParameter('separator', sep) 
p.parse(varargin{:});

switch p.Results.dataType
case 'suite2p'
    data = load(['Fall.mat']);
    
    % check if nb of frames per plane computed is true in suite2p files
    nFrameThisPlane = size(data.Fneu,2);
    if nFrameThisPlane==sum(nFrames_oneplane(:,i))
        disp(['Correct nb of frame for plane ' num2str(i) '. Good to go!']);
    else
        error(['Bad nb of frame computed for plane ' num2str(i) ', check it out!']);
    end
    
    ly = data.ops.Ly;
    lx = data.ops.Lx;

    fileID = fopen('data.bin','r'); % open binary file   

case 'sbx'
    disp('Work in progress!')
end

switch p.Results.roiType

case 'imageJ'
    roiName = [mouse '_roi' int2str(i-1) '.zip'];
    rois = ReadImageJROI(roiName); %read imagej rois

case 'suite2p'
    disp('Work in progress!')
end


end