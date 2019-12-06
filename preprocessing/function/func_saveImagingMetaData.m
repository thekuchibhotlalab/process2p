function  func_saveImagingMetaData(mouse,varargin)
%myFun - Description
%
% Syntax:  func_saveImagingMetaData = myFun(input)
%
% Long description
    
p = inputParser;
%p.addParameter('filename','all')
p.parse(varargin{:});

[nFrames_oneplane,nFrames, nFrames_add nPlanes] = func_getSbxFrames('mouse',mouse);
imgConfigPath = func_loadMasterConfig('mouse',mouse, 'returnData', 'imagingConfig');
behConfigPath = func_loadMasterConfig('mouse',mouse, 'returnData', 'behaviorConfig');
h5path = func_loadMasterConfig('mouse',mouse, 'returnData', 'h5path');
sbxpath = func_loadMasterConfig('mouse',mouse, 'returnData', 'sbxpath');
suite2ppath = func_loadMasterConfig('mouse',mouse, 'returnData', 'suite2ppath');
behavpath = func_loadMasterConfig('mouse',mouse, 'returnData', 'behavpath');
filename = func_loadMouseConfig(mouse,'returnData','ImagingFile');
% roiname = func_loadMouseConfig(mouse,'returnData','ROIFile');
behavname = func_loadMouseConfig(mouse,'returnData','BehavFile');
nFiles = length(filename);

imgConfig = func_loadImagingConfig(imgConfigPath);
%behConfig = func_getBehaviorConfig(behConfigPath);
nChannels = length(imgConfig.allChannel);

% Code for the date of the imaging sessions 
%goodlastday='';a=0;
%while ~strcmp(goodlastday,lastday_behavior)
%    goodlastday = sbxfiles(end-a).date(1:end-8);
%    a=a+1;
%end



save([mouse '_ImagingMetaData.mat']);
% add functional channel and stuff

end