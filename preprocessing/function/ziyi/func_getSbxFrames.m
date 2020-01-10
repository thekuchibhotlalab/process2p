function [nFrames_oneplane,nFrames,nFrames_add, nPlanes] = func_getSbxFrames(varargin)

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('mouse',[])
p.addParameter('filename','all')
p.addParameter('root',pwd)
p.addParameter('sbxpath',[])
p.parse(varargin{:});

currPath = pwd;
global info
if strcmp(p.Results.filename,'all')
    filename = func_loadMouseConfig(p.Results.mouse,'ReturnData','ImagingFile','root',p.Results.root);
else
    filename = p.Results.filename;
    if ~iscell(filename)
        filename = {filename};

    end
end
%sbxpath = p.Results.sbxpath;
sbxpath = func_loadMasterConfig('mouse',p.Results.mouse, 'ReturnData', 'sbxpath','root',p.Results.root);
cd(sbxpath);
nFiles = length(filename);
nFrames = nan(nFiles,1);nFrames_add = nan(nFiles,1);
sbxread(filename{1},1,1);
nPlanes = info.otparam(3);


for i=1:nFiles
    sbxread(filename{i},1,1);
    nFrames(i) = info.max_idx;    
    if i>1
        nFrames_add(i) = nFrames(i)+nFrames_add(i-1); 
    else
        nFrames_add(i) = nFrames(i); 
    end
    
    if mod(nFrames(i),2)
        if i>1
            nFrames_oneplane(i,:) = [round(nFrames(i)/nPlanes)+nFrames_oneplane(i-1,1) round(nFrames(i)/nPlanes)-1+nFrames_oneplane(i-1,2)];
        else
            nFrames_oneplane(i,:) = [round(nFrames(i)/nPlanes) round(nFrames(i)/nPlanes)-1];
        end
    else
        nFrames_oneplane(i,:) = [nFrames(i)/nPlanes+nFrames_oneplane(i-1,1) nFrames(i)/nPlanes+nFrames_oneplane(i-1,2)];
    end        
end
nFrames_oneplane = [zeros(1,nPlanes);nFrames_oneplane];

cd(currPath)
end