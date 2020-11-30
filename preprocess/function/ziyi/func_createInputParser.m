function p = func_createInputParser()
%myFun - Description
%
% Syntax: output = myFun(input)
%
% Long description

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('mouse', [])
p.addParameter('root', pwd)
p.addParameter('h5path', [])
p.addParameter('suite2ppath', [])
p.addParameter('sbxpath', [])
p.addParameter('savepath', [])
p.addParameter('datapath', [])
p.addParameter('behavpath', [])
p.addParameter('filename', [])
p.addParameter('nPlanes', 2)
p.addParameter('functionalChannel', 'green')
p.addParameter('nFrames_oneplane', [])
p.addParameter('nFrames_oneplane_TC', [])
p.addParameter('nFrames_oneplane_bin', [])
p.addParameter('dataType', 'suite2p')
p.addParameter('roiType', 'cell')
p.addParameter('roiMethod', 'suite2p')
p.addParameter('tcMethod', 'suite2p')
p.addParameter('roiFile', [])
p.addParameter('tcFile', [])
p.addParameter('spikeFile', [])
p.addParameter('filenameTCFlag', [])
p.addParameter('filenameTC', [])
p.addParameter('filenameBinFlag', [])
p.addParameter('filenameBin', [])
p.addParameter('day', [])
p.addParameter('xlen', 697)
p.addParameter('ylen', 403)
p.addParameter('redrawFile', false)
p.addParameter('neuropil', 'false')
p.addParameter('neuropilMethod', 'mean')
p.addParameter('meanImgFlag', [])
p.addParameter('funcParam', 'None')
p.addParameter('alignMethod', [])
p.addParameter('sep', '\')

end