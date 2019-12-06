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
p.addParameter('filename', [])
p.addParameter('nPlanes', 2)
p.addParameter('functionalChannel', 'green')
p.addParameter('nFrames_oneplane', [])
p.addParameter('roiType', 'cell')
p.addParameter('roiMethod', 'suite2p')
p.addParameter('tcMethod', 'suite2p')
p.addParameter('roiFile', [])
p.addParameter('tcFile', [])

end