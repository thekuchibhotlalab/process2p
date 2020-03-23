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
p.addParameter('dataType', 'suite2p')
p.addParameter('roiType', 'cell')
p.addParameter('roiMethod', 'suite2p')
p.addParameter('tcMethod', 'suite2p')
p.addParameter('roiFile', [])
p.addParameter('tcFile', [])
p.addParameter('filenameTCFlag', [])
p.addParameter('filenameTC', [])
p.addParameter('sep', '\')

end