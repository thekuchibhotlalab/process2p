function [nFuncChannel, functionalChannel, roiType] = func_getFuncChanRoiType(varargin)

p = func_createInputParser();
p.parse(varargin{:});
sep = '\';

if iscell(p.Results.functionalChannel)
    nFuncChannel = length(p.Results.functionalChannel);
    functionalChannel = p.Results.functionalChannel;
    roiType = p.Results.roiType{1};
else
    nFuncChannel = 1;
    functionalChannel = {p.Results.functionalChannel};
    roiType = p.Results.roiType;
end
