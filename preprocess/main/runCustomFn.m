function runCustomFn(varargin)
% runCustomFn(varargin)
%
% Input Arguments -- 
% varargin: The set of imaging and sesion parameters from GUI
%
% GUI Input--
% Select one or multiple sessions to load the TC of those sessions.
% Specify your custom functions using fn1, fn2, and fn3
%
% Usage --
% 1. TC is loaded as a cell, each entry is the TC of each session
% 2. Custom functions are ran directly after loading the TC. 
% 3. Custom functions can directly call any variables in the main script
%    (e.g. TC, roisCoord, p.Results.nPlanes (from input parsers))
% 4. Custom functions are ran in the order of fn1, fn2, fn3

p = func_createInputParser();
p.parse(varargin{:});
%---------CHECK NUMBER OF FRAMES IN SBX FILE-----------
global nPlanes

nPlanes = str2double(p.Results.nPlanes);
filename = p.Results.filename;

%---------GET RELEVANT PARAMETERS-----------
[nFuncChannel, functionalChannel, roiType] = func_getFuncChanRoiType(varargin{:});

%---------GET CALCIUM TRACES-----------
[TC, neuronEachPlane, roisCoord] = func_loadTCRoi(varargin{:});

%---------CHECK IF NUMBER OF CHANNELS IS CORRECT-----------
if size(TC,1) ~= nFuncChannel
    disp('ERROR - Number of functional channels not correct')
end
%---------EVALUATE FUNCTIONS -----------
if ~isempty(p.Results.fn1); disp('Function 1 detected and running!'); feval(p.Results.fn1);end
if ~isempty(p.Results.fn2); disp('Function 2 detected and running!'); feval(p.Results.fn2); end
if ~isempty(p.Results.fn3); disp('Function 3 detected and running!'); feval(p.Results.fn3); end

end