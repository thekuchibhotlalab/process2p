function f = fn_figureSmartDim(varargin)
%figureSmatrDim - Description
%
%Syntax: f = figureSmatrDim(varargin)
%Input Arguments (varagin):
%   hSize, wSize, hStart, wStart
%   widthHeightRatio, screenRatio
%   Units, visible
%Example: 
%   fn_figureSmartDim('hSize',0.4, 'widthHeightRatio', 3);
%Long description:
%   Create figures of a particular width-height ratio (e.g. for multiple subplots)
    
screenRatio = get(0, 'ScreenSize'); screenRatio = squeeze(screenRatio);
screenRatio = screenRatio(3) / screenRatio(4);

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('widthHeightRatio', 1);
p.addParameter('screenRatio', screenRatio);
p.addParameter('wSize', []);
p.addParameter('hSize', []);
p.addParameter('wStart', 0.0);
p.addParameter('hStart', 0.0);
p.addParameter('Units', 'Normalized');
p.addParameter('visible', 'on');
p.parse(varargin{:});

widthHeightRatio =  screenRatio / p.Results.widthHeightRatio;

if isempty(p.Results.wSize) && isempty(p.Results.hSize)
    hSize = 0.5; wSize = hSize / widthHeightRatio;
elseif isempty(p.Results.hSize)
    wSize = p.Results.wSize; hSize = wSize * widthHeightRatio;
elseif isempty(p.Results.wSize)
    hSize = p.Results.hSize; wSize = hSize / widthHeightRatio;
end

if strcmp(p.Results.Units,'Normalized')
    hStart = 1 - p.Results.hStart - hSize;
else
    hStart = p.Results.hStart; disp('Dimension for non-normalized not implemented');
end

screenSize = [p.Results.wStart, hStart, wSize, hSize];
f = figure('visible', p.Results.visible);
set(f, 'Units',p.Results.Units, 'OuterPosition', screenSize);

end