function dff = fn_getDff(rawTC,varargin)
%fn_getDff - Description
%
% Syntax: dff = fn_getDff(rawTC,varargin)
%
% Long description


p = inputParser;
p.KeepUnmatched = true;
p.addParameter('method', 'rollingMedian')
p.addParameter('dffWindow', 1000)
p.addParameter('baselineCorrectionPostDff', false)
p.addParameter('baselineCorrectionWindow', 1000)
p.parse(varargin{:});

switch p.Results.method
case 'movMean'
    baseline = smoothdata(rawTC,1,'movmean',p.Results.dffWindow);
    dff = rawTC./ baseline - 1;

case 'mean'
    %baseline = mean(rawTC,1);
    %dff = rawTC ./ repmat(baseline,[size(rawTC,1) 1]);
    baseline = mean(rawTC,1);
    dff = rawTC ./ repmat(baseline,[size(rawTC,1) 1]) - 1;
end

if p.Results.baselineCorrectionPostDff
    dff = dff - smoothdata(dff,1,'movmean',p.Results.baselineCorrectionWindow);
end

    
end