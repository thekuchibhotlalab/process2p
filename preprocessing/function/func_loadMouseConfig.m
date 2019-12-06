function T = func_loadMouseConfig(mouse,varargin)
    
    p = inputParser;
    p.addParameter('returnData','all');
    p.addParameter('root',pwd);
    p.addParameter('filename',[mouse '_config.csv']);
    p.parse(varargin{:});
    
    sep = '\';
    T = readtable([p.Results.root sep 'mouse' sep p.Results.filename]);
    if ~strcmp(p.Results.returnData, 'all')
        T = T{:,p.Results.returnData};
    end
end