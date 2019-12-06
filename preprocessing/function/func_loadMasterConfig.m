function returnData = func_loadMasterConfig(varargin)
    % default parameters if none given
    p = inputParser;
    p.addParameter('root', pwd)
    p.addParameter('mouse', 'all')
    p.addParameter('returnData', 'all')
    p.parse(varargin{:});
    
    filename = 'mouse\mouseMaster.txt';

    tempTxt = importdata([p.Results.root '\' filename]);
    infoData = struct();
    allMouse = {};
    for i = 1:length(tempTxt)
        tempSplit = strsplit(tempTxt{i});
        if strcmp(tempSplit{1}, 'mouse')
            tempMouseName = tempSplit{2}; 
            infoData.(tempMouseName).mouse = tempSplit{2};
            allMouse = [allMouse tempMouseName];
            % we can do this because mouse is always the first name
        else 
            infoData.(tempMouseName).(tempSplit{1}) = tempSplit{2};
        end
    
    end

    % return only the fields that was asked by the input
    if strcmp(p.Results.mouse,'all')
        mouseList = allMouse;
    elseif ~iscell(p.Results.mouse)
        mouseList = {p.Results.mouse};
    else
        mouseList = p.Results.mouse;
        %infoData = infoData.(p.Results.mouse);
    end

    returnData = struct();

    for i = 1:length(mouseList)
        if strcmp(p.Results.returnData,'all')
            returnList =  fieldnames(infoData.(mouseList{i}));
        elseif ~iscell(p.Results.returnData)
            returnList = {p.Results.returnData};
        else
            returnList = p.Results.returnData;
            %infoData = infoData.(p.Results.mouse);
        end
        for j = 1:length(returnList)
            returnData.(mouseList{i}).(returnList{j}) = infoData.(mouseList{i}).(returnList{j});
        end
    end 
    
    if ~iscell(p.Results.mouse) && ~iscell(p.Results.returnData) && ~strcmp(p.Results.returnData,'all')
        returnData = returnData.(p.Results.mouse).(p.Results.returnData);
    end
end