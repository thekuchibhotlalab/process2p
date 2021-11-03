function paramValue = fn_paramPreviewAddTableParam(configTable,filenameIdx,paramValue,imagingConfig)

varNames = configTable.Properties.VariableNames;

for i = 1:length(varNames)
    tempParam = configTable.(varNames{i})(filenameIdx);
    
    if strcmp(varNames{i},'nFrames_oneplane')
        saveParam = [];
        if ~iscell(tempParam); tempParam = {tempParam}; end
        for k = 1:length(tempParam)
            nFrames_oneplaneSplit = strsplit(tempParam{k});
            for j = 1:str2double(imagingConfig.nPlanes)
                saveParam(k,j) = str2double(nFrames_oneplaneSplit{j});
            end
        end  
    else
        saveParam = tempParam;
    end
    
    if iscell(saveParam) && length(saveParam) > 1
        
        if isequal(saveParam{:}); saveParam = saveParam{1};
        else;saveParam = strjoin(saveParam);
            if ~strcmp(varNames{i},'BehavFile') && ~strcmp(varNames{i},'ImagingFile') && ~strcmp(varNames{i},'BehavType') 
                disp(['WARNING: Not all files have same ' varNames{i} '!'])
            end
        end
    elseif isnumeric(saveParam)
        saveParam = {saveParam};
    end
    paramValue = [paramValue;varNames{i};saveParam];
end





end