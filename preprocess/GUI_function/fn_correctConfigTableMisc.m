function [configTable] = fn_correctConfigTableMisc(configTable)
    %FN_CORRECTMISCELL Summary of this function goes here
    %   Detailed explanation goes here

    % Correct the exception of oneplane where matlab recognize as number, not string.
    

    
    if isnumeric(configTable.nFrames_oneplane)
        tempInt = configTable.nFrames_oneplane;
        tempStr = cell(size(tempInt));
        for i = 1:length(tempInt)
            tempStr{i} = int2str(tempInt(i));
        end
        configTable.nFrames_oneplane = tempStr;
    end

end

