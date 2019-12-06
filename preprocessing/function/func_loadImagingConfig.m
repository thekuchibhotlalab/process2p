function imgConfig = func_loadImagingConfig(filename,varargin)

    p = inputParser;
    p.addParameter('root', pwd)
    p.addParameter('split', true)
    p.parse(varargin{:});

    %filename = ['config\imaging\' filename];
    %tempTxt = importdata([p.Results.Root '\' filename]);

    fidi = fopen( filename );
    tempTxt = textscan(fidi, '%s %s');
    %tempTxt = tempTxt{1};
    fclose( fidi );
    
    imgConfig = struct();
    for i= 1:length(tempTxt{1})
        %tempSplit = strsplit(tempTxt{i});
        rowName = tempTxt{1}{i};
        rowValue = tempTxt{2}{i};
        multipleEntry = {'allChannel, functionalChannel, roiType'};
        if p.Results.split && contains(multipleEntry,rowName)
            tempSplit = strsplit(rowValue,',');
            imgConfig.(rowName) = tempSplit;
        else 
            imgConfig.(rowName) = rowValue;
        end
    end
    
end