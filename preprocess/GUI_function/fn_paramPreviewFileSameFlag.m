function [filenameSelFlag,filenameFrames] = fn_paramPreviewFileSameFlag(configTable,filenames,imagingConfig,targEntry)

filenameIdx = cellfun(@(x)(find(strcmp(configTable.ImagingFile,x))),filenames);

targEntrySelected = configTable.(targEntry)(filenameIdx);
%if length(tcFileSelected)>1 && ~isequal(tcFileSelected{:}); disp('ERROR: Selected files does not have same tcFile!'); end
%targEntrySelected_removeNone = targEntrySelected(~strcmp(targEntrySelected,'None'));
[nUnique,~,~]=unique(targEntrySelected);
if length(nUnique) == 1
    disp('Selected files have the SAME TC.')
    filenameFlagCommon = strcmp(configTable.(targEntry),nUnique{1});
    % nFrames_oneplane_load is used to load TC file corresponding to selected files 
    nFrames_oneplane_Str = configTable.nFrames_oneplane(filenameFlagCommon);
    % filename_loadFlag is used to identify where selected files 
    filenameSelFlag = zeros(size(filenameFlagCommon)); filenameSelFlag(filenameIdx) = 1; 
    %filenameSel = configTable.ImagingFile(logical(filenameSelFlag));
    filenameSelFlag = logical(filenameSelFlag(filenameFlagCommon));

    if ~iscell(nFrames_oneplane_Str); nFrames_oneplane_Str = {nFrames_oneplane_Str}; end
    for i = 1:length(nFrames_oneplane_Str)
        nFrames_oneplane_StrSplit = strsplit(nFrames_oneplane_Str{i});
        for j = 1:str2double(imagingConfig.nPlanes)
            filenameFrames(i,j) = str2double(nFrames_oneplane_StrSplit{j}); %#ok<*AGROW>
        end
    end

else
    disp('Selected files have the DIFFERENT TC. Loading different TC for each individual session.');
    msgbox('UNDER CONSTRUCTION!');
end



end