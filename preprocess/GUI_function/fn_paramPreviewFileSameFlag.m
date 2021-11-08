function [filenameSelIdx,filenameFrames,fileMultiLoadFlag] = fn_paramPreviewFileSameFlag(configTable,filenames,imagingConfig,targEntry)

filenameIdx = cellfun(@(x)(find(strcmp(configTable.ImagingFile,x))),filenames);

targEntryThisFile = configTable.(targEntry)(filenameIdx);
%if length(tcFileSelected)>1 && ~isequal(tcFileSelected{:}); disp('ERROR: Selected files does not have same tcFile!'); end
%targEntrySelected_removeNone = targEntrySelected(~strcmp(targEntrySelected,'None'));
[nUnique,~,~]=unique(targEntryThisFile);
if length(nUnique) == 1
    disp('Selected files have the SAME TC.')
    filenameFlagCommon = strcmp(configTable.(targEntry),nUnique{1});
    % nFrames_oneplane_load is used to load TC file corresponding to selected files 
    filenamesCommon = configTable.ImagingFile(filenameFlagCommon);
    nFrames_oneplaneCommon = configTable.nFrames_oneplane(filenameFlagCommon);
    
    filenameSelIdx = cellfun(@(x)(find(strcmp(filenamesCommon,x))),filenames);
    filenameFrames = fn_readnFramesStr(nFrames_oneplaneCommon,imagingConfig.nPlanes);
    fileMultiLoadFlag = false;
else
    disp('Selected files have the DIFFERENT TC. Loading different TC for each individual session.');
    for i = 1:length(filenameIdx)
        targEntryThisFile = configTable.(targEntry)(filenameIdx(i));
        filenameFlagCommon = strcmp(configTable.(targEntry),targEntryThisFile);
        % nFrames_oneplane_load is used to load TC file corresponding to selected files 
        filenamesCommon = configTable.ImagingFile(filenameFlagCommon);
        filenameSelIdx(i) = find(strcmp(filenamesCommon,filenames{i}));
        nFrames_oneplaneCommon = configTable.nFrames_oneplane(filenameFlagCommon);
        filenameFrames{i} = fn_readnFramesStr(nFrames_oneplaneCommon,imagingConfig.nPlanes);  %#ok<*AGROW>
    end
    fileMultiLoadFlag = true;
end

end