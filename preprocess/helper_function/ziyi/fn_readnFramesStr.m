function nFrames_oneplane = fn_readnFramesStr(nFrames_oneplaneStr,nPlanes)

if ~iscell(nFrames_oneplaneStr); nFrames_oneplaneStr = {nFrames_oneplaneStr}; end
for i = 1:length(nFrames_oneplaneStr)
    nFrames_oneplane_StrSplit = strsplit(nFrames_oneplaneStr{i});
    for j = 1:str2double(nPlanes)
        nFrames_oneplane(i,j) = str2double(nFrames_oneplane_StrSplit{j}); %#ok<*AGROW>
    end
end

end