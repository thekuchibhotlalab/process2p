function TCnew = func_attachNanFrames(TC, nFrames_oneplane,varargin)
nPlanes = size(nFrames_oneplane,2);
nFiles = size(nFrames_oneplane,1);
TCnew = cell(1,nPlanes);
nFrames_oneplane_cumsum = cumsum(nFrames_oneplane,1);
nFrames_oneplane_cumsum = [zeros(1,nPlanes);nFrames_oneplane_cumsum];
for i = 1:nPlanes
    tempTC = nan(size(TC{i},1),max(sum(nFrames_oneplane,1)));
    startFrame = 0;
    for j = 1:nFiles
        tempTC(:,startFrame+1:startFrame+nFrames_oneplane(j,i)) = ...
            TC{i}(:,nFrames_oneplane_cumsum(j,i)+1:...
            nFrames_oneplane_cumsum(j+1,i));
        
        startFrame = startFrame + max(nFrames_oneplane(j,:));
    end
    TCnew{i} = tempTC;
end

end