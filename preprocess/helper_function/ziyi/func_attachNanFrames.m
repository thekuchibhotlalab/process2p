function TC = func_attachNanFrames(TC)
nPlanes = size(TC,2);
nFiles = size(TC,1);
for i = 1:nFiles
    nFrames_plane1 = size(TC{i,1},2);
    for j = 1:nPlanes
        nFrames_thisPlane = size(TC{i,j},2);
        if nFrames_thisPlane < nFrames_plane1
            TC{i,j} = cat(2,TC{i,j},nan(size(TC{i,j},1),nFrames_plane1-nFrames_thisPlane));      
        end
    end
end
 

end