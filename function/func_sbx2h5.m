function func_sbx2h5(sbxpath,targetPath)
% convert sbx to h5 for all files in the directory
cd(sbxpath);
allFileStruct = dir([sbxpath '\*.sbx']);
allFile = {allFileStruct.name};
for i = 1:length(allFile)
    disp(['processing' allFile{i}]); tic;
    sbxname = strsplit(allFile{i},'.');
    sbx2h5cropdir(sbxname{1},targetPath); 
    toc;
end

end