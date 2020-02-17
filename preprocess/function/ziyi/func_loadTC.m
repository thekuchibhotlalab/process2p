function TC = loadTC(mouse)

p = inputParser;
p.addParameter('session','all')
p.parse(varargin{:});

if iscell(p.Results.Session)


elseif strcmp(p.Results.Session,'all')
    cd(h5path);
    files = dir('*.h5');
    nFiles = length(files);
    names = cell(nFiles,1);

    for i=1:nFiles, names{i} = files(i).name(1:end-3); end

    cd(sbxpath);
    nFrames = nan(nFiles,1);nFrames_add = nan(nFiles,1);
    sbxread(names{1},1,1);
    nPlanes = info.otparam(3);
    for i=1:nFiles
        sbxread(names{i},1,1);
        nFrames(i) = info.max_idx;    
        if i>1
            nFrames_add(i) = nFrames(i)+nFrames_add(i-1); 
        else
            nFrames_add(i) = nFrames(i); 
        end
        
        if mod(nFrames(i),2)
            if i>1
                nFrames_oneplane(i,:) = [round(nFrames(i)/nPlanes)+nFrames_oneplane(i-1,1) round(nFrames(i)/nPlanes)-1+nFrames_oneplane(i-1,2)];
            else
                nFrames_oneplane(i,:) = [round(nFrames(i)/nPlanes) round(nFrames(i)/nPlanes)-1];
            end
        else
            nFrames_oneplane(i,:) = [nFrames(i)/nPlanes+nFrames_oneplane(i-1,1) nFrames(i)/nPlanes+nFrames_oneplane(i-1,2)];
        end        
    end
    nFrames_oneplane = [zeros(1,nPlanes);nFrames_oneplane];

    
else 

    
end


end