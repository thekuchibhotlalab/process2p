function Bin2sessionAvg(mouse)%,disk)

% Compute Mean Image per session (from binary file)

global info;

if strcmp(mouse,'cd017')
    suite2ppath = 'H:\celine\cd017\suite2p\';
    h5path = 'U:\LabData4\celine\cd017\h5\'; % h5path = 'H:\celine\cd017\h5\';
    sbxpath = 'U:\LabData4\celine\cd017\';
end

if strcmp(mouse,'cd036')
    suite2ppath = 'H:\celine\cd036\suite2p\';
    h5path = 'H:\celine\cd036\'; 
    sbxpath = 'T:\LabData5\cd036\';
end

if strcmp(mouse,'cd040')
    suite2ppath = 'H:\celine\cd040\suite2p\';
    h5path = 'H:\celine\cd040\'; 
    sbxpath = 'T:\LabData5\cd040\';
end

if strcmp(mouse,'cd037')
    suite2ppath = 'T:\LabData5\cd037\\suite2p\';
    h5path = 'T:\LabData5\cd037\'; 
    sbxpath = 'T:\LabData5\cd037\';
end


cd(h5path);
files = dir('*.h5');
nFiles = length(files);
names = cell(nFiles,1);
for i=1:nFiles, names{i} = files(i).name(1:end-3); end

cd(sbxpath);
% sbxread(names{1},1,1);
% nPlanes = info.otparam(3);
% nFrames_oneplane = nan(nFiles,nPlanes);
nFrames = nan(nFiles,1);
for i=1:nFiles
    sbxread(names{i},1,1);
    nFrames(i) = info.max_idx;
%     if mod(nFrames(i),2)
%         nFrames_oneplane(i,:) = [round(nFrames(i)/nPlanes) round(nFrames(i)/nPlanes)-1];
%     else
%         nFrames_oneplane(i,:) = [nFrames(i)/nPlanes nFrames(i)/nPlanes];
%     end    
end
% nFrames_oneplane = [[0 0];nFrames_oneplane];
nPlanes = info.otparam(3);
% Here, the nb of frames / plane is NOT cumulative
nFrames_oneplane = nan(nFiles,nPlanes);
nFrames_oneplane(logical(mod(nFrames,2)),:) = [round(nFrames(logical(mod(nFrames,2)))/nPlanes) round(nFrames(logical(mod(nFrames,2)))/nPlanes)-1];
nFrames_oneplane(~mod(nFrames,2),:) = [nFrames(~mod(nFrames,2))/nPlanes nFrames(~mod(nFrames,2))/nPlanes];
nFrames_oneplane = [[0 0];nFrames_oneplane];


imshow = true;
avg_sessions = cell(nFiles,2);
for i=1:nPlanes
    cd([suite2ppath 'plane' num2str(i-1)]);
    data = load('Fall.mat');
    ly = data.ops.Ly;
    lx = data.ops.Lx;
    fileID = fopen('data.bin','r'); % open binary file
    for j=1:nFiles
        k=0; a=1;
        nimg = nFrames_oneplane(j+1,i);
        blksize = 7000;%2000; % nb of frames loaded at a time (depend on RAM)
        to_read = min(blksize,nimg-k);  
        avgA = []; avgA = nan(lx,ly);
        while to_read>0
            A = fread(fileID,ly*lx*to_read,'*int16');
            A = reshape(A,lx,ly,[]);
            avgA(:,:,a) = mean(A,3);
            a=a+1;
            k = k+to_read;
            to_read = min(blksize,nimg-k);
        end
        % to check the bloc averages:   figure;for l=1:9,subplot(3,3,l);hold on;imagesc(avgA(:,:,l));colormap gray;title(num2str(l));end
        avg_sessions{j,i} = mean(avgA,3)';
    end   
    % Save the mean img of this plane   
    tosave = avg_sessions(:,i);
    save([mouse '_MeanImgPerSessions_Plane' num2str(i) '.mat'],'tosave');
    fclose all;  
    % Plot mean image each sessions if asked
    if imshow
        try
            figure; 
            for j=1:49
                subplot(7,7,j);
                imagesc(avg_sessions{j,i});colormap gray; 
                title(num2str(j));
            end
        catch
            continue;
        end
    end
end

% %% To check after
% 
% for i=1:nPlanes
%     cd([suite2ppath 'plane' num2str(i-1)]);
%     avgs = load(['MeanImgPerSessions_Plane' num2str(i) '.mat']);
%     avgs = avgs.tosave;
%     nAvgs = size(avgs,1);
%     figure;
%     for j=1:49
%         subplot(7,7,j);
%         imagesc(avgs{j});colormap gray; 
%         title(num2str(j));
%     end
% end


