mouse = 'cd017';
% suite2ppath = 'H:\celine\cd017\suite2p\';
suite2ppath = 'D:\labData\GNG_excit\data\cd017\';
nPlanes = 2;
% 
% signals = cell(nPlanes,1);
for i=1:nPlanes
   cd([suite2ppath 'plane' num2str(i-1)]);
   tc = load([mouse '_TC_plane' num2str(i-1) '.mat']);
   signals{i} = tc.tempTC;
   for j=1:nFiles
       submat = signals{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i));
       signals{i}(:,nFrames_oneplane(j,i)+1:nFrames_oneplane(j+1,i),2) = ...
           (submat-median(submat,2))./median(submat,2);%./repmat(max(submat,[],2)-min(submat,[],2),1,size(submat,2));
   end
end

%%

for p=1:nPlanes
    signals_thisplane = signals{p};
    nCells = size(signals_thisplane,1); 
%     for i=86:nCells
    for i=[97]'
        trace = signals_thisplane(i,:,2);
        figure;plot(trace);title(['cell ' num2str(i)]);
        load([suite2ppath '\plane' int2str(p-1) '\ishere_plane'  int2str(p-1) '.mat'])
        here = ishere{p}(i,:);
      
        nFrames = size(trace,2);
        restrict = 180000:280000; 
        figure;plot(trace(restrict));title(['cell ' num2str(i)]);
        [c, s, options] = deconvolveCa(trace(restrict),'foopsi'); % time x one cell
% outputs:
%   c: T x 1 vector, denoised trace
%   s: T x 1 vector, deconvolved signal
%   b: fluorescence baseline
%   kernel: struct variable containing the parameters for the selected
%       convolution model
%   lambda: Optimal Lagrange multiplier for noise constraint under L1 penalty
%     """olves the noise constrained sparse nonnegat

        figure;
        
        subplot(3,1,1); hold on;
        plot(trace(restrict));
        title('Raw trace (median centered)');
        xlabel('frames'); ylabel('ddf');
        
        subplot(3,1,2); hold on;
        plot(c);
        title('Denoised trace');
        xlabel('frames'); ylabel('ddf');
        
        subplot(3,1,3); hold on;
        plot(s);
        title('Deconvolved trace');
        xlabel('frames'); ylabel('spike rate');
        
    end
end






