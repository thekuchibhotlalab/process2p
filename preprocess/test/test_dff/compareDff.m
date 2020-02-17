spikeITI = 20:10:100;
noise = 0.2;
[CaTrace, allParam] = func_getFakeCaTrace(spikeITI,noise);
plotCompareDff(CaTrace);
%%
spikeITI = 80;
noise = 0.2:0.1:1.0;
[CaTrace, allParam] = func_getFakeCaTrace(spikeITI,noise);
peakDff = plotCompareDff(CaTrace,'plot',true);

%%
spikeITI = 45./[0.1:0.1:0.9 1:0.2:1.8 2:0.25:3 3.5:0.5:5];
noise = 0.05:0.05:1.5;
rep = 10;
peakTrue = zeros(length(spikeITI),length(noise),rep);
peakMed = zeros(length(spikeITI),length(noise),rep);
peakDeconv = zeros(length(spikeITI),length(noise),rep);
for j = 1:rep
    tic;
    for i = 1:length(spikeITI)
        [CaTrace, allParam] = func_getFakeCaTrace(spikeITI(i),noise);
        peakDff = plotCompareDff(CaTrace,'plot',false,'smoothWindow',5);
        peakTrue(i,:,j) = peakDff(:,4);
        peakMed(i,:,j) = peakDff(:,1);
        peakDeconv(i,:,j) = peakDff(:,3);
    end
    toc;
end

%%
cmap = redbluecmap;
newCmap = imresize(cmap, [64, 3]);  % original color map contain just 11 colors, this increase it to 64
newCmap = min(max(newCmap, 0), 1);

peakTrueMean = mean(peakTrue,3);
peakMedMean = mean(peakMed,3);
peakDeconvMean = mean(peakDeconv,3);
%[x,y] = meshgrid(noise,spikeITI);

figure(1);
colorlim1 = max(abs(peakTrueMean(:)-peakDeconvMean(:))./peakTrueMean(:));
imagesc((peakTrueMean-peakDeconvMean)./peakTrueMean)
colormap(newCmap)
colorbar
yticks(1:3:length(spikeITI))
yticklabels(round(3./spikeITI(1:3:length(spikeITI))*15,2))
xticks(1:3:length(noise))
xticklabels(noise(1:3:length(noise)))
ylabel('spikes/s'); xlabel('noise(std of baseline)');

figure(2);
colorlim2 = (max(peakMedMean(:))-min(peakMedMean(:)))/max(peakMedMean(:));
imagesc(abs(peakMedMean-max(peakMedMean(:)))/max(peakMedMean(:)))
colormap(newCmap)
colorbar
yticks(1:3:length(spikeITI))
yticklabels(round(3./spikeITI(1:3:length(spikeITI))*15,2))
xticks(1:3:length(noise))
xticklabels(noise(1:3:length(noise)))
ylabel('spikes/s'); xlabel('noise(std of baseline)');

colorlim = max(colorlim1,colorlim2);
figure(1);caxis([-colorlim colorlim])
figure(2);caxis([-colorlim colorlim])
%%
function peakDff = plotCompareDff(CaTrace,varargin)

p = inputParser;
p.addParameter('mixGuassian',false)
p.addParameter('smoothWindow',0)
p.addParameter('plot',true)
p.parse(varargin{:});

nNeuron = size(CaTrace,1);
nFrames = size(CaTrace,2);
nFramesPerTrial= 200;
nTrial = 200;

mu = zeros(nNeuron,2);
sig = zeros(nNeuron,2);

if p.Results.mixGuassian
tic;
for i = 1:nNeuron
    [tempmu,tempsig,t,iter] = fit_mix_gaussian(CaTrace(i,:)',2);
    mu(i,:) = tempmu;
    sig(i,:) = tempsig;
end
toc;
end

%% try deconvolve
tau=0.7;
dt = 1/15;
decayConstant = 0.933;
trace = CaTrace;
if ~(p.Results.smoothWindow==0)
    trace = smoothdata(trace,2,'gaussian',p.Results.smoothWindow);
end

deconvTrace = diff(trace,1,2)+trace(:,1:end-1)*dt/tau;
deconvTrace = [nan(size(trace,1),1) deconvTrace];
deconvTraceZ = (deconvTrace - repmat(nanmean(deconvTrace,2),...
    [1 size(deconvTrace,2)])) ./ repmat(nanstd(deconvTrace,0,2),...
    [1 size(deconvTrace,2)]);
deconvTraceFilter = deconvTraceZ;
deconvTraceFilter(deconvTraceZ < 1.5) = 0;
deconvTraceFrame = round(log(0.1 ./deconvTraceFilter) / log(decayConstant));
deconvTraceFrame(isinf(deconvTraceFrame)) = nan;
framesFlag = ones(size(deconvTrace));
for i = 1:size(deconvTrace,2)
    for j = 1:size(deconvTrace,1)
        if ~isnan(deconvTraceFrame(j,i))
            framesFlag(j,i:i+deconvTraceFrame(j,i)) = 0;
        end
    end
end
framesFlag = framesFlag(:,1:nFrames);
baseline_median = zeros(nNeuron,1);
for i = 1:nNeuron
    temp = trace(i,:);
    temp(~framesFlag(i,:)) = nan;
    baseline_median(i) = nanmedian(temp);
end
%tempTrace = getFakeCaTrace(15,5,30,0,900,500);
%deconvTrace= diff(tempTrace)+tempTrace(1:end-1)*(1-decayConstant);
%deconvTrace = [0; deconvTrace];

%% plot the dff results
%load('synDff.mat')
%figure;
%scatter(median(CaTrace,2),mu(:,1),40,'.');

% get 9 cells with different spontaneous firing
%cellIndex = find(allParam(1,:) == 9 & allParam(2,:) == 1 ... %abs(allParam(4,:) - 1.5) < 1e-4...
%    & allParam(5,:) == 900 & allParam(6,:) == 500);
baseline = 1200;
tempMedian = median(CaTrace,2);
tempPerc = prctile(CaTrace,25,2);
tc_reshape = reshape(CaTrace,nNeuron,...
    nFramesPerTrial,nTrial);
tc_psth = mean(tc_reshape,3);
tc_psth = tc_psth - repmat(mean(tc_psth(:,1:4),2),[1 200]);

dff = (CaTrace-repmat(tempMedian,1,nFrames)) ./...
    repmat(tempMedian,1,nFrames);
dff_reshape = reshape(dff,nNeuron,nFramesPerTrial,nTrial);
dff_psth = mean(dff_reshape,3);
dff_psth = dff_psth - repmat(mean(dff_psth(:,1:4),2),[1 200]);


dff_perc = (CaTrace-repmat(tempPerc,1,nFrames)) ./...
    repmat(tempPerc,1,nFrames);
dff_reshape_perc = reshape(dff_perc,nNeuron,nFramesPerTrial,nTrial);
dff_psth_perc = mean(dff_reshape_perc,3);
dff_psth_perc = dff_psth_perc - repmat(mean(dff_psth_perc(:,1:4),2),[1 200]);

dff_deconv = (CaTrace-repmat(baseline_median,1,nFrames)) ./...
    repmat(baseline_median,1,nFrames);
dff_reshape_deconv = reshape(dff_deconv,nNeuron,nFramesPerTrial,nTrial);
dff_psth_deconv = mean(dff_reshape_deconv,3);
dff_psth_deconv = dff_psth_deconv - repmat(mean(dff_psth_deconv(:,1:4),2),[1 200]);

dff_true = (CaTrace-baseline) ./ baseline;
dff_reshape_true = reshape(dff_true,nNeuron,nFramesPerTrial,nTrial);
dff_psth_true = mean(dff_reshape_true,3);
dff_psth_true = dff_psth_true - repmat(mean(dff_psth_true(:,1:4),2),[1 200]);

if p.Results.mixGuassian
    peakDff = zeros(nNeuron,5);
else
    peakDff = zeros(nNeuron,4);
end
peakDff(:,1) = max(dff_psth,[],2);
peakDff(:,2) = max(dff_psth_perc,[],2);
peakDff(:,3) = max(dff_psth_deconv,[],2);
peakDff(:,4) = max(dff_psth_true,[],2);


if p.Results.mixGuassian
    mu_selected = mu(:,1);
    dff_mg = (CaTrace-repmat(mu_selected,1,nFrames)) ./...
        repmat(mu_selected,1,nFrames);
    dff_reshape_mg = reshape(dff_mg,nNeuron,nFramesPerTrial,nTrial);
    dff_psth_mg = mean(dff_reshape_mg,3);
    dff_psth_mg = dff_psth_mg - repmat(mean(dff_psth_mg(:,1:4),2),[1 200]);
    peakDff(:,5) = max(dff_psth_mg,[],2);
end

%colormat = repmat(0.1:0.1:0.9,3,1);
if p.Results.plot

    for i = 1:nNeuron
        figure(1);
        subplot(3,ceil(nNeuron/3),i)
        temp = trace(i,:);
        temp(~framesFlag(i,:)) = nan;
        plot(trace(i,:)); hold on; plot(temp)
        figure(2);
        subplot(3,3,i)
        histogram(temp)
        baseline_median(i) = nanmedian(temp);
    end


    figure(3);
    if p.Results.mixGuassian
        for i = 1:nNeuron
            subplot(6,2,1);hold on;plot(CaTrace(i,1:1000)');
            subplot(6,2,2);hold on;plot(tc_psth(i,:)');
            subplot(6,2,3);hold on;plot(dff(i,1:1000)');
            subplot(6,2,4);hold on;plot(dff_psth(i,:)');
            subplot(6,2,5);hold on;plot(dff_mg(i,1:1000)');
            subplot(6,2,6);hold on;plot(dff_psth_mg(i,:)');
            subplot(6,2,7);hold on;plot(dff_perc(i,1:1000)');
            subplot(6,2,8);hold on;plot(dff_psth_perc(i,:)');
            subplot(6,2,9);hold on;plot(dff_deconv(i,1:1000)');
            subplot(6,2,10);hold on;plot(dff_psth_deconv(i,:)');
            subplot(6,2,11);hold on;plot(dff_true(i,1:1000)');
            subplot(6,2,12);hold on;plot(dff_psth_true(i,:)');
        end
    else 
        for i = 1:nNeuron
            subplot(5,2,1);hold on;plot(CaTrace(i,1:1000)');
            subplot(5,2,2);hold on;plot(tc_psth(i,:)');
            subplot(5,2,3);hold on;plot(dff(i,1:1000)');
            subplot(5,2,4);hold on;plot(dff_psth(i,:)');
            subplot(5,2,5);hold on;plot(dff_perc(i,1:1000)');
            subplot(5,2,6);hold on;plot(dff_psth_perc(i,:)');
            subplot(5,2,7);hold on;plot(dff_deconv(i,1:1000)');     
            subplot(5,2,8);hold on;plot(dff_psth_deconv(i,:)');        
            subplot(5,2,9);hold on;plot(dff_true(i,1:1000)');      
            subplot(5,2,10);hold on;plot(dff_psth_true(i,:)');
        end

    end
end


end