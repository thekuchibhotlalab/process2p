function [CaTrace, allParam] = func_getFakeCaTrace(spikeITI,noise)
global convFilter nFrames nTrial nFramesPerTrial spikeDelay

nTrial = 200;
nFramesPerTrial = 200;
nFrames = nTrial * nFramesPerTrial;
spikeDelay = 5;

%evokedSpike  = [3 9 15];
%sponSpike = [1 3 5];
%spikeITI = [10 20 30];
%noise = [0.5 1.0 1.5];
%baseline = [300 600 900];
%fluroGain = [100 300 500];

evokedSpike  = 9;
sponSpike = 3;
%spikeITI = [50 80];
%noise = [0.1 0.2 0.3 0.4];
baseline = 1200;
fluroGain = 500;
allParam = combvec(evokedSpike,sponSpike,spikeITI,noise,baseline,fluroGain);
nNeuron = size(allParam,2);

decayConstant = 0.933; % for 10 frames half life
convFilter = zeros(1,nFramesPerTrial);
convFilter(1) = 1;
for i = 2:nFramesPerTrial
    convFilter(i) =  convFilter(i-1) * decayConstant;
end

CaTrace = zeros(nNeuron,nFrames);
for i = 1:nNeuron
    tempCell = num2cell(allParam(:,i));
    tempTrace = getFakeCaTrace(tempCell{:});
    CaTrace(i,:) = tempTrace(1:nFrames);
end


end

function CaTrace = getFakeCaTrace(evokedSpike,sponSpike,spikeITI,noise,baseline,fluroGain)
    global convFilter nFrames nTrial nFramesPerTrial spikeDelay
    spikeTrain = zeros(nFramesPerTrial, nTrial);
    spikeTrain(spikeDelay,:) = evokedSpike;
    spikeTrain = reshape(spikeTrain,nFrames,1);
    sponSpikeITI = exprnd(spikeITI,10000,1);
    sponSpikeTime = ceil(cumsum(sponSpikeITI));
    sponSpikeTime = (sponSpikeTime(sponSpikeTime<nFrames));

    spikeTrain(sponSpikeTime) = spikeTrain(sponSpikeTime) + sponSpike;
    CaTrace = conv(spikeTrain,convFilter);
    CaTrace = CaTrace * fluroGain + baseline + baseline * normrnd(0,noise,size(CaTrace,1),size(CaTrace,2)) ;
    %CaTrace = baseline + baseline * normrnd(0,noise,size(CaTrace,1),size(CaTrace,2)) ;
end