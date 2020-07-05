global nPlanes
nTones = 17;
nTrials = 10;
nFramesPerTone = 100/nPlanes; % 50
nFramesPerTrial = nFramesPerTone * nTones; % 850
startTrial = 2; % the first tone is on frame 0
nFrames = nFramesPerTrial*nTrials; % 4250
frameRate = 30.98/nPlanes; % in the future, do not hard code this.
pretoneFrames = 10;
baselineFrames = 5;

smoothWindow = 5;

saveSingleNeuronFlag = true;
toneOnset = 0/nPlanes;