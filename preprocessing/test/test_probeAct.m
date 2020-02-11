%% hre is for cd036
%beh1v3 = importdata('cd036_1v3.txt');
beh1v4 = importdata('cd036_1v4.txt');
beh2v2 = importdata('cd036_2v2.txt');

%frameSum1v2 = 8500*2+2500+13193; %+11469+10289+5000+13408;
%frameSum2v2 = 8500*2+2500+13193+11469+10289+5000+13408;
%frameSum1v3 = 2500*2 + 8500*2 + 16774+9663;
frameSum1v4 = 2500*2 + 8500*2 + 16774+9663 + 11745+5000;
frameSum2v2 = 2500*2 + 8500*2 + 16774+9663 + 11745+5000 + 14834+13185;

%toneFrame = round(beh1v3(:,12)/2) + frameSum1v3;
%plotStuffReinfProbe(toneFrame,tempTC,beh1v3);

toneFrame = round(beh1v4(:,12)/2) + frameSum1v4;
plotStuffReinfProbe(toneFrame,tempTC,beh1v4);

toneFrame = round(beh2v2(:,12)/2) + frameSum2v2;
plotStuffReinfProbe(toneFrame,tempTC,beh2v2);


%% here is for cd017

%h5list = {'cd017_000_000','cd017_000_001','cd017_001_000','cd017_001_001',...
%'cd017_001_002','cd017_001_003','cd017_002_000','cd017_002_001','cd017_002_002',...
%'cd017_002_003','cd017_002_005','cd017_003_001','cd017_003_002','cd017_003_003',...
%'cd017_003_004','cd017_003_005','cd017_004_001','cd017_004_002'};
%for i = 1:length(h5list); sbxread(h5list{i},1,1); allFrames(i) = round(info.max_idx/2); end
%allFramesCumsum = cumsum(allFrames);

beh2v3 = importdata('cd017_2v3.txt');
beh4v2 = importdata('cd017_4v2.txt');

frameSum2v3 = 69780;
frameSum4v2 = 128624;


toneFrame = round(beh2v3(:,12)/2) + frameSum2v3;
plotStuffReinfProbe(toneFrame,tempTC,beh2v3);

toneFrame = round(beh4v2(:,12)/2) + frameSum4v2;
plotStuffReinfProbe(toneFrame,tempTC,beh4v2);


function  toneAct = plotStuffReinfProbe(toneFrame,tempTC,beh)
neuron = size(tempTC,1);
trials = size(beh,1);
toneAct = zeros(neuron,60,trials);
for i = 1:trials; toneAct(:,:,i) = ...
        tempTC(:,toneFrame(i)-9:toneFrame(i)+50); end
meanFlu = mean(reshape(toneAct,neuron,[]),2);
stdFlu = std(reshape(toneAct,neuron,[]),0,2);
toneZ = (toneAct-repmat(meanFlu,1,60,trials)) ./ repmat(stdFlu,1,60,trials);
toneDff = (toneAct-repmat(meanFlu,1,60,trials)) ./ repmat(meanFlu,1,60,trials);
%figure; plot(mean(mean(toneDff,3),1))

toneDffProbe = toneDff(:,:,~logical(beh(:,13)));
toneDffReinf = toneDff(:,:,logical(beh(:,13)));
figure; plot(mean(mean(toneDffProbe,3),1));
hold on ; plot(mean(mean(toneDffReinf,3),1));
legend('probe','reinf')
xlabel('frames'); ylabel('avg dff of all neuron')
figure; subplot(2,1,1); title('probe')
imagesc(mean(toneDffProbe,3)); colorbar; caxis([-0.05 0.1])
subplot(2,1,2); title('reinf')
imagesc(mean(toneDffReinf,3)); colorbar; caxis([-0.05 0.1])
end

