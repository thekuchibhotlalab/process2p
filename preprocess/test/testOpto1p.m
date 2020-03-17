datapath = 'Q:\celine\cd049';
cd(datapath)
global info;
a = sbxread('cd049_000_002',1,1);
maxidx = info.max_idx;

mov = sbxread('cd049_000_002',1,maxidx);
mov = squeeze(mov);
avg_mov = mean(mov,3);
mov = reshape(mov,512*796,[]);
%%
figure; imagesc(avg_mov); colormap gray; a0 = imfreehand();
a0 = a0.createMask();
a0mask = reshape(a0,1,[]);
a0TC = mean(mov(a0mask,:));
onoffset = diff(a0TC);
onset = onoffset < -1000;
onsetFrames = find(onset);

figure; imagesc(avg_mov); colormap gray; a1 = imfreehand();
a1 = a1.createMask();
a1mask = reshape(a1,1,[]);
a1TC = mean(mov(a1mask,:));

figure; imagesc(avg_mov); colormap gray; a2 = imfreehand();
a2 = a2.createMask();
a2mask = reshape(a2,1,[]);
a2TC = mean(mov(a2mask,:));

selectFrames = -20:100;
%% for cd049_000_002
onsetFlag = [0 0 -1 1 1 1 0 0 0 2 2 2 0 0 3 3 3 ...
    0 0 0 4 4 4 0 0 0 5 5 5 0 0 0 6 6 6 0 0 0];

temp = 0:6;
onsetTC = {};
for j = 1:7
for i = 1:length(onsetFrames(onsetFlag==temp(j)))
tempFrames = onsetFrames(onsetFlag==temp(j));
onsetTC{j}(:,i) = a2TC(tempFrames(i)+selectFrames) - mean(a1TC(tempFrames(i)+(-20:0)));
end
end
figure;for i = 1:7
subplot(2,4,i)
plot(onsetTC{i}); if i==1; title('no light');else title(['scale' num2str(i-1)]);end
end
subplot(2,4,8);
for i = 1:7
plot(mean(onsetTC{i},2),'color',1-[0.1 0.1 0.1]*i);hold on; title('all scale')
end
%% plot cd049_000_002 in better way
allTC = [];
for i = 1:7
    temp = (smoothdata(onsetTC{i},1,'gaussian',8) - mean(mean(onsetTC{i}([1:20],:),2),1)) / mean(mean(onsetTC{i}([1:20],:),2),1);
    allTC(:,i) = median(temp,2);
    allTCSem(:,i) = std(onsetTC{i},0,2) / sqrt(size(onsetTC{i},2));
    allToneOnset{i} = mean(temp(21:35,:),1);
    meanToneOnset(i) = mean(allToneOnset{i});
    semToneOnset(i) = std(allToneOnset{i})./sqrt(length(allToneOnset{i}));
end
figure; hold on;
for i = 1:7
    if i==1
        plot(0:120,allTC(:,i),'color',[0.3 0.3 0.3],'LineWidth',2)
    elseif i==2
        plot(0:120,allTC(:,i),'color',[0.3010    0.7450    0.9330],'LineWidth',2)
    elseif i==4
        plot(0:120,allTC(:,i),'color',[0 0.4470 0.7410],'LineWidth',2)
    end
    
     
end
plot([19 19],[-0.6 0.6],'color',[0.8 0.8 0.8]);
plot([0 120],[0 0],'color',[0.8 0.8 0.8]);
xlim([0 65]);ylim([-0.6 0.6])
xticks([5 20 35 50 65]-1); xticklabels({'-1','0','1','2','3'})
xlabel('Time after tone onset (s)')
ylabel('Df/f')
legend('No Stimulation','Stimulation Power=1','Stimulation Power=2')

figure; bar(0:6,meanToneOnset); hold on; errorbar(0:6,meanToneOnset,semToneOnset)
tempMean = [meanToneOnset([1 2 4]) mean([allToneOnset{6} allToneOnset{7}])];
tempStd = [semToneOnset([1 2 4]) std([allToneOnset{6} allToneOnset{7}])/sqrt(6)];
figure; bar(0:3,tempMean,'FaceColor',[0.3010 0.7450 0.9330]); hold on; 
errorbar(0:3,tempMean,tempStd,'.','LineWidth',0.5,'Color',[0.00 0.0 0.0])
xticklabels({'0','1','2','3'})
xlabel('Stimulation Power'); ylabel('Mean tone-evoked df/f')
ylim([-1 1])



%% for cd049_000_003

onsetFlag = [0 0 0 1 1 1 0 0 0 2 2 2 0 0 0 0 0 0 3 3 3 ...
    0 0 0 4 4 4 0 0 0];

temp = 0:4;
onsetTC = {};
for j = 1:5
for i = 1:length(onsetFrames(onsetFlag==temp(j)))
tempFrames = onsetFrames(onsetFlag==temp(j));
onsetTC{j}(:,i) = a2TC(tempFrames(i)+selectFrames) - mean(a1TC(tempFrames(i)+(-20:0)));
end
end
figure;for i = 1:5
subplot(2,3,i)
plot(onsetTC{i}); if i==1; title('no light');else title(['scale' num2str(i-1)]);end
end
subplot(2,3,6);
for i = 1:5
plot(mean(onsetTC{i},2),'color',1-[0.1 0.1 0.1]*i);hold on; title('all scale')
end