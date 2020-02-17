global info;
a = sbxread('cd049_000_003',1,1);
maxidx = info.max_idx;

mov = sbxread('cd049_000_003',1,maxidx);
mov = squeeze(mov);
avg_mov = mean(mov,3);
mov = reshape(mov,512*796,[]);
%%
figure; imagesc(avg_mov); colormap gray; a0 = imfreehand();
a0 = a0.createMask();
a0mask = reshape(a0,1,[]);
a0TC = mean(mov(a0mask,:));
onoffset = diff(a0TC);
onset = onoffset < -2000;
onsetFrames = find(onset);

onsetFlag = [0 0 0 1 1 1 0 0 0 2 2 2 0 0 0 0 0 0 3 3 3 ...
    0 0 0 4 4 4 0 0 0];

figure; imagesc(avg_mov); colormap gray; a1 = imfreehand();
a1 = a1.createMask();
a1mask = reshape(a1,1,[]);
a1TC = mean(mov(a2mask,:));

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

%% for cd049_000_003

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
for i = 1:7
plot(mean(onsetTC{i},2),'color',1-[0.1 0.1 0.1]*i);hold on; title('all scale')
end