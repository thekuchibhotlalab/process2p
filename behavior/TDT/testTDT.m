% new code (revised resetData) started from 3v2
% 3v2-4 TTL working, aligned well
% 3v4 has one water delivery that is off by one frame. 
% note 3v5 is grabbed manually. others are grabbed by TTL
% weird thing is that when finished grabbing 3v5 manually, if do TTL, then
% TTL started without the session being started. 
% 3v5 not aligned, with a shift of 314
% 3v6 is grabbed with TTL, aligned well
% 3v7 is grabbed with TTL, not aligned, constant shift for only first few
% licks
% 3v8 grabbed manually. TTL did not work, but aligned greatly! However,
% when TTL is clicked after aborting the session, acquisition started. 
% 3v9 grabbed with TTL, not aligned, shift of -272
% 3v10 grabbed with TTL, aligned
% changed other duplicated variables, started again
% 4v1 grabbed with TTL, aligned. TTL was not starting, but started after
% unclicking it, and when clicked back on, it started scanning immediately
% after scanning finished
% fixed TTL trigger using zBusTrig
% 4v4 has a shift of 197 (using zBus.zBusTrigB(0,0,2);)
% 4v5 has a shift of -1 (using zBus.zBusTrigB(0,0,5);)
% 4v6 has a shift of 0 (using zBus.zBusTrigB(0,0,5);)
% 4v7 has a shift of -702 (using zBus.zBusTrigB(0,0,5);)
% moved ResetDataFrame etc. before start of scanning
% 4v8 has a shift of 0 (using zBus.zBusTrigB(0,0,5);)
% 4v9 has a shift of 0 (using zBus.zBusTrigB(0,0,5);)
% 4v10 has a shift of 354 (using zBus.zBusTrigB(0,0,5);)
% added lines to output frameindex (index that scanning start at)
% 4v11 has a shift of 0 (with same frame, lick, water index)
% 4v12 has a shift of 0 (with same frame, lick, water(+1) index)
% 4v13 has a shift of 0 (with same frame, lick, water index)
% 4v14 has a shift of 105 (with same frame, lick, water index)
% 4v15 has a shift of 0 (all trials good)
% 4v16 has a shift of 0 (all trials good)
% 4v17 has a shift of 0 (all trials good)
% 4v18 has a shift of 0 for water, but ~1000 index shift for lick. lick was
% not aligned good. hypo: resetting is the culprit
% remove resetting!
% 4v19 has a shift of 0 (all trials good)
% 4v20 has a shift of 0 (all trials good)
% 4v21 has a shift of 0 (all trials good)
% 4v22 has a shift of 0 (all trials good)
% 4v23 has a shift of 0 (all trials good)
% 4v24 has a shift of 0 (all trials good)
% 4v25 has a shift of 0 (all trials good)
% 4v26 has a shift of 0 (all trials good)
% 4v27 has a shift of 0 (all trials good)
% 4v28 has a shift of 0 (all trials good)
% 4v29 has a shift of 0 (all trials good) Long session (100 trials). Great!
% 11/21 revised unpair control code. test it!
% 5v1 good
% 5v2 good
% 5v3 good
% TTL only work when shutter is on
% manually focusing and stopping gives multiple ramp up and downs
% 5v4 good (except extra probe water delivery)
% 5v5 good
clear;
framerate = 31.25;
temp = 'test_1v1';
frames=load([temp '.txtframes.txt']);
lick=load([temp '.txtlicks.txt']);
water=load([temp '.txtwaterdelivery.txt']);

frameDiff=diff(frames);
frameDiff=frameDiff>0;
frameNum = cumsum(frameDiff);

lickDiff = diff(lick);
lickDiff = lickDiff>0;
lickFrame = frameNum(lickDiff)';

waterDiff = diff(water);
waterDiff = waterDiff>0;
waterFrame = frameNum(waterDiff)';

beh = load([temp '.txt']);

temp = beh(:,10);
temp(temp==1000000) = [];
%disp(waterFrame - unique(temp))

[~,maxidx] = max(frameNum);

figure;
plot(frameNum, 'LineWidth',2)
hold on;
plot(waterDiff * max(frameNum), 'Color', [0 0 0], 'LineWidth',2)
plot(lickDiff * max(frameNum), 'Color', [0.7 0.7 0.7])
xlim([0 maxidx + 10000])

%%

disp([ int2str(sum(~isnan(beh(:,3)))) ' tones out of ' int2str((beh(end,2))) 'trials'] )
%tempToneRewardOffset = (beh(:,12) - beh(:,16))./framerate;
%tempToneRewardOffset = tempToneRewardOffset(~isnan(tempToneRewardOffset));
%histogram(tempToneRewardOffset)
rewardedTrials = beh(:,10)<100000 & beh(:,13)==1;
disp(['tones in rewarded trials: ' int2str(sum(~isnan(beh(rewardedTrials,12))))])

figure;
plot(beh(:,12),ones(1,size(beh,1)),'.'); hold on;
plot(beh(rewardedTrials,10),ones(1,sum(rewardedTrials)),'.');

rewardTime = beh(rewardedTrials,10);
for i = 1:length(rewardTime)
    minvalue(i) = min(abs(rewardTime(i) - beh(:,12)));
end