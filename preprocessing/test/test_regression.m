clear all; close all
stimB  = 0.3; % 30 HZ
lickB = 0.1; % 5HZ
movB = 0.4; % 3HZ
stateB = 0.3; % 3 HZ

lickITImean = 40;
lickITImin = 15;
nSamplePerTrial = 300;
nTrial = 100;
nFrames = nSamplePerTrial * nTrial;

stimEvent = zeros(nSamplePerTrial, nTrial);
stimEvent(101:120,:) = 0.8;
stimEvent = reshape(stimEvent,nSamplePerTrial*nTrial,1);

lickEvent = zeros(nSamplePerTrial*nTrial,1);
moveEvent = randn(nSamplePerTrial*nTrial,1);
moveEvent = smoothdata(moveEvent,'gaussian',30);
stateEvent = randn(nSamplePerTrial*nTrial,1);
stateEvent = smoothdata(stateEvent,'gaussian',100);

lickITI = lickITImin + exprnd(lickITImean-lickITImin,10000,1);
lickTime = round(cumsum(lickITI));
lickTime = (lickTime(lickTime<nFrames));
lickEvent(lickTime) = 1;
for i = 1:5
    lickEvent(lickTime+i) = 1;
end

timeAxis = 0:1:10;
lickKernel = normpdf(timeAxis,3,3);
lickAct = conv(lickEvent, lickKernel);
lickAct = lickAct(1:nFrames);

activity = (stimB * stimEvent + lickB * lickAct + moveEvent * movB + stateEvent * stateB + randn(nFrames,1)*0.03 )*100;
%activity = smoothdata(activity,'movmean',10);
%randnumgen = rand(nFrames,1);
%spikes = double(randnumgen < (stim*0.2) );
regressor = [stimEvent, lickEvent, moveEvent];
resultsLM = fitlm(regressor,activity);
resultsRidge = ridge(activity,regressor,0.001,0);

plotTime = 2000;
figure;subplot(4,1,1);plot(smoothdata(activity(1:plotTime),'gaussian',10),'LineWidth',3, 'Color', [0 0.4470 0.7410]); axis off;
subplot(4,1,2);plot(stimEvent(1:plotTime),'LineWidth',3, 'Color', [0.8500 0.3250 0.0980]); axis off;
subplot(4,1,3);plot(lickEvent(1:plotTime),'LineWidth',3, 'Color', [0.9290 0.6940 0.1250]); axis off; ylim([0 2])
subplot(4,1,4);plot(moveEvent(1:plotTime),'LineWidth',3, 'Color', [0.4940 0.1840 0.5560]); axis off; ylim([-1.5 1.5])
figure;
subplot(2,1,1)
plot(smoothdata(activity(1:plotTime),'gaussian',10),'LineWidth',2, 'Color', [0 0.4470 0.7410]); axis off; hold on
plot(smoothdata(resultsLM.Fitted(1:plotTime),'gaussian',10),'LineWidth',2, 'Color', [0.8500 0.3250 0.0980]);

regressor = [stimEvent, lickAct, moveEvent];
resultsLM2 = fitlm(regressor,activity);
resultsRidge2 = ridge(activity,regressor,0.001,0);

subplot(2,1,2) 
plot(activity(1:plotTime),'LineWidth',2, 'Color', [0 0.4470 0.7410]); axis off; hold on
plot(resultsLM2.Fitted(1:plotTime),'LineWidth',2, 'Color', [0.8500 0.3250 0.0980]);

newAct = exp(activity/40*3.7);
figure; plot(activity(1:plotTime)); hold on; plot(newAct(1:plotTime))
regressor = [stimEvent, lickAct, moveEvent];
resultsLM3 = fitlm(regressor,newAct);
resultsRidge3 = ridge(newAct,regressor,0.001,0);
resultGLM = fitglm(regressor, newAct,'Link','log','Distribution','poisson');

figure;
subplot(2,1,1)
plot(newAct(1:plotTime),'LineWidth',2, 'Color', [0 0.4470 0.7410]); axis off; hold on
plot(smoothdata(resultsLM3.Fitted(1:plotTime),'gaussian',10),'LineWidth',2, 'Color', [0.8500 0.3250 0.0980]);
subplot(2,1,2)
plot(newAct(1:plotTime),'LineWidth',2, 'Color', [0 0.4470 0.7410]); axis off; hold on
plot(smoothdata(resultGLM.Fitted.Response(1:plotTime),'gaussian',10),'LineWidth',2, 'Color', [0.8500 0.3250 0.0980]);

%%

a = mvnrnd([1,1],[0.05 0.02; 0.02, 0.05],4); a(5,:) = [1.2 2.5];
b = mvnrnd([1,1],[0.05 0.02; 0.02, 0.05],5);

%%

figure; scatter(a(:,1),a(:,2),'filled'); hold on; scatter(b(:,1),b(:,2), [], [0.5 0.5 0.5],'filled');

%%

c = mvnrnd([1,1],[0.05 0.035; 0.035, 0.05],50);
figure; scatter(c(:,1),c(:,2)+0.3,'filled');
xlim([0.4 1.6])
ylim([0.7 1.9])