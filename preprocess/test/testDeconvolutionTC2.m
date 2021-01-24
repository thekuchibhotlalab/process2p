%% load data
clear;
mouse = 'cd017';
load(['C:\Users\zzhu34\Documents\tempdata\excitData\' mouse '\TC\' mouse '_TC_plane0.mat'],'tempTC');
configPath = 'C:\Users\zzhu34\Documents\tempdata\excitData\config\mouse\';
datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mouse '\ar1_eachday_smooth3\'];
mkdir(datapath);
deconvolveByDay(tempTC,mouse,configPath,datapath);
%% load data
clear;
mouse = 'cd036';
load(['C:\Users\zzhu34\Documents\tempdata\excitData\' mouse '\TC\' mouse '_TC_plane0.mat'],'tempTC');
configPath = 'C:\Users\zzhu34\Documents\tempdata\excitData\config\mouse\';
datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mouse '\ar1_eachday_smooth3\'];
mkdir(datapath);
deconvolveByDay(tempTC,mouse,configPath,datapath);
%% load data
clear;
mouse = 'cd037';
load(['C:\Users\zzhu34\Documents\tempdata\excitData\' mouse '\TC\' mouse '_TC_plane0_final.mat'],'TC');
configPath = 'C:\Users\zzhu34\Documents\tempdata\excitData\config\mouse\';
datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mouse '\ar1_eachday_smooth3\'];
mkdir(datapath);
deconvolveByDay(TC,mouse,configPath,datapath); %deconvolveMaster(TC,nFrames,mouse,datapath)

%% load data
clear;
mouse = 'cd041';
load(['C:\Users\zzhu34\Documents\tempdata\excitData\' mouse '\TC\' mouse '_TC_plane0.mat'],'tempTC');
configPath = 'C:\Users\zzhu34\Documents\tempdata\excitData\config\mouse\';
datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mouse '\ar1_eachday_smooth3\'];
mkdir(datapath);
deconvolveByDay(tempTC,mouse,configPath,datapath);

%% load data
clear;
mouse = 'cd042';
load(['C:\Users\zzhu34\Documents\tempdata\excitData\' mouse '\TC\' mouse '_TC_plane0_final.mat'],'TC');
configPath = 'C:\Users\zzhu34\Documents\tempdata\excitData\config\mouse\';
datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mouse '\ar1_eachday_smooth3\'];
mkdir(datapath);
deconvolveByDay(TC,mouse,configPath,datapath);

%% load data
clear;
mouse = 'cd044';
load(['C:\Users\zzhu34\Documents\tempdata\excitData\' mouse '\TC\' mouse '_TC_plane0_final.mat'],'TC');
configPath = 'C:\Users\zzhu34\Documents\tempdata\excitData\config\mouse\';
datapath = ['C:\Users\zzhu34\Documents\tempdata\deconv_test_' mouse '\ar1_eachday_smooth3\'];
mkdir(datapath);
deconvolveByDay(TC,mouse,configPath,datapath);


%% all functions
function deconvolveMaster(TC,nFrames,mouse,datapath)
tic;
% Parse Dff into each session
for i = 1:size(nFrames,1)-1
    sessionTC = TC(:,nFrames(i)+1:nFrames(i+1),1);
    tempSessionDff = (sessionTC - repmat(median(sessionTC,2),1,size(sessionTC,2))) ./ repmat(median(sessionTC,2),1,size(sessionTC,2));
    tempSessionDff = tempSessionDff - smoothdata(tempSessionDff,2,'movmedian',3000);
    tempSessionDff(isnan(tempSessionDff)) = 0;
    sessionDff{i} = tempSessionDff(:,2:end);
end
% use the calmAn AR2 methods
[s,c,p] = deconvolveDff(sessionDff);
save([datapath '\' mouse '_calman_ar2_foo95_optb_nosmin_c.mat'],'c','-v7.3');
save([datapath '\' mouse '_calman_ar2_foo95_optb_nosmin_s.mat'],'s','-v7.3');
save([datapath '\' mouse '_calman_ar2_foo95_optb_nosmin_p.mat'],'p');
toc;
end

function deconvolveAllDay(TC,mouse,configpath,datapath,nPlane)
tic;
% Parse Dff into each session
[nFrames, ~] = func_readnFrames(mouse,'root',configpath);
nFrames = [0 0; nFrames];
configTable =  readtable([configpath '\' mouse '_config.csv']);
sessionDff = cell(1,size(nFrames,1)-1);
% Get Dff of each session
for i = 1:size(nFrames,1)-1
    sessionTC = TC(:,nFrames(i,nPlane)+1:nFrames(i+1,nPlane));
    tempSessionDff = (sessionTC - repmat(median(sessionTC,2),1,size(sessionTC,2))) ./ repmat(median(sessionTC,2),1,size(sessionTC,2));
    tempSessionDff = tempSessionDff - smoothdata(tempSessionDff,2,'movmedian',3000);
    tempSessionDff(isnan(tempSessionDff)) = 0;
    sessionDff{i} = tempSessionDff;
    if nPlane == 1; sessionDff{i} =  sessionDff{i}(:,2:end); end
end
% Estimate Decay Parameter for all days
sessionTCAll = zeros(size(TC)); pars = [];
for i = 1:length(sessionDff)
    sessionTCAll(:,nFrames(i,nPlane)+(1:size(sessionDff{i},2)) ) = sessionDff{i};
end
    
% estimate time constant for each neuron separately
for j = 1:size(sessionTCAll,1)
    y = reshape(sessionTCAll(j,:), [], 1); 
    pars(j,:) = estimate_time_constant(y, 1, GetSn(y));
end

[temps,tempc,tempp] = deconvolveDff(sessionDff,'pars',pars);
s = temps; c = tempc; p = tempp;

save([datapath '\' mouse '_calman_ar1_foo90_pars_allday_c.mat'],'c','-v7.3');
save([datapath '\' mouse '_calman_ar1_foo90_pars_allday_s.mat'],'s','-v7.3');
save([datapath '\' mouse '_calman_ar1_foo90_pars_allday_p.mat'],'p','pars');
toc;
end

function deconvolveByDay(TC,mouse,configpath,datapath)
tic;
% Parse Dff into each session
[nFrames, ~] = func_readnFrames(mouse,'root',configpath);
nFrames = [0 0; nFrames];
configTable =  readtable([configpath '\' mouse '_config.csv']);
% Get Dff of each session
for i = 1:size(nFrames,1)-1
    sessionTC = TC(:,nFrames(i,1)+1:nFrames(i+1,1));
    tempSessionDff = (sessionTC - repmat(median(sessionTC,2),1,size(sessionTC,2))) ./ repmat(median(sessionTC,2),1,size(sessionTC,2));
    tempSessionDff = tempSessionDff - smoothdata(tempSessionDff,2,'movmedian',3000);
    tempSessionDff(isnan(tempSessionDff)) = 0;
    sessionDff{i} = tempSessionDff(:,2:end);
end
% Estimate Decay Parameter for each day,a nd 
nDays = sort(unique(configTable.Day),'ascend');
s = cell(1,size(configTable,1));c = cell(1,size(configTable,1));p = cell(1,size(configTable,1));
count = 0;
for i = nDays'
    count = count + 1;
    disp(['Day ' int2str(i) ' started!'])
    sessionNum = find(configTable.Day == i);
    sessionTC = []; tempPars = [];
    for j = 1:length(sessionNum)
        sessionTC = cat(2, sessionTC, sessionDff{sessionNum(j)});
    end
    
    % estimate time constant for each neuron separately
    for j = 1:size(sessionTC,1)
        y = reshape(sessionTC(j,:), [], 1); 
        y = smoothdata(y,'movmean',3);
        tempPars(j,:) = estimate_time_constant(y, 1, GetSn(y));
    end
    pars{count} = tempPars;
    sessionTCCell = sessionDff(sessionNum);
    [temps,tempc,tempp] = deconvolveDff(sessionTCCell,'pars',tempPars);
    s(sessionNum) = temps; c(sessionNum) = tempc; p(sessionNum) = tempp;
end

save([datapath '\' mouse '_calman_ar2_foo90_pars_day_c.mat'],'c','-v7.3');
save([datapath '\' mouse '_calman_ar2_foo90_pars_day_s.mat'],'s','-v7.3');
save([datapath '\' mouse '_calman_ar2_foo90_pars_day_p.mat'],'p','pars');
toc;
end

% BELOW ARE OLD DEPRECATE FUNCTIONS
function [c,s,p] = deconvolveData(selectDff, params)
nNeuron = size(selectDff{1},1);
c = {};
s = {};
p = {};
tic;
for i=1:length(selectDff)
    c{i} = zeros(nNeuron, size(selectDff{i},2));
    s{i} = zeros(nNeuron, size(selectDff{i},2));
   for j = 1:nNeuron
        [c_deconv, s_deconv, options] = deconvolveCa(selectDff{i}(j,:),params{:});% time x one cell
        c{i}(j,:) = c_deconv;
        s{i}(j,:) = s_deconv;
        p{i,j} = options;
   end  
end
t = toc;
disp(['Time Elapsed = ' num2str(t,'%.2f') ' secs.']);
end

function [c,s,p] = deconvolveDataCalmAn(selectDff,method,lam_pr,est_smin)
nNeuron = size(selectDff{1},1);
c = {};
s = {};
p = {};
options.fr = 15.63;options.decay_time = 0.7;
options.spk_SNR = 0.5; options.lam_pr = lam_pr;
model_ar = 'ar2';

for i=1:length(selectDff)
    tic;
    disp(['Session' int2str(i) 'started'])
    c{i} = zeros(nNeuron, size(selectDff{i},2));
    s{i} = zeros(nNeuron, size(selectDff{i},2));
    for j = 1:nNeuron
        tempDff = selectDff{i}(j,:);
        if est_smin; spkmin = options.spk_SNR*GetSn(tempDff); else; spkmin = 0; end 
        if strcmp(method,'foopsi') || strcmp(method,'thresholded')
            lam = choose_lambda(exp(-1/(options.fr*options.decay_time)),GetSn(tempDff),options.lam_pr);
            [c_deconv,s_deconv,ops] = deconvolveCa(tempDff,model_ar,'method',method,'optimize_pars',true,...
                'optimize_b',true,'maxIter',20,'window',150,'lambda',lam,'smin',spkmin);
        else 
            [c_deconv,s_deconv,ops] = deconvolveCa(tempDff,model_ar,'method',method,'optimize_pars',true,...
                'optimize_b',true,'maxIter',20,'window',150,'smin',spkmin);
        end
        c{i}(j,:) = c_deconv;
        s{i}(j,:) = s_deconv;
        p{i,j} = ops;
    end  
    t = toc;
    disp(['Time Elapsed = ' num2str(t,'%.2f') ' secs.']);
end


end
