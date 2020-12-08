%% load data
clear;
%parobj = parpool;
global tempTC ishere nFrames; 
load('C:\Users\zzhu34\Documents\tempdata\excitData\cd036\TC\cd036_TC_plane0.mat','tempTC');
load('C:\Users\zzhu34\Documents\tempdata\excitData\cd036\roi\ishere_plane0.mat','ishere');
ishere = ishere{1};
[nFrames, nFrames_oneplane] = func_readnFrames('cd036','root',...
    'C:\Users\zzhu34\Documents\tempdata\excitData\config\mouse');
nFrames = [0 0; nFrames];

%% Parse Dff into each session
for i = 1:size(nFrames,1)-1
    sessionTC = tempTC(:,nFrames(i)+1:nFrames(i+1),1);
    tempSessionDff = (sessionTC - repmat(median(sessionTC,2),1,size(sessionTC,2))) ./ repmat(median(sessionTC,2),1,size(sessionTC,2));
    tempSessionDff = tempSessionDff - smoothdata(tempSessionDff,2,'movmedian',3000);
    sessionDff{i} = tempSessionDff(:,2:end);
end

%% 4 - test the calmAn AR2 methods
datapath = 'C:\Users\zzhu34\Documents\tempdata\deconv_test_cd036\allSessions\';
[c,s,p] = deconvolveDataCalmAn(sessionDff,'foopsi',0.95,false);
save([datapath 'cd036_calman_ar2_foo95_optb_nosmin_c.mat'],'c','-v7.3');
save([datapath 'cd036_calman_ar2_foo95_optb_nosmin_s.mat'],'s','-v7.3');
save([datapath 'cd036_calman_ar2_foo95_optb_nosmin_p.mat'],'p');

[c,s,p] = deconvolveDataCalmAn(sessionDff,'foopsi',0.95,true);
save([datapath 'cd036_calman_ar2_foo95_optb_c.mat'],'c','-v7.3');
save([datapath 'cd036_calman_ar2_foo95_optb_s.mat'],'s','-v7.3');
save([datapath 'cd036_calman_ar2_foo95_optb_p.mat'],'p');

%% stop parallel pooling

% delete(parobj);
%% all functions


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
