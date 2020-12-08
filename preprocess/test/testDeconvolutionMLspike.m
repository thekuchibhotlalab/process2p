%% load data
clear;
parobj = parpool;
global tempTC ishere nFrames; 
load('C:\Users\zzhu34\Documents\tempdata\excitData\cd036\TC\cd036_TC_plane0.mat','tempTC');
load('C:\Users\zzhu34\Documents\tempdata\excitData\cd036\roi\ishere_plane0.mat','ishere');
ishere = ishere{1};
[nFrames, nFrames_oneplane] = func_readnFrames('cd036','root',...
    'C:\Users\zzhu34\Documents\tempdata\excitData\config\mouse');
nFrames = [0 0; nFrames];
%% 1 - select the sessions (baseline)
tbpsessionSelected = [13 18 22 26 30 34 38]; 
day = 3:9;
[selectDff, selectNeuron] = selectSession(sessionSelected,day);
disp('Selection Done')

% confoo - 1. estimate baseline, 2. give baseline
params = {'constrained','ar1','fit_params', 'optimize_b'};[c,s,p] = deconvolveData(selectDff, params);%params = {'constrained','ar1', 0.9381, 'optimize_b'};[c,s,p] = deconvolveData(selectDff, params);
save('baseline\cd036_confoo.mat','c','s','p','day','selectNeuron','params');
params = {'constrained','ar1', 'fit_params', 'b',0};[c,s,p] = deconvolveData(selectDff, params);
save('baseline\cd036_confoo_b.mat','c','s','p','day','selectNeuron','params');
params = {'constrained','ar1','fit_params', 'optimize_b','s_min',-3};[c,s,p] = deconvolveData(selectDff, params);%params = {'constrained','ar1', 0.9381, 'optimize_b'};[c,s,p] = deconvolveData(selectDff, params);
save('baseline\cd036_confoo_smin-3.mat','c','s','p','day','selectNeuron','params');
%params = {'constrained','ar2', 'fit_params', 'optimize_b'};[c,s,p] = deconvolveData(selectDff, params);
%save('baseline\cd036_confoo_ar2.mat','c','s','p','day','selectNeuron','params');
%params = {'constrained','ar2', 'fit_params', 'optimize_b','s_min',-3};[c,s,p] = deconvolveData(selectDff, params);
%save('baseline\cd036_confoo_ar2_smin-3.mat','c','s','p','day','selectNeuron','params');

% foo - 1. estimate baseline, 2. give baseline, 3. different lambda = 0.0,
% 4. diff lambda = 0.5, give baseline, 5,6. different lambda = 1.0
params = {'foopsi','ar1', 'fit_params', 'optimize_b','lambda', 0};[c,s,p] = deconvolveData(selectDff, params);
save('baseline\cd036_foo_lambda0.mat','c','s','p','day','selectNeuron','params');
params = {'foopsi','ar1', 'fit_params', 'b',0,'lambda', 0};[c,s,p] = deconvolveData(selectDff, params);
save('baseline\cd036_foo_lambda0_b.mat','c','s','p','day','selectNeuron','params');

params = {'foopsi','ar1', 'fit_params', 'optimize_b','lambda', 0.5};[c,s,p] = deconvolveData(selectDff, params);
save('baseline\cd036_foo_lambda05.mat','c','s','p','day','selectNeuron','params');
params = {'foopsi','ar1', 'fit_params', 'b',0,'lambda', 0.5};[c,s,p] = deconvolveData(selectDff, params);
save('baseline\cd036_foo_lambda05_b.mat','c','s','p','day','selectNeuron','params');

params = {'foopsi','ar1', 'fit_params', 'optimize_b','lambda', 1.0};[c,s,p] = deconvolveData(selectDff, params);
save('baseline\cd036_foo_lambda10.mat','c','s','p','day','selectNeuron','params');
params = {'foopsi','ar1', 'fit_params', 'b',0,'lambda', 1.0};[c,s,p] = deconvolveData(selectDff, params);
save('baseline\cd036_foo_lambda10_b.mat','c','s','p','day','selectNeuron','params');

% thresholded - 1. optimize b, threshold = 0.99, 2. thre = 0.95, 3,4 . smin =
% -3
params = {'thresholded','ar1','fit_params','optimize_smin', 'optimize_b','thresh_factor', 0.99};
[c,s,p] = deconvolveData(selectDff, params);
save('baseline\cd036_thre_fact99.mat','c','s','p','day','selectNeuron','params');
params = {'thresholded','ar1','fit_params','optimize_smin', 'optimize_b','thresh_factor', 0.95};
[c,s,p] = deconvolveData(selectDff, params);
save('baseline\cd036_thre_fact95.mat','c','s','p','day','selectNeuron','params');

params = {'thresholded','ar1',0.9381,'smin',0.5, 'optimize_b','thresh_factor', 0.99};
[c,s,p] = deconvolveData(selectDff, params);
save('baseline\cd036_thre_fact99_smin05.mat','c','s','p','day','selectNeuron','params');
params = {'thresholded','ar1',0.9381,'smin',0.5, 'optimize_b','thresh_factor', 0.95};
[c,s,p] = deconvolveData(selectDff, params);
save('baseline\cd036_thre_fact95_smin05.mat','c','s','p','day','selectNeuron','params');


%% 2 - select the sessions - behavioral sessions
sessionSelected = [10 15 19 23 27 31 35]; 
day = 3:9;
[selectDff, selectNeuron] = selectSession(sessionSelected,day);
disp('Selection Done')

% confoo - 1. estimate baseline, 2. give baseline
params = {'constrained','ar1','fit_params', 'optimize_b'};[c,s,p] = deconvolveData(selectDff, params);%params = {'constrained','ar1', 0.9381, 'optimize_b'};[c,s,p] = deconvolveData(selectDff, params);
save('behavior\cd036_confoo.mat','c','s','p','day','selectNeuron','params');
params = {'constrained','ar1', 'fit_params', 'b',0};[c,s,p] = deconvolveData(selectDff, params);
save('behavior\cd036_confoo_b.mat','c','s','p','day','selectNeuron','params');
params = {'constrained','ar1','fit_params', 'optimize_b','s_min',-3};[c,s,p] = deconvolveData(selectDff, params);%params = {'constrained','ar1', 0.9381, 'optimize_b'};[c,s,p] = deconvolveData(selectDff, params);
save('behavior\cd036_confoo_smin-3.mat','c','s','p','day','selectNeuron','params');

% foo - 1. estimate baseline, 2. give baseline, 3. different lambda = 1.2,
% 4. diff lambda, give baseline, 5,6. different lambda = 3.6
params = {'foopsi','ar1', 'fit_params', 'optimize_b','lambda', 0};[c,s,p] = deconvolveData(selectDff, params);
save('behavior\cd036_foo_lambda0.mat','c','s','p','day','selectNeuron','params');
params = {'foopsi','ar1', 'fit_params', 'b',0,'lambda', 0};[c,s,p] = deconvolveData(selectDff, params);
save('behavior\cd036_foo_lambda0_b.mat','c','s','p','day','selectNeuron','params');

params = {'foopsi','ar1', 'fit_params', 'optimize_b','lambda', 0.5};[c,s,p] = deconvolveData(selectDff, params);
save('behavior\cd036_foo_lambda05.mat','c','s','p','day','selectNeuron','params');
params = {'foopsi','ar1', 'fit_params', 'b',0,'lambda', 0.5};[c,s,p] = deconvolveData(selectDff, params);
save('behavior\cd036_foo_lambda05_b.mat','c','s','p','day','selectNeuron','params');

params = {'foopsi','ar1', 'fit_params', 'optimize_b','lambda', 1.0};[c,s,p] = deconvolveData(selectDff, params);
save('behavior\cd036_foo_lambda10.mat','c','s','p','day','selectNeuron','params');
params = {'foopsi','ar1', 'fit_params', 'b',0,'lambda', 1.0};[c,s,p] = deconvolveData(selectDff, params);
save('behavior\cd036_foo_lambda10_b.mat','c','s','p','day','selectNeuron','params');

% thresholded - 1. optimize b, threshold = 0.99, 2. thre = 0.95, 3,4 . smin =
% -3
params = {'thresholded','ar1','fit_params','optimize_smin', 'optimize_b','thresh_factor', 0.99};
[c,s,p] = deconvolveData(selectDff, params);
save('behavior\cd036_thre_fact99.mat','c','s','p','day','selectNeuron','params');
params = {'thresholded','ar1','fit_params','optimize_smin', 'optimize_b','thresh_factor', 0.95};
[c,s,p] = deconvolveData(selectDff, params);
save('behavior\cd036_thre_fact95.mat','c','s','p','day','selectNeuron','params');

params = {'thresholded','ar1',0.9381,'smin',0.5, 'optimize_b','thresh_factor', 0.99};
[c,s,p] = deconvolveData(selectDff, params);
save('behavior\cd036_thre_fact99_smin05.mat','c','s','p','day','selectNeuron','params');
params = {'thresholded','ar1',0.9381,'smin',0.5, 'optimize_b','thresh_factor', 0.95};
[c,s,p] = deconvolveData(selectDff, params);
save('behavior\cd036_thre_fact95_smin05.mat','c','s','p','day','selectNeuron','params');

%% 3 - test the calmAn methods
sessionSelected = [13 18 22 26 30 34 38]; day = 3:9;
[selectDff, selectNeuron] = selectSession(sessionSelected,day);
[c,s,p] = deconvolveDataCalmAn(selectDff,'constrained_foopsi',[],true);
save('baseline\cd036_calman_confoo_optb.mat','c','s','p','day','selectNeuron');

[c,s,p] = deconvolveDataCalmAn(selectDff,'constrained_foopsi',[],false);
save('baseline\cd036_calman_confoo_optb_nosmin.mat','c','s','p','day','selectNeuron');


% behavioral sessions
sessionSelected = [10 15 19 23 27 31 35];  day = 3:9;
[selectDff, selectNeuron] = selectSession(sessionSelected,day);
[c,s,p] = deconvolveDataCalmAn(selectDff,'constrained_foopsi',[],true);
save('behavior\cd036_calman_confoo_optb.mat','c','s','p','day','selectNeuron');

[c,s,p] = deconvolveDataCalmAn(selectDff,'constrained_foopsi',[],false);
save('behavior\cd036_calman_confoo_optb_nosmin.mat','c','s','p','day','selectNeuron');

%% 4 - test the calmAn AR2 methods
sessionSelected = [13 18 22 26 30 34 38]; day = 3:9;
[selectDff, selectNeuron] = selectSession(sessionSelected,day);
[c,s,p] = deconvolveDataCalmAn(selectDff,'foopsi',0.90,true);
save('baseline\cd036_calman_ar2_foo90_optb.mat','c','s','p','day','selectNeuron');
[c,s,p] = deconvolveDataCalmAn(selectDff,'foopsi',0.90,false);
save('baseline\cd036_calman_ar2_foo90_optb_nosmin.mat','c','s','p','day','selectNeuron');

[c,s,p] = deconvolveDataCalmAn(selectDff,'foopsi',0.95,true);
save('baseline\cd036_calman_ar2_foo95_optb.mat','c','s','p','day','selectNeuron');
[c,s,p] = deconvolveDataCalmAn(selectDff,'foopsi',0.95,false);
save('baseline\cd036_calman_ar2_foo95_optb_nosmin.mat','c','s','p','day','selectNeuron');


% behavioral sessions
sessionSelected = [10 15 19 23 27 31 35];  day = 3:9;
[selectDff, selectNeuron] = selectSession(sessionSelected,day);
[c,s,p] = deconvolveDataCalmAn(selectDff,'foopsi',0.90,true);
save('behavior\cd036_calman_ar2_foo90_optb.mat','c','s','p','day','selectNeuron');
[c,s,p] = deconvolveDataCalmAn(selectDff,'foopsi',0.90,false);
save('behavior\cd036_calman_ar2_foo90_optb_nosmin.mat','c','s','p','day','selectNeuron');

[c,s,p] = deconvolveDataCalmAn(selectDff,'foopsi',0.95,true);
save('behavior\cd036_calm an_ar2_foo95_optb.mat','c','s','p','day','selectNeuron');
[c,s,p] = deconvolveDataCalmAn(selectDff,'foopsi',0.95,false);
save('behavior\cd036_calman_ar2_foo95_optb_nosmin.mat','c','s','p','day','selectNeuron');

%% stop parallel pooling
delete(parobj)

%% all functions
function [selectDff,selectNeuron] = selectSession(sessionSelected,day)
    global tempTC ishere nFrames;
    selectDff = cell(1,length(sessionSelected));
    selectNeuron = sum(ishere(:,day),2)==length(day);

    for i = 1:length(sessionSelected)
        sessionTC = tempTC(:,(nFrames(sessionSelected(i),1)+1):nFrames(sessionSelected(i)+1,1));
        sessionDff = (sessionTC - repmat(median(sessionTC,2),1,size(sessionTC,2))) ./ repmat(median(sessionTC,2),1,size(sessionTC,2));
        sessionDff = sessionDff - smoothdata(sessionDff,2,'movmedian',3000);
        selectDff{i} = sessionDff(selectNeuron,2:end);
    end
end

function topDff = selectTopNeurons(selectDff,nSelect)
    topActivity = prctile(selectDff{1},98,2);
    [~,sortIndex] = sort(topActivity,1,'ascend');
    neuronIndex = sortIndex(linspace(1,length(sortIndex),nSelect));
    for i = 1:length(sessionSelected)
        topDff{i} = selectDff{i}(neuronIndex,:);
    end
end

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


pax = spk_autocalibration('par');
pax.dt = 1/15.63;
% (set limits for A and tau)
pax.amin = amin;
pax.amax = amax;
pax.taumin = taumin;
pax.taumax = taumax;
% (set saturation parameter)
pax.saturation = 0.1;
% (give real values to the algorithm for display purposes - they obviously won't be used in the estimation!)
pax.realspikes = spikes;
pax.reala = a;
pax.realtau = tau;
% (when running MLspike from spk_autocalibratio, do not display graph summary)
pax.mlspikepar.dographsummary = false;

% perform auto-calibration
[tauest aest sigmaest] = spk_autocalibration(calcium,pax)


tic;
for i=1:length(selectDff)
    c{i} = zeros(nNeuron, size(selectDff{i},2));
    s{i} = zeros(nNeuron, size(selectDff{i},2));
    for j = 1:nNeuron
        tempDff = selectDff{i}(j,:);
       
        [c_deconv,s_deconv,ops] = deconvolveCa(tempDff,model_ar,'method',method,'optimize_pars',true,...
            'optimize_b',true,'maxIter',20,'window',150,'lambda',lam,'smin',spkmin);

        c{i}(j,:) = c_deconv;
        s{i}(j,:) = s_deconv;
        p{i,j} = ops;
    end  
end
t = toc;
disp(['Time Elapsed = ' num2str(t,'%.2f') ' secs.']);
end
