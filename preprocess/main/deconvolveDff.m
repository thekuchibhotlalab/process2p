function [s,c,p] = deconvolveDff(sessionDff,varargin)
% [s,c,p] = deconvolveDff(sessionDff,varargin)
% Input Arguments -- 
% sessionDff: cell array (multi-session) or matrix (single-session). 
%   Each session shaped in nNeuronXframes. 
%
% Varargin --
% 'pars': optinoal decay parameter of each neuron, shape of [nNeuron,1] for
% dacay or [nNeuron,2] for both rise and decay. If given, shared between
% all sessions.
% 'fr': frame rate. default = 15.63.
% 'decay_time': decay time (sec) is only used to set sparness parameter. 
%   default = 0.7
% 'method': algorithm. 'foopsi'(default), 'constrained', 'threshold'
% 'model_ar': order of the model. 'ar1', 'ar2'(default), 'exp2' or 'kernel'
% 'lam_pr': prob of spk in each frame, sparseness parameter. default = 0.9
% 'spk_SNR': minimum signal strength (in std unit) to detect a spike. default = 0
%   reasonable choice are from 0.1 - 1 (std)
%
% Output Arguments --
% s: deconvolved spks. c: denoised calcium trace. p: parameters.
% If input is cell/matrix, output is corresponding cell/matrix. 
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('pars', [])
p.parse(varargin{:});
pars = p.Results.pars;

tic;
if iscell(sessionDff)
    c = cell(1,length(sessionDff));
    s = cell(1,length(sessionDff));
    p = cell(1,length(sessionDff));
    for i=1:length(sessionDff)
        nNeuron = size(sessionDff{i},1);
        c{i} = zeros(size(sessionDff{i}));s{i} = zeros(size(sessionDff{i}));p{i} = cell(1,nNeuron);
        for j = 1:nNeuron
            neuronDff = sessionDff{i}(j,:);
            if isempty(pars)
                [c_deconv,s_deconv,ops] = deconvolveNeuron(neuronDff,[],varargin{:});
            else 
                [c_deconv,s_deconv,ops] = deconvolveNeuron(neuronDff,pars(j,:),varargin{:});
            end
            c{i}(j,:) = c_deconv;
            s{i}(j,:) = s_deconv;
            p{i}{j} = ops;
        end  
        t = toc;
        disp(['Session ' int2str(i) ' Done. Total Time Elapsed = ' num2str(t,'%.2f') ' secs.']);
    end
else
    nNeuron = size(sessionDff,1);
    c = zeros(size(sessionDff));s = zeros(size(sessionDff));p = cell(1,nNeuron);
    for j = 1:nNeuron
        neuronDff = sessionDff(j,:);
        if isempty(pars)
            [c_deconv,s_deconv,ops] = deconvolveNeuron(neuronDff,[],varargin{:});
        else 
            [c_deconv,s_deconv,ops] = deconvolveNeuron(neuronDff,pars(j,:),varargin{:});
        end
        c(j,:) = c_deconv;
        s(j,:) = s_deconv;
        p{j} = ops;
    end  
    t = toc;
    disp(['Session Done. Time Elapsed = ' num2str(t,'%.2f') ' secs.']);
end

end


function [c_deconv,s_deconv,ops] = deconvolveNeuron(neuronDff,pars,varargin)
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('fr', 15.63)
p.addParameter('decay_time', 0.7)
p.addParameter('model_ar', 'ar1')
p.addParameter('method', 'foopsi')
p.addParameter('lam_pr', 0.9)
p.addParameter('spk_SNR', 0)
p.parse(varargin{:});

fr = p.Results.fr; decay_time = p.Results.decay_time;
spk_SNR = p.Results.spk_SNR; lam_pr = p.Results.lam_pr;
model_ar = p.Results.model_ar; method = p.Results.method;
spkmin = spk_SNR*GetSn(neuronDff); 

if ~strcmp(method,'constrained')
    lam = choose_lambda(exp(-1/(fr*decay_time)),GetSn(neuronDff),lam_pr);
    if isempty(pars)
        [c_deconv,s_deconv,ops] = deconvolveCa(neuronDff,model_ar,'method',method,'optimize_pars',true,...
            'optimize_b',true,'maxIter',20,'lambda',lam,'smin',spkmin);
    else
        [c_deconv,s_deconv,ops] = deconvolveCa(neuronDff,model_ar,'method',method,'pars',pars,...
            'optimize_b',true,'maxIter',20,'lambda',lam,'smin',spkmin); 
    end
else 
    if isempty(pars)
        [c_deconv,s_deconv,ops] = deconvolveCa(neuronDff,model_ar,'method',method,'optimize_pars',true,...
            'optimize_b',true,'maxIter',20,'smin',spkmin);
    else
        [c_deconv,s_deconv,ops] = deconvolveCa(neuronDff,model_ar,'method',method,'pars',pars,...
            'optimize_b',true,'maxIter',20,'smin',spkmin);
    end
end

end