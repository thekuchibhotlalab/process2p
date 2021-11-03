function func_spikeDeconv(TC, varargin)
% OUTDATED -- use deconvolveDff instead
if isempty(varargin)
    [c, s, options] = deconvolveCa(submat(k,:),'constrained','ar1', 0.9381,...
        'optimize_b','optimize_pars');% time x one cell
    c_deconvolve_by_trial = c;
    s_deconvolve_by_trial = s;
    params{i,j,k} = options;
    return;
end
p = func_createInputParser();
p.addParameter('file', [])
p.parse(varargin{:});

[c, s, options] = deconvolveCa(submat(k,:),'constrained','ar1', 0.9381,...
 'optimize_b','optimize_pars');% time x one cell
c_deconvolve_by_trial = c;
s_deconvolve_by_trial = s;
params{i,j,k} = options;



end