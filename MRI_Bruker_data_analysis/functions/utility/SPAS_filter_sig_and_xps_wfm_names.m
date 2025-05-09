function [sig, xps] = SPAS_filter_sig_and_xps_wfm_names(sig, xps, select_wfm_names)
% return s and xps for selected wfm names 

wfm_ind = find(contains(xps.wfm_names, select_wfm_names));
ind = find(ismember(xps.s_ind, wfm_ind ));

% make new xps
f = fieldnames(xps);
exclude_fields = {'n','wfm_names','wfm_src'};
f = f(~ismember(f,exclude_fields));
xps_sub.n = length(ind);

for i = 1:numel(f)
    xps_sub.(f{i}) = xps.(f{i})(ind, :,:);
end
% write back excluded fields containing strings 
if isfield(xps,'wfm_names') xps_sub.wfm_names = xps.wfm_names(wfm_ind); end
if isfield(xps,'wfm_src') xps_sub.wfm_src = xps.wfm_src(wfm_ind); end

sig = sig(ind);
xps = xps_sub;
