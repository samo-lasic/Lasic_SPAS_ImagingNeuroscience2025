function out_nii_fn = SPAS_average_repetitions(nii_fn)
% average repetitions from Bruker scans

xps = mdm_xps_load(mdm_fn_nii2xps(nii_fn));
[sig, h] = mdm_nii_read(nii_fn);
sig = double(sig);

rep_a_ind = unique(xps.rep_a_ind);
s_ind = unique(xps.s_ind);
sig_a = sig;

if length(rep_a_ind) > 1
    ind = [];
    for n_wfm = 1:length(s_ind)
        for n_rep = 1:length(rep_a_ind)
            ind_tmp = find(xps.rep_a_ind == rep_a_ind(n_rep) & xps.s_ind == s_ind(n_wfm));
%             xps.rep_dim_ind(ind_tmp)
            ind = [ind; ind_tmp(1)];
            sig_a(:,:,:,ind_tmp(1)) = mean(sig_a(:,:,:,ind_tmp),4);
        end
    end
end

% reduce to the indices containing mean values
sig_a = sig_a(:,:,:,ind);

% make new xps
f = fieldnames(xps);
exclude_fields = {'n','wfm_names','wfm_src'};
f = f(~ismember(f,exclude_fields));
xps_a.n = length(ind);

for i = 1:numel(f)
    xps_a.(f{i}) = xps.(f{i})(ind, :,:);
end
xps_a.wfm_names = xps.wfm_names;
if isfield(xps,'wfm_src')
    xps_a.wfm_src = xps.wfm_src;
end


out_nii_fn = append_nii_fn(nii_fn, 'a');
mdm_nii_write(single(sig_a),out_nii_fn, h);
mdm_xps_save(xps_a, mdm_fn_nii2xps(out_nii_fn));


