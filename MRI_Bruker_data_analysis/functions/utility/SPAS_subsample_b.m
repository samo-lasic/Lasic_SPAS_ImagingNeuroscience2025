function out_nii_fn = SPAS_subsample_b(nii_fn, b_ind)
% subsample depending on b-value, xps can contain string fields 

xps = mdm_xps_load(mdm_fn_nii2xps(nii_fn));
sig = double(mdm_nii_read(nii_fn));

ub = uniquetol(xps.b,1e-5);
ub = ub(b_ind)

ind = [];
for n = 1:length(ub)
    ind = [ind; find(abs(xps.b - ub(n)) < 1e5)];
end

% make new xps
f = fieldnames(xps);
exclude_fields = {'n','wfm_names','wfm_src'};
f = f(~ismember(f,exclude_fields));
xps_sub.n = length(ind);

for i = 1:numel(f)
    xps_sub.(f{i}) = xps.(f{i})(ind, :,:);
end
% give back excluded fields containing strings 
if isfield(xps,'wfm_names') xps_sub.wfm_names = xps.wfm_names; end
if isfield(xps,'wfm_src') xps_sub.wfm_src = xps.wfm_src; end


sig = sig(:,:,:,ind);

%out_nii_fn = append_nii_fn(nii_fn, 'sub');
out_nii_fn = nii_fn;
mdm_nii_write(single(sig),out_nii_fn);
mdm_xps_save(xps_sub, mdm_fn_nii2xps(out_nii_fn));
