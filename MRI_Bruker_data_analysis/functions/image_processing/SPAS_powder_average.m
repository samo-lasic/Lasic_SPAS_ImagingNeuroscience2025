function out_nii_fn = SPAS_powder_average(nii_fn)
% powder average

% read image
xps = mdm_xps_load(mdm_fn_nii2xps(nii_fn));
[sig, h] = mdm_nii_read(nii_fn);


% arithmetic average over directions
uwfm = unique(xps.s_ind);
udir = unique(xps.rot_ind);
ub = unique(xps.a_ind);
cnt = 0;
for n_wfm = 1:length(uwfm)
    for n_b = 1:length(ub)
        ind = find(uwfm(n_wfm) == xps.s_ind & ub(n_b) == xps.a_ind);
        
        if length(ind) > 1
            cnt = cnt+1;
            sig_a(:,:,:,cnt) = mean(sig(:,:,:,ind),4);
            ind_a(cnt) = ind(1);
        end
        
    end
end


f = fieldnames(xps);
exclude_fields = {'n','wfm_names','wfm_src','v','DwShapeDir','DwR','DwR3x3','u','rot_ind'};
f = f(~ismember(f,exclude_fields));
xps_a.n = length(ind_a);

for i = 1:numel(f)
    xps_a.(f{i}) = xps.(f{i})(ind_a, :,:);
end

xps_a.wfm_names = xps.wfm_names;

out_nii_fn = append_nii_fn(nii_fn, 'pa');
mdm_nii_write(single(sig_a), out_nii_fn, h);
mdm_xps_save(xps_a, mdm_fn_nii2xps(out_nii_fn));




