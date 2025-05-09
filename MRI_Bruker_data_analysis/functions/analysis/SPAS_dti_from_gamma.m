function out = SPAS_dti_from_gamma(nii_fn,fit_fn,wfm_name)
% DT from gamma fits saved in fit_fn      
% need data with LTE, e.g. geoSPAS
% outputs:
% out.md_nii_fn
% out.fa_nii_fn
% out.vd_nii_fn
% out.ap_dADC_nii_fn
% out.ap_rdADC_nii_fn


if nargin < 3
    wfm_name = 'geoSPAS'; %use waveform with wfm_name for DT fitting
end

% read fits
xps = mdm_xps_load(mdm_fn_nii2xps(nii_fn));

load(fit_fn)

% find waveform
wfm_ind = find(contains(res.wfm_names,wfm_name));
% find decay indices
decay_ind = find([res.decay.n_wfm] == wfm_ind)';

% find rotation vectors
u_ind = [res.decay(decay_ind).b_ind];
u_ind = u_ind(1,:)';
u = xps.u(u_ind,:);

% find antipodal indices
ap_ind = xps.ap_ind(u_ind);
uap_ind = unique(ap_ind);

% find ADCs for all voxels
ADC = res.gamma_fits(:,decay_ind,2);

% change units
ADC = ADC*1e9;

% average antipodal ADCs
for n_ap = 1:length(uap_ind);
    tmp_ind = find(ap_ind == uap_ind(n_ap));
    tmp = ADC(:,tmp_ind);
    ap_mADC(:,n_ap) = mean(tmp,2);
    ap_dADC(:,n_ap) = tmp(:,end)-tmp(:,1);
    u_ap(n_ap,:) = u(tmp_ind(1),:); % collapse antipodal directions
end


for n_vox = 1:size(ADC,1)

    Din = ap_mADC(n_vox,:);

    DT = tensor_from_projections(u_ap,Din);

    % check if eigenvalues are negative
    if any(DT.l<0)
        MD(n_vox,1) = 0;
        FA(n_vox,1) = 0;
    else
        MD(n_vox,1) = DT.m;
        FA(n_vox,1) = DT.FA;
    end

    V(n_vox,:) = DT.V(:,1);
end


md_map = zeros(res.img_size);
fa_map = md_map;

md_map(res.vox_mask_ind) = single(MD);
fa_map(res.vox_mask_ind) = single(FA);

out_nii_fn = append_nii_fn(nii_fn, 'md');
mdm_nii_write(md_map,out_nii_fn);
out.md_nii_fn = out_nii_fn;

out_nii_fn = append_nii_fn(nii_fn, 'fa');
mdm_nii_write(fa_map,out_nii_fn);
out.fa_nii_fn = out_nii_fn;


z = zeros(res.img_size);
for n = 1:3
    vd_map1 = z;
    vd_map1(res.vox_mask_ind) = single(V(:,n));
    vd_map(:,:,:,n) = vd_map1;
end

out_nii_fn = append_nii_fn(nii_fn, 'vd');
mdm_nii_write(vd_map,out_nii_fn);
out.vd_nii_fn = out_nii_fn;


for n = 1:size(ap_dADC,2)
    ap_dADC_map1 = z;
    ap_dADC_map1(res.vox_mask_ind) = single(ap_dADC(:,n));
    ap_dADC_map(:,:,:,n) = ap_dADC_map1;

    ap_mADC_map1 = z;
    ap_mADC_map1(res.vox_mask_ind) = single(ap_mADC(:,n));
    ap_mADC_map(:,:,:,n) = ap_mADC_map1;

    ap_rdADC_map1 = z;
    tmp = ap_dADC_map1./ap_mADC_map1;
    tmp(isnan(tmp)) = 0;
    ap_rdADC_map(:,:,:,n) = tmp; %ap_dADC_map1./ap_mADC_map1;

end

out_nii_fn = append_nii_fn(nii_fn, 'ap_mADC');
mdm_nii_write(ap_mADC_map,out_nii_fn);
out.ap_dADC_nii_fn = out_nii_fn;

out_nii_fn = append_nii_fn(nii_fn, 'ap_dADC');
mdm_nii_write(ap_dADC_map,out_nii_fn);
out.ap_dADC_nii_fn = out_nii_fn;

out_nii_fn = append_nii_fn(nii_fn, 'ap_rdADC');
mdm_nii_write(ap_rdADC_map,out_nii_fn);
out.ap_rdADC_nii_fn = out_nii_fn;





