function out = SPAS_sig_dif_to_nifti(nii_fn, mask_nii_fn, opt)
% extract signal difference from SPAS1-SPAS3 and/or geoSPAS-STE at max b
% extract average signals at min/max b

if nargin < 3
    opt.normalize_case = 2; % 0-no normalization, 1-divide by smin, 2-divide by smax

    opt.save.TDD_smin = 1;
    opt.save.TDD_smax = 1;
    opt.save.TDD_ds = 1;
    opt.save.TDD_dlogs = 1;

    opt.save.muA_smin = 1;
    opt.save.muA_smax = 1;
    opt.save.muA_ds = 1;
    opt.save.muA_dlogs = 1;
end
% read image
xps = mdm_xps_load(mdm_fn_nii2xps(nii_fn));
[sig, h_sig] = mdm_nii_read(nii_fn);
sig = double(sig);
% read mask
mask = single(mdm_nii_read(mask_nii_fn));
mask = double(mask);

% find wfm indices for 1-SPAS1, 2-SPAS2, 3-SPAS3, 4-STE, 5-tLTE, 6-geoSPAS in this order
SPAS_wfms = SPAS_order(xps);

bmax_ind = xps.a_ind == max(xps.a_ind); % a_ind is linked to b value
bmin_ind = xps.a_ind == min(xps.a_ind);

if SPAS_wfms.contains_SPAS
    ind1_bmin = xps.s_ind == SPAS_wfms.ordered_wfm_ind(1) & bmin_ind; % SPAS1
    ind2_bmin = xps.s_ind == SPAS_wfms.ordered_wfm_ind(3) & bmin_ind; % SPAS3

    ind1_bmax = xps.s_ind == SPAS_wfms.ordered_wfm_ind(1) & bmax_ind; % SPAS1
    ind2_bmax = xps.s_ind == SPAS_wfms.ordered_wfm_ind(3) & bmax_ind; % SPAS3

    TDD_smin = (sig(:,:,:,ind1_bmin) + sig(:,:,:,ind2_bmin) ) / 2;
    TDD_smax = (sig(:,:,:,ind1_bmax) + sig(:,:,:,ind2_bmax) ) / 2;

    TDD_ds = sig(:,:,:,ind1_bmax) - sig(:,:,:,ind2_bmax);
    TDD_dlogs = log(sig(:,:,:,ind1_bmax)) - log(sig(:,:,:,ind2_bmax));

    % normalize
    if opt.normalize_case > 0
        if opt.normalize_case == 1
            TDD_ds = TDD_ds ./ TDD_smin;
        elseif opt.normalize_case == 2
            TDD_ds = TDD_ds ./ TDD_smax;
        end
    end

    % mask
    TDD_smin = mask .* TDD_smin;
    TDD_smax = mask .* TDD_smax;

    TDD_ds = mask .* TDD_ds;
    TDD_dlogs = mask .* TDD_dlogs;

    % NaN -> 0, Inf -> 0
    TDD_ds(isnan(TDD_ds)) = 0;
    TDD_ds(abs(TDD_ds) == Inf) = 0;
    TDD_dlogs(isnan(TDD_dlogs)) = 0;
    TDD_dlogs(abs(TDD_dlogs) == Inf) = 0;


    % save
    out_nii_fn = append_nii_fn(nii_fn, 'raw_TDD_smin');
    mdm_nii_write(single(TDD_smin), out_nii_fn, h_sig);
    out.raw_TDD_smin_nii_fn = out_nii_fn;

    out_nii_fn = append_nii_fn(nii_fn, 'raw_TDD_smax');
    mdm_nii_write(single(TDD_smax),out_nii_fn, h_sig);
    out.raw_TDD_smax_nii_fn = out_nii_fn;

    out_nii_fn = append_nii_fn(nii_fn, 'raw_TDD_ds');
    mdm_nii_write(single(TDD_ds),out_nii_fn, h_sig);
    out.raw_TDD_ds_nii_fn = out_nii_fn;

    out_nii_fn = append_nii_fn(nii_fn, 'raw_TDD_dlogs');
    mdm_nii_write(single(TDD_dlogs),out_nii_fn, h_sig);
    out.raw_TDD_dlogs_nii_fn = out_nii_fn;

end

if SPAS_wfms.contains_STE & SPAS_wfms.contains_geoSPAS

    ind1_bmin = xps.s_ind == SPAS_wfms.ordered_wfm_ind(6) & bmin_ind; % geoSPAS
    ind2_bmin = xps.s_ind == SPAS_wfms.ordered_wfm_ind(4) & bmin_ind; % STE

    ind1_bmax = xps.s_ind == SPAS_wfms.ordered_wfm_ind(6) & bmax_ind; % geoSPAS
    ind2_bmax = xps.s_ind == SPAS_wfms.ordered_wfm_ind(4) & bmax_ind; % STE

    muA_smin = (sig(:,:,:,ind1_bmin) + sig(:,:,:,ind2_bmin) ) / 2;
    muA_smax = (sig(:,:,:,ind1_bmax) + sig(:,:,:,ind2_bmax) ) / 2;

    muA_ds = sig(:,:,:,ind1_bmax) - sig(:,:,:,ind2_bmax);
    muA_dlogs = log(sig(:,:,:,ind1_bmax)) - log(sig(:,:,:,ind2_bmax));

    % normalize
    if opt.normalize_case > 0
        if opt.normalize_case == 1
            muA_ds = muA_ds ./ muA_smin;
        elseif opt.normalize_case == 2
            muA_ds = muA_ds ./ muA_smax;
        end
    end

    % NaN -> 0
    muA_ds(isnan(muA_ds)) = 0;
    muA_dlogs(isnan(muA_dlogs)) = 0;

    % Inf -> 0
    muA_ds(isinf(muA_ds)) = 0;
    muA_dlogs(isinf(muA_dlogs)) = 0;

    % mask
    muA_smin = mask .* muA_smin;
    muA_smax = mask .* muA_smax;

    muA_ds = mask .* muA_ds;
    muA_dlogs = mask .* muA_dlogs;


    % save
    out_nii_fn = append_nii_fn(nii_fn, 'raw_muA_smin');
    mdm_nii_write(single(muA_smin),out_nii_fn, h_sig);
    out.raw_muA_smin_nii_fn = out_nii_fn;

    out_nii_fn = append_nii_fn(nii_fn, 'raw_muA_smax');
    mdm_nii_write(single(muA_smax),out_nii_fn, h_sig);
    out.raw_muA_smax_nii_fn = out_nii_fn;

    out_nii_fn = append_nii_fn(nii_fn, 'raw_muA_ds');
    mdm_nii_write(single(muA_ds),out_nii_fn, h_sig);
    out.raw_muA_ds_nii_fn = out_nii_fn;

    out_nii_fn = append_nii_fn(nii_fn, 'raw_muA_dlogs');
    mdm_nii_write(single(muA_dlogs),out_nii_fn, h_sig);
    out.raw_muA_dlogs_nii_fn = out_nii_fn;


end


