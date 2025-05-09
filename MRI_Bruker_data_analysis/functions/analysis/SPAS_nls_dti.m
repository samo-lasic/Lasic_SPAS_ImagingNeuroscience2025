
function SPAS_nls_dti(nii_fn, mask_nii_fn, opt)
% dti fit for all waveforms in opt.wfm_name
if ~iscell(opt.wfm_name)
    wfm_name = {opt.wfm_name};
else
    wfm_name = opt.wfm_name;
end

for n = 1:numel(wfm_name)
    opt.wfm_name = wfm_name{n};
    out = SPAS_nls_dti1(nii_fn, mask_nii_fn, opt);

    if isempty(out)
        display(sprintf('%s','missing data'))
    else
        f = fields(out);
        for n = 1:numel(f)
            display(sprintf('%s',out.(f{n})))
        end
    end
end

end

function out = SPAS_nls_dti1(nii_fn, mask_nii_fn, opt)
% outputs
% out.norm_nii_fn
% out.norm_fit_nii_fn
% out.s0_nii_fn
% out.adc_nii_fn
% out.mu2_nii_fn
% out.norm_par_mat_fn

if ~isfield(opt,'do_thresh_vol')
    opt.do_thresh_vol = 0;
end

% read image
xps = mdm_xps_load(mdm_fn_nii2xps(nii_fn));
[sig, h_sig] = mdm_nii_read(nii_fn);

% read mask
mask = mdm_nii_read(mask_nii_fn);

% filter volumes
if opt.do_thresh_vol & isfield(xps,'good_vol')
    ind = xps.good_vol;
else
    ind = ones(xps.n,1);
end

% collect indices for selected waveform
n_wfm = find(contains(xps.wfm_names,opt.wfm_name));

all_ind = xps.s_ind == n_wfm;
ind = find(all_ind & ind);

n_all_ind = length(find(all_ind));
n_good_ind = length(ind);

% eliminate single b (single rotations)
rot_ind = xps.rot_ind(ind);
for n = 1:length(rot_ind)
    cnt(n) = length(find(rot_ind == rot_ind(n)));
end
ind = ind(cnt > 1);

% select b-values
if opt.b_lim == 0 % select first 2
    opt.b_lim = uniquetol(xps.b(ind),1e-6);
    if opt.b_lim < 3
        opt.b_lim = opt.b_lim(2)+1e6;
    else
        opt.b_lim = opt.b_lim(3);
    end
end

ind = ind(xps.b(ind) < opt.b_lim);
display(sprintf('uses %d b-selected volumes out of %d good ones and %d total', length(ind), n_good_ind, n_all_ind))

% subsample xps
f = fieldnames(xps);
exclude_fields = {'n','wfm_names','wfm_src'};
f = f(~ismember(f,exclude_fields));
xps_sub.n = length(ind);
for i = 1:numel(f)
    xps_sub.(f{i}) = xps.(f{i})(ind, :,:);
end
xps_sub.wfm_names = opt.wfm_name; %xps.wfm_names;
xps = xps_sub;

% check that minimum fitting requirements are satified
% i.e. min 3 non collinear directions

% reshape signal into 2D
sig_sz = size(sig);
n_voxels = prod(sig_sz(1:3));
sig = reshape(sig,n_voxels,sig_sz(4));
sig = sig(:,ind);

vox_mask_ind = find(mask);
n_voxels2fit = length(vox_mask_ind);

nls_dti_fit_parameters = zeros(n_voxels2fit,7);
nls_dti_maps = zeros(n_voxels2fit,8);

% Verify the xps
dti_nls_check_xps(xps);
opt = dti_nls_opt(opt);

% parfor
tic
parfor n_vox_mask_ind = 1:n_voxels2fit
    %for n_vox_mask_ind = 1:n_voxels2fit
    vox_ind = vox_mask_ind(n_vox_mask_ind);

    signal = sig(vox_ind,:)';

    % testing inversion
    % u = [1 0 0; ...
    %     0 1 0;
    %     0 0 1];

    %     rotL = squeeze(rotate_diagonal_tensors_around_XYZ([0 0 1],[0 1.1*pi/4 1*pi/4;]));
    %     b = [1 6]*1e8;
    %     c = 0;
    %     for nb = 1:2
    %     for i = 1:size(u,1)
    %         c = c+1;
    %         D(c) = u(i,:)*rotL*u(i,:)' * 1e-9;
    %         xps.bt(c,:) = b(nb) * tm_3x3_to_1x6(u(i,:)'*u(i,:));
    %         signal(c) = exp(-b(nb)*D(i));
    %     end
    %     end

    %opt.dti_nls.do_plot = 1;

    m = dti_nls_1d_data2fit(signal, xps, opt);

    %     signal_fit = dti_nls_1d_fit2data(m, xps);
    %     figure(1),clf
    %     semilogy(xps.b,signal,'.',xps.b,signal_fit,'x');

    DT = tm_1x6_to_3x3(m(2:end));
    T = tensor_info(DT);
    v = abs(T.V(1,:)); % main eigen vector

    nls_dti_fit_parameters(n_vox_mask_ind,:) = m;
    nls_dti_maps(n_vox_mask_ind,:) = [m(1) T.m T.FA T.Delta T.Eta v];

end
toc

% save maps

if opt.save_s0_map
    s0_map = zeros(n_voxels,1);
    s0_map(vox_mask_ind,1) = single(nls_dti_maps(:,1));

    s0_map = reshape(s0_map,sig_sz(1:3)); %reshape back to image dimensions

    out_nii_fn = append_nii_fn(nii_fn, ['dti_s0_' opt.wfm_name]);
    mdm_nii_write(s0_map, out_nii_fn, h_sig);

    out.s0_map_nii_fn = out_nii_fn;
end
if opt.save_md_map
    md_map = zeros(n_voxels,1);
    md_map(vox_mask_ind,1) = single(nls_dti_maps(:,2));

    md_map = reshape(md_map,sig_sz(1:3)); %reshape back to image dimensions

    out_nii_fn = append_nii_fn(nii_fn, ['dti_md_' opt.wfm_name]);
    mdm_nii_write(md_map,out_nii_fn, h_sig);

    out.md_map_nii_fn = out_nii_fn;
end
if opt.save_fa_map
    fa_map = zeros(n_voxels,1);
    fa_map(vox_mask_ind,1) = single(nls_dti_maps(:,3));

    fa_map = reshape(fa_map,sig_sz(1:3)); %reshape back to image dimensions

    out_nii_fn = append_nii_fn(nii_fn, ['dti_fa_' opt.wfm_name]);
    mdm_nii_write(fa_map,out_nii_fn, h_sig);

    out.fa_map_nii_fn = out_nii_fn;
end
if opt.save_v_map
    v_map = zeros(n_voxels,3);
    v_map(vox_mask_ind,:) = single(nls_dti_maps(:,6:end));

    v_map = reshape(v_map,[sig_sz(1:3) 3]); %reshape back to image dimensions

    out_nii_fn = append_nii_fn(nii_fn, ['dti_v_' opt.wfm_name]);
    mdm_nii_write(v_map,out_nii_fn, h_sig);

    out.v_map_nii_fn = out_nii_fn;
end

% save normalized signals
% out_nii_fn = append_nii_fn(nii_fn, 'dti_norm');
% mdm_nii_write(single(sig_norm),out_nii_fn, h_sig);
% mdm_xps_save(xps, mdm_fn_nii2xps(out_nii_fn));
% out.norm_nii_fn = out_nii_fn;

% save fit parameters (mat file)
if opt.save_dti_fit_parameters

    % add vox_mask_ind info
    res.dti_fits = nls_dti_fit_parameters;
    res.vox_mask_ind = vox_mask_ind;
    res.img_size = sig_sz(1:3);

    mat_fn = append_mat_fn(nii_fn, ['nls_dti_p_' opt.wfm_name]);
    save(mat_fn,'res')
    out.mat_fn = mat_fn;

    mdm_xps_save(xps, mdm_fn_nii2xps(mat_fn));

end
end

