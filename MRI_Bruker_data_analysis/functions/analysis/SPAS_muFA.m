function out = SPAS_muFA(nii_fn, mask_nii_fn, opt)
% fit muFA from STE and geoSPAS data

if nargin < 3

    opt.data_is_normalized = 0; % if 1 then guess s0 = 1 else estimate from data
    opt.b_low_lim = 600 *1e6; % for initial slope estimation, if 0 then use first two decay points
    opt.weighted = 0;
    opt.thresh = 1; % 1-all pass, .8 - only 20% best fits pass, .85 with NRMSE gives similar threshold as .7 with RSS
    opt.thresh_with_NRMSE = 1;
    opt.muFA_lim = []; %[] or e.g. [0 1]; % clamp maps if array is not empty

    opt.save.muFA = 1;
    opt.save.NRMSE = 0;
    opt.save.RSS = 0;
    opt.save.mask = 0;

end

if ~isfield(opt,'LTE_wfm_name')
    opt.LTE_wfm_name = {'geoSPAS'};
end

if ~iscell(opt.LTE_wfm_name) % ensure cell
    opt.LTE_wfm_name = {opt.LTE_wfm_name};
end


% read image
xps = mdm_xps_load(mdm_fn_nii2xps(nii_fn));
[sig, h_sig] = mdm_nii_read(nii_fn);

% read mask
mask = mdm_nii_read(mask_nii_fn);

% find wfm indices for 1-SPAS1, 2-SPAS2, 3-SPAS3, 4-STE, 5-tLTE, 6-geoSPAS in this order
SPAS_wfms = SPAS_order(xps);

% muFA gamma fit
if SPAS_wfms.contains_STE & numel(find(contains(SPAS_wfms.names,opt.LTE_wfm_name))) == 1 %SPAS_wfms.contains_geoSPAS

    % reshape signal into 2D
    sig_sz = size(sig);
    im_sz = sig_sz(1:3);
    n_voxels = prod(im_sz);
    sig = reshape(sig,n_voxels,sig_sz(4));
    vox_mask_ind = find(mask);

    % fit options
    % opt.fit.data_is_normalized determines if this guess is going to be used or we estimate from data
    if opt.weighted
        fit_opt = fit_muFA_weighted_opt(opt);
    else
        fit_opt = fit_muFA_opt(opt);
    end


    for nLTE = 1:numel(opt.LTE_wfm_name)
        out(nLTE) = make_muFA_maps(nii_fn, h_sig, sig, n_voxels, vox_mask_ind, im_sz, xps, ...
            'STE', opt.LTE_wfm_name{nLTE}, opt, fit_opt);
    end


else
    display('missing STE or geoSPAS data')
    out = [];
end
end

% fitting
function out = make_muFA_maps(nii_fn, h_sig, sig, n_voxels, vox_mask_ind, im_sz, xps, STE_wfm_name, LTE_wfm_name, opt, fit_opt)
n_voxels2map = length(vox_mask_ind);

% select STE
ind = xps.s_ind == find(contains(xps.wfm_names,STE_wfm_name));
voxSig2fit(:,:,2) = double(sig(:,ind)); % STE
b1 = xps.b(ind);

% choose LTE to use for fitting
ind = xps.s_ind == find(contains(xps.wfm_names,LTE_wfm_name));
voxSig2fit(:,:,1) = double(sig(:,ind));
b2 = xps.b(ind);

if any(abs(b1-b2) > 1e5) % check if b are different
    display('different b values!')
end

b = [b1; b2];

% parfor lists
muFA_fit_par_list = zeros(n_voxels2map,4);

NRMSE_list = zeros(n_voxels2map,1); % normalised root mean square error
RSS_list = zeros(n_voxels2map,1); % residual sum of squares
mask_list = ones(n_voxels2map,1);

tic
% for n_vox_mask_ind = 1:n_voxels2map % debugging
parfor n_vox_mask_ind = 1:n_voxels2map

    vox_ind = vox_mask_ind(n_vox_mask_ind);
    sig2fit = squeeze(voxSig2fit(vox_ind,:,:));
    sig2fit = sig2fit(:);

    %uses fit_muFA, alternatively use fit_muFA_dmu2

    if opt.weighted
        Pout = fit_muFA_weighted_par(b,sig2fit(:),fit_opt);
    else
        Pout = fit_muFA_par(b,sig2fit(:),fit_opt);
    end

    sig_fit = fit_muFA(Pout,b);
    NRMSE_list(n_vox_mask_ind) = sqrt(mean((sig_fit-sig2fit).^2))/(max(sig2fit)-min(sig2fit)); % normalised
    RSS_list(n_vox_mask_ind) = sum((sig_fit-sig2fit).^2); % residual sum of squares

    muFA_fit_par_list(n_vox_mask_ind,:) = Pout;

end
toc

% find thresholds
if opt.thresh_with_NRMSE % NRMSE
    error_list = NRMSE_list;
else % RSS
    error_list = RSS_list;
end

norm_sort_error_list = sort(error_list)/max(error_list);
thresh_error = max(error_list) * norm_sort_error_list(round(n_voxels2map * opt.thresh));

dmu2tilda = muFA_fit_par_list(:,4)./muFA_fit_par_list(:,2).^2;
muFA_list = sqrt(3/2)*(1+2/5*1./dmu2tilda).^(-1/2);

% mask bad fits
ind = error_list > thresh_error;
mask_list(ind) = 0;
RSS_list(ind) = 0;
NRMSE_list(ind) = 0;


muFA_fit_par_list(ind,:) = 0;
muFA_list(ind) = 0;

if strcmp(LTE_wfm_name,'geoSPAS') % if it is geoSPAS, not appending name
    LTE_wfm_name_append = '';
else
    LTE_wfm_name_append = [LTE_wfm_name '_'];
end

if opt.save.mask
    out.mask_nii_fn = append_nii_fn(nii_fn, [LTE_wfm_name_append 'muFA_mask']);
    img_map = make_map(mask_list, vox_mask_ind, n_voxels, im_sz, []);
    mdm_nii_write(int8(img_map), out.mask_nii_fn, h_sig);
end

if opt.save.NRMSE
    out.NRMSE_nii_fn = append_nii_fn(nii_fn, [LTE_wfm_name_append 'muFA_NRMSE']);
    img_map = make_map(NRMSE_list, vox_mask_ind, n_voxels, im_sz, []);
    mdm_nii_write(single(img_map), out.NRMSE_nii_fn, h_sig);
end

if opt.save.RSS
    out.RSS_nii_fn = append_nii_fn(nii_fn, [LTE_wfm_name_append 'muFA_RSS']);
    img_map = make_map(RSS_list, vox_mask_ind, n_voxels, im_sz, []);
    mdm_nii_write(single(img_map), out.RSS_nii_fn, h_sig);
end

if opt.save.s0
    out.s0_nii_fn = append_nii_fn(nii_fn, [LTE_wfm_name_append 'muFA_s0']);
    img_map = make_map(muFA_fit_par_list(:,1), vox_mask_ind, n_voxels, im_sz, opt.s0_lim);
    mdm_nii_write(single(img_map), out.s0_nii_fn, h_sig);
end

if opt.save.md
    out.md_nii_fn = append_nii_fn(nii_fn, [LTE_wfm_name_append 'muFA_md']);
    img_map = make_map(muFA_fit_par_list(:,2), vox_mask_ind, n_voxels, im_sz, opt.md_lim);
    mdm_nii_write(single(img_map), out.md_nii_fn, h_sig);
end

if opt.save.muFA
    out.uFA_nii_fn = append_nii_fn(nii_fn, [LTE_wfm_name_append 'muFA']);
    img_map = make_map(muFA_list, vox_mask_ind, n_voxels, im_sz, opt.muFA_lim);
    mdm_nii_write(single(img_map), out.uFA_nii_fn, h_sig);
end

if opt.save.muFA2
    out.uFA2_nii_fn = append_nii_fn(nii_fn, [LTE_wfm_name_append 'muFA2']);
    img_map = make_map(muFA_list.^2, vox_mask_ind, n_voxels, im_sz, opt.muFA2_lim);
    mdm_nii_write(single(img_map), out.uFA2_nii_fn, h_sig);
end

end

function img_map = make_map(img_list, vox_mask_ind, n_voxels, im_sz, lim)
% unpack img_list
img_map = zeros(n_voxels,1);
img_map(vox_mask_ind) = img_list;

% limit img_map
if ~isempty(lim)
    img_map(img_map < lim(1)) = lim(1);
    img_map(img_map > lim(2)) = lim(2);
end

%reshape back to image dimensions
img_map = reshape(img_map,im_sz);

end


% function used for debugging
function check_some_fits(n, thresh_error, error_list, b, b1, b2, voxSig2fit, Pout_list, vox_mask_ind)
% check some good fits
ind = find(error_list < thresh_error);
[~,sort_ind] = sort(error_list(ind));

if n == 0
    ind = ind(sort_ind(end));
else
    ind = ind(sort_ind(n));
end

vox_ind = vox_mask_ind(ind);
Pout = Pout_list(ind,:);
sig_fit = fit_muFA(Pout,b);
sig2fit = squeeze(voxSig2fit(vox_ind,:,:));
sig2fit = sig2fit(:);
figure(1),clf
semilogy(b1,sig2fit(1:length(b1)),'xr',b1,sig_fit(1:length(b1)),'-r',...
    b2,sig2fit(1+length(b1):end),'ob',b2,sig_fit(1+length(b1):end),'-b')


% check some bad fits
ind = find(error_list > thresh_error);
[~,sort_ind] = sort(error_list(ind));

if n == 0
    ind = ind(sort_ind(end));
else
    ind = ind(sort_ind(n));
end

vox_ind = vox_mask_ind(ind);
Pout = Pout_list(ind,:);
sig_fit = fit_muFA(Pout,b);
sig2fit = squeeze(voxSig2fit(vox_ind,:,:));
sig2fit = sig2fit(:);
figure(2),clf
semilogy(b1,sig2fit(1:length(b1)),'xr',b1,sig_fit(1:length(b1)),'-r',...
    b2,sig2fit(1+length(b1):end),'ob',b2,sig_fit(1+length(b1):end),'-b')

end