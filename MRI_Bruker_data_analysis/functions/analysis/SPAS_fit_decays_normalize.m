function out = SPAS_fit_decays_normalize(nii_fn,mask_nii_fn, opt)
% used to fit decays and normalize data (output also gamma fit results and fit decays)
% expects averaged repetitions

% outputs
% out.norm_nii_fn
% out.norm_fit_nii_fn
% out.s0_nii_fn
% out.adc_nii_fn
% out.mu2_nii_fn
% out.norm_par_mat_fn

% read image
xps = mdm_xps_load(mdm_fn_nii2xps(nii_fn));
[sig, h_nii] = mdm_nii_read(nii_fn);

% read mask
[mask, h_mask] = mdm_nii_read(mask_nii_fn);

Nwfm = max(xps.s_ind);

% collect decay indices for each waveform / rotation / repetition 
decay = {};
cnt = 0;
for n_wfm = 1:Nwfm
    for n_rot = unique(xps.rot_ind)'
        b_ind = find(xps.rot_ind == n_rot & xps.s_ind == n_wfm);
        cnt = cnt +1;
        decay(cnt).b_ind = b_ind;
        decay(cnt).len = length(b_ind);
        decay(cnt).n_rot = n_rot;
        decay(cnt).n_wfm = n_wfm;
    end
end
% remove single rotations
decay = decay(find([decay.len] ~= 1));
Ndecays = numel(decay);


% reshape signal into 2D
sig_sz = size(sig);
n_voxels = prod(sig_sz(1:3));
sig = reshape(sig,n_voxels,sig_sz(4));

% ensure data_is_normalized = 1 because we will normalize based on initial slope first 
opt.data_is_normalized = 1;


if opt.weighted
    fit_opt = fit_gamma_weighted_opt(opt);
else
    fit_opt = fit_gamma_opt(opt);
end

% fit_opt.lsqcurvefit.Display = 'on';

% b-range limit for initial slope
fit_opt.b_low_lim = opt.b_low_lim;

vox_mask_ind = find(mask);
n_voxels2fit = length(vox_mask_ind);

gamma_fit_parameters = zeros(n_voxels2fit,Ndecays,3);

tic
parfor n_vox_mask_ind = 1:n_voxels2fit
%    for n_vox_mask_ind = 1:n_voxels2fit
    vox_ind = vox_mask_ind(n_vox_mask_ind);
    
    for n_decay = 1:Ndecays
        
        b_ind = decay(n_decay).b_ind;
        sig1 = double(sig(vox_ind,b_ind)');
        
        % signal must be positive
        if isempty(find(sig1<=0))
            
            b1 = double(xps.b(b_ind));
            
            % quick fit / check
            b_low_ind = find(b1 < fit_opt.b_low_lim);
            if isempty(b_low_ind)
                b_low_ind = [1 2]';
            end
            sig_low = sig1(b_low_ind);
            b_low = b1(b_low_ind);
            
            X = [ones(length(b_low),1) -b_low];
            P = X\log(sig_low);
            
            %         figure(1),clf
            %         semilogy(b_low,sig_low,'-o')
            %
            %         sig_fit = exp(X*P);
            %         figure(2),clf
            %         semilogy(b_low,sig_low,'o',b_low,sig_fit,'-')
            
            if P(2) > 0
                s0 = exp(P(1));

                if opt.weighted
                    Pout = fit_gamma1_weighted_par(b1,sig1/s0, fit_opt);
                else
                    Pout = fit_gamma1_par(b1,sig1/s0, fit_opt);
                end
                               
                gamma_fit_parameters(n_vox_mask_ind,n_decay,:) = [s0*Pout(1) Pout(2:3)];
                
            end
        end
    end
end
toc

% rearange signals and maps
sig_norm = 0*sig;
if opt.save.gamma_fit_signal
    sig_norm_fit = sig_norm;
end

if opt.save.s0_map
    s0_map = zeros(n_voxels,Ndecays);
end
if opt.save.adc_map
    adc_map = zeros(n_voxels,Ndecays);
end
if opt.save.mu2_map
    mu2_map = zeros(n_voxels,Ndecays);
end

% save normalized signal and maps

for n_vox_mask_ind = 1:n_voxels2fit
    vox_ind = vox_mask_ind(n_vox_mask_ind);
    
    for n_decay = 1:Ndecays
        b_ind = decay(n_decay).b_ind;
        gamma_fit_p = squeeze(gamma_fit_parameters(n_vox_mask_ind,n_decay,:));
        if prod(gamma_fit_p(1:2) ~= 0) % if fit was ok (s0, ADC)
            
            sig_norm(vox_ind,b_ind) = sig(vox_ind,b_ind)/gamma_fit_p(1);
            
            if opt.save.gamma_fit_signal
                b1 = xps.b(b_ind);
                sig_norm_fit(vox_ind,b_ind) = fit_gamma1([1; gamma_fit_p(2:end)],b1);
                %
                %                 figure(1),clf
                %                 semilogy(b1,sig_norm(vox_ind,b_ind),'o',b1,sig_norm_fit(vox_ind,b_ind),'-')
            end
            
            if exist('s0_map')
                s0_map(vox_ind,n_decay) = gamma_fit_p(1);
            end
            if exist('adc_map')
                adc_map(vox_ind,n_decay) = gamma_fit_p(2);
            end
            if exist('mu2_map')
                mu2_map(vox_ind,n_decay) = gamma_fit_p(3);
            end
            
        end
    end
end

%reshape back to image dimensions
sig_norm = reshape(sig_norm,sig_sz);

if opt.save.gamma_fit_signal
    sig_norm_fit = reshape(sig_norm_fit,sig_sz);
end

% make a subset of xps (Ndecays)
fds_remove = {'n','b','bt','rot_ax','rota_ax_str','rot_ax','wfm_names', 'wfm_src','g_ind',''};
xps_sub = xps;
for n = 1:length(fds_remove)
    if isfield(xps_sub,fds_remove{n})
        xps_sub = rmfield(xps_sub,fds_remove{n});
    end
end

for n_decay = 1:Ndecays
    sub_ind(n_decay) = decay(n_decay).b_ind(1);
end


fds = fields(xps_sub);
for n = 1:length(fds)
    tmp = xps_sub.(fds{n});
    sz = size(tmp);
    xps_sub.(fds{n}) = reshape(tmp(sub_ind,:),[Ndecays sz(2:end)]);
end
xps_sub.n = Ndecays;

% save normalized signals
out_nii_fn = append_nii_fn(nii_fn, 'norm');
mdm_nii_write(single(sig_norm),out_nii_fn);
mdm_xps_save(xps, mdm_fn_nii2xps(out_nii_fn));
out.norm_nii_fn = out_nii_fn;

% save normalized fitted signals
if opt.save.gamma_fit_signal
    out_nii_fn = append_nii_fn(nii_fn, 'norm_fit');
    mdm_nii_write(single(sig_norm_fit),out_nii_fn);
    mdm_xps_save(xps, mdm_fn_nii2xps(out_nii_fn));
    out.norm_fit_nii_fn = out_nii_fn;
end

% save maps
if exist('s0_map')
    s0_map = single(s0_map);
    s0_map = reshape(s0_map,[sig_sz(1:3) Ndecays]);
    out_nii_fn = append_nii_fn(nii_fn, 's0');
    mdm_nii_write(s0_map,out_nii_fn);
    mdm_xps_save(xps_sub, mdm_fn_nii2xps(out_nii_fn));
    out.s0_nii_fn = out_nii_fn;
end
if exist('adc_map')
    adc_map = single(adc_map);
    adc_map = reshape(adc_map,[sig_sz(1:3) Ndecays]);
    out_nii_fn = append_nii_fn(nii_fn, 'adc');
    mdm_nii_write(adc_map,out_nii_fn);
    mdm_xps_save(xps_sub, mdm_fn_nii2xps(out_nii_fn));
    out.adc_nii_fn = out_nii_fn;
end
if exist('mu2_map')
    mu2_map = single(mu2_map);
    mu2_map = reshape(mu2_map,[sig_sz(1:3) Ndecays]);
    out_nii_fn = append_nii_fn(nii_fn, 'mu2');
    mdm_nii_write(mu2_map,out_nii_fn);
    mdm_xps_save(xps_sub, mdm_fn_nii2xps(out_nii_fn));
    out.mu2_nii_fn = out_nii_fn;
end

% save fit parameters (mat file)
if opt.save.gamma_fit_parameters
    
    % add vox_mask_ind info
    res.gamma_fits = gamma_fit_parameters;
    res.vox_mask_ind = vox_mask_ind;
    res.img_size = sig_sz(1:3);
    
    % associate decays to rotations and waveforms
    res.decay = decay;
    res.wfm_names = xps.wfm_names;
    
    out_mat_fn = append_mat_fn(nii_fn, 'norm_par');
    save(out_mat_fn,'res')
    
    out.norm_par_mat_fn = out_mat_fn;

end




