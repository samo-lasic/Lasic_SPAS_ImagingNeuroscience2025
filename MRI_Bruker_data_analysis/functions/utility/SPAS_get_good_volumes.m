function good_volumes = SPAS_get_good_volumes(nii_fn, mask_nii_fn, thresh, do_plot)
% uses mask to get relative signal differences within and outside mask
% thresh represents number of standard deviations for largest data with largest b value
if nargin < 4
    do_plot = 0;
end

roi = mdm_nii_read(mask_nii_fn);
not_roi = ~roi;

% read nifti
nii = mdm_nii_read(nii_fn);
xps = mdm_xps_load(mdm_fn_nii2xps(nii_fn));

mean_s = nii.*single(roi);
mean_s = mean(mean(mean(mean_s)));
mean_s = squeeze(mean_s);

mean_not_s = nii.*single(not_roi);
mean_not_s = mean(mean(mean(mean_not_s)));
mean_not_s = squeeze(mean_not_s);

rel_dif = (mean_s-mean_not_s)./(mean_s+mean_not_s);%tot_mean_sig;

% figure, hold on
% plot(rel_dif,'-bo')
% plot(mean_s,'-bo')
% plot(mean_not_s,'-rx')

ind = find(xps.a_ind == max(xps.a_ind)); % find min b
d = rel_dif(ind);
md = mean(d);
stdd = std(d);
ind = ind(find(abs((d-md)/stdd) > thresh));

n_volumes = length(mean_s);
vol_num = 1:n_volumes;


good_volumes = ones(n_volumes,1);
good_volumes(ind) = 0;

if do_plot

    n_rejected = length(ind);
    disp_str = sprintf('%d volumes out of %d may be rejected', n_rejected,n_volumes);
    display(sprintf('%s',disp_str));

    figure
    hold on
    plot(vol_num,rel_dif,'-bo')
    plot(vol_num(ind),rel_dif(ind),'rx')
    title(disp_str)
end

end

function [step, thresh] = t_step(N, rel_dif)
t = linspace(min(rel_dif),max(rel_dif),N);
for n = 1:N
    if numel(find(rel_dif < t(n))) == 1; 
        step = n/N;
        thresh = t(n);
        continue
    end
end
end