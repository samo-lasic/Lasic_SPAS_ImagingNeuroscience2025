% Generates tiled slice images from normalized signal differences and/or log-signal differences at maximum b-value. 
% Maps labeled `_raw_TDD_` and `_raw_muA_` are saved in the `fig` folder. 
% Optionally overlays ROI contours (colored outlines) to highlight selected regions.

clear all
close all

setup_code_path()

data_opt.append_folder_name = '_out';
data_path_struct = setup_data_path(data_opt);

warning ('off','all');

nii_base_str = 'merged_';


% ------------------
% 1 - ex vivo (manual rois)
% 2 - ex vivo (automatic rois)
% 3 - in vivo
% 4 - in vivo (automatic rois)


c = 1;

% select slices and tiling 
opt.fig = SPAS_tiling_fig_opt(3); % 4 slices (2-5)
% opt.fig = SPAS_tiling_fig_opt(4); % 4 slices (1-4)
% opt.fig = SPAS_tiling_fig_opt(5); % 6 slices (1-6)
opt.fig = SPAS_tiling_fig_opt(6); % 4 slices in line (2-5)
% opt.fig = SPAS_tiling_fig_opt(7); % 6 slices in line (1-6)



opt.fig.save = 1;
opt.fig.num = 10;
opt.fig.show_rois = 0; % if not present, don't show roi


if (1) % general - data label prefix
    % opt.map_name_prepend = '_a_g_pa_raw_';
    % opt.map_name_prepend = '_mc_a_g_pa_raw_';
    opt.map_name_prepend = '_dn_mc_a_g_pa_raw_';
    opt.map_name_prepend = '_dn_mc_s_a_g_pa_raw_';

    opt.map_names = {'muA_ds','TDD_ds', 'muA_dlogs', 'TDD_dlogs',};
    opt.title_strings = opt.map_names;
    % selects maps for display
    opt.select_map_numbers = [3 4]; % [3 4]; 

    opt.gamma_cor = 1; % I = I^gamma_cor - if 1 no correction (linear contrast)
end




% ---------------------------------------------------
% opt.map_names
opt.autoMIN_array = 0 * [1 1 1 1]; % to set limits automatically, set this to 1
opt.autoMAX_array = 0 * [1 1 1 1];

if c == 1 || c == 2
    opt.MIN_array = [0.00 0.00 0.00 0.00]; % full brain
    opt.MAX_array = [0.35 0.45 0.35 0.45];

    % opt.MIN_array = [0.00 0.00 0.02 0.02]; % cortex
    % opt.MAX_array = [0.35 0.45 0.2 0.25];

elseif c == 3 || c == 4
    opt.MIN_array = [0.00 0.00 0.00 0.00]; % rat in vivo
%     opt.MAX_array = [0.35 .3 0.7 .3];
    opt.MAX_array = [0.35 .3 1 .25];
end

opt.colormap = 'gray'; % gray inferno magma turbo jet cividis viridis plasma

opt.data_zoom_thresh = 0; % .005 .015 signal projection is larger than threshold

opt.do_zoom_range = 0; % this is for extra zooming
opt.zoom_range_x = [0 1]; % 0-1 relative to full range
opt.zoom_range_y = [0 1]; % 0-1 relative to full range

% opt.zoom_range_x = [.2 .8]; % 0-1 relative to full range
% opt.zoom_range_y = [0.6 .9]; % 0-1 relative to full range


opt.modulate_intensity.type = 0; % 0 - no modulation, 1 - intensity (not color) is modulated by smin, 2 - intensity (not color) is modulated by smax

opt.modulate_intensity.thresh = 0;
opt.modulate_intensity.gamma_cor.type = 3; % 3 % 0 - only A-scale, 1. x^p [p > 0], 2. sin(x)^p [p > 0], 3. tanh(x/p)/tanh(1/p) [0 < p < 1]
opt.modulate_intensity.gamma_cor.p = .3; %.3; % "gamma correction", shape factor (>0, 1 for no effect)
opt.modulate_intensity.gamma_cor.A = 1; % "gamma correction", amplification factor

% ------------------

if c == 1 || c == 2 % ex vivo
    opt.fig.roi_colors = viridis(numel(data_path_struct.roi_names));
    opt.fig.roi_alpha = 0.75 * ones(1,numel(data_path_struct.roi_names));

    opt.fig.roi_detail = 1; % relative to 1
    opt.fig.roi_smoothing = 2; % 2
end

if c == 3 || c == 4 % in vivo
    
    opt.fig.roi_colors = [];
    Nrois = numel(data_path_struct.roi_names);
    if Nrois == 1
        opt.fig.roi_colors = [1 1 1];
    elseif Nrois < 5
        opt.fig.roi_colors = viridis(numel(data_path_struct.roi_names));
    elseif Nrois == 5
        opt.fig.roi_colors = viridis(numel(data_path_struct.roi_names)-1);
        opt.fig.roi_colors(end+1,:) = [1 1 1];
    end

    
    %     opt.fig.roi_alpha = 0.75 * ones(1,numel(data_path_struct.roi_names));
    %     opt.fig.roi_alpha(end) = .3;
    opt.fig.roi_alpha = [1 0.75 0.75 0.75 0.3];
    opt.fig.roi_plot_order = [1 4 2 3 5];

    opt.fig.roi_detail = 1.5; % relative to 1
    opt.fig.roi_smoothing = 1.5; % 2

    opt.fig.XLIM = [0.0441913, 1.15639]; % automatic if this field is not existing or empty []
    opt.fig.YLIM = [0.000554465, 0.320087];
end

if c == 2 || c == 4
    opt.fig.roi_colors = [1 0.4 0; 0 .8 1];
    opt.fig.roi_alpha = 1*[1 1];
end


for n_data_path = 1:numel(data_path_struct)
    root_data_path = data_path_struct(n_data_path).root_data_path;
    exp_folder_name = data_path_struct(n_data_path).exp_folder_name;
    root_data_path = fullfile(root_data_path,exp_folder_name);
    select_subfolders = data_path_struct(n_data_path).select_subfolders;
    roi_names = data_path_struct(n_data_path).roi_names;
    roi_names = ensure_fn_ext(roi_names, '.nii.gz');

    if isfield(data_path_struct(n_data_path),'post_mask')
        opt.post_mask = data_path_struct(n_data_path).post_mask;
        opt.post_mask = ensure_fn_ext(opt.post_mask, '.nii.gz');
    end

    merged_names = SPAS_numbers_to_names(select_subfolders, nii_base_str);

    % generate and save normalized/log signal difference maps
    SPAS_make_sig_dif_maps(root_data_path, merged_names, roi_names, opt);

end






