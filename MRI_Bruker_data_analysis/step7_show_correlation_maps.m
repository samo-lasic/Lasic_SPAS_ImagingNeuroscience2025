% Display and/or save qualitative correlation maps.
% These are joint contrast maps combining two parameter maps (e.g., TDD and µA) with optional modulation.

% Combine maps such as TDD-muA, TDD-muA-OP, TDD-muA-FA, or TDD-muA-MD.

% `correlate_name` defines two maps to be mixed.
% `modulate_name` defines a map to control brightness or color modulation.
% - If `modulate_name` has multiple entries, their average will be used.
% - If `modulate_name` is empty, no modulation will be applied.

% `prepend_name` and `append_name` can be empty ('') but may be used for cleaner titles.
% If `title_str` is not empty, it is used directly as the title; otherwise, a title is generated automatically:
% title format: cor. correlate_append from correlate_name{1} and correlate_name{2}, mod. by modulate_append from modulate_name{1} and modulate_name{2} (alpha)

% `modulate_name` can be empty (`{}`) for no modulation. If two names are given, their average intensity is used.
% `color_order` specifies two or three elements: R(1), G(2), B(3), or AlphaMask(0), associated with the correlate and modulate maps.




clear all
close all

setup_code_path()

data_opt.append_folder_name = '_out';
data_path_struct = setup_data_path(data_opt);

warning ('off','all');

% ------------------
% cases:
% 1 - ex vivo (manual rois)
% 2 - ex vivo (automatic rois)
% 3 - in vivo
% 4 - in vivo (automatic rois)

c = 3; % cases

% select slices and tiling 
opt.fig = SPAS_tiling_fig_opt(3); % 4 slices (2-5)
% opt.fig = SPAS_tiling_fig_opt(4); % 4 slices (1-4)
% opt.fig = SPAS_tiling_fig_opt(5); % 6 slices (1-6)
opt.fig = SPAS_tiling_fig_opt(6); % 4 slices in line (2-5)
opt.fig = SPAS_tiling_fig_opt(7); % 6 slices in line (1-6)

opt.fig.save = 1;
opt.fig.num = 10;
opt.fig.show_rois = 0; % if not present, don't show roi

opt.fig.show_roi_scatter = 0; % if not present, don't show roi
opt.fig.show_roi_contours = 1;
opt.fig.remove_roi_outliers = 1; % remove pixels with negative uA or TDD contrast (some in vivo pixels have negative TDD)
opt.fig.show_colormap = 1;

% ---------------------------------------------------
% preprocessing label prefix 

opt.map_name_prepend = '_s_a_g'; % common to all maps
% opt.map_name_prepend = '_mc_s_a_g'; % common to all maps
opt.map_name_prepend = '_a_g'; % common to all maps
% opt.map_name_prepend = '_mc_a_g'; % common to all maps
opt.map_name_prepend = '_dn_mc_a_g'; % common to all maps
opt.map_name_prepend = '_dn_mc_s_a_g'; % common to all maps




% ---------- some defaults ----------------
correlate_prepend = 'pa_raw_'; % e.g. 'pa_raw'
modulate_prepend = correlate_prepend;

% map range prior to normalization to range 0-1
correlate_range.relative = 1; % if 1 use relative limits to [min max] range


if c == 1 || c == 2
    % ex vivo
    correlate_range.LIM1 = [0 .5]; % .6 auto-limit if empty or NaN, e.g. [] or [0 NaN]
    correlate_range.LIM2 = [0 .4]; % .4 auto-limit if empty or NaN, e.g. [] or [0 NaN]
elseif c == 3 || c == 4
    % in vivo
    correlate_range.LIM1 = [0 .65];%[0 .65]; % .6 auto-limit if empty or NaN, e.g. [] or [0 NaN]
    correlate_range.LIM2 = [0 .3];%[0 .3]; % .4 auto-limit if empty or NaN, e.g. [] or [0 NaN]
end

modulate_range.relative = 1; % if 1 use relative limits to [min max] range
modulate_range.LIM = [0 .3]; % .3 2 auto-limit if empty or NaN, e.g. [] or [0 NaN]

title_str = ''; % custom title for correlation map 

color_order = [3 1 0]; % [cor1, cor2, mod] 1-red, 2-green, 3-blue, mod can be 0-alpha or any positive number

correlate_gamma_cor.type = 0; % 3 %%% 0 no correction only A scale , 1. x^p [p > 0], 2. sin(x)^p [p > 0], 3. tanh(x/p)/tanh(1/p) [0 < p < 1]
correlate_gamma_cor.p = .7; % "gamma correction", shape factor (>0, 1 for no effect)
correlate_gamma_cor.A = 1; % "gamma correction", amplification factor

modulate_gamma_cor.type = 0; % 3 %%% 0. no correction only A scale, 1. x^p [p > 0], 2. sin(x)^p [p > 0], 3. tanh(x/p)/tanh(1/p) [0 < p < 1]
modulate_gamma_cor.p = .3; % .3 "gamma correction", shape factor (>0, 1 for no effect)
modulate_gamma_cor.A = 1; % "gamma correction", amplification factor

% ----------

correlate_maps_struct = SPAS_make_correlate_maps_struct({}, ...
    correlate_prepend, {'muA','TDD'}, 'dlogs', correlate_range, correlate_gamma_cor, ...
    modulate_prepend, {}, 'smax', modulate_range, modulate_gamma_cor, color_order, title_str);

correlate_maps_struct = SPAS_make_correlate_maps_struct({}, ...
    correlate_prepend, {'muA','TDD'}, 'dlogs', correlate_range, correlate_gamma_cor, ...
    modulate_prepend, {'muA','TDD'}, 'smax', modulate_range, modulate_gamma_cor, color_order, title_str);

% ----------

% Alternative correlation map configurations

% correlate_maps_struct = SPAS_make_correlate_maps_struct({}, ...
%     correlate_prepend, {'muA','TDD'}, 'dlogs', correlate_range, correlate_gamma_cor, ...
%     'dti', {''}, 'fa_geoSPAS', modulate_range, modulate_gamma_cor, color_order, title_str);

% correlate_maps_struct = SPAS_make_correlate_maps_struct(correlate_maps_struct, ...
%     correlate_prepend, {'muA','TDD'}, 'dlogs', correlate_range, correlate_gamma_cor, ...
%     '', {''}, 'OP', modulate_range, modulate_gamma_cor, color_order, title_str);

%     correlate_maps_struct = SPAS_make_correlate_maps_struct(correlate_maps_struct, ...
%         correlate_prepend, {'muA','TDD'}, 'dlogs', correlate_range, correlate_gamma_cor, ...
%         '', {}, '', modulate_range, modulate_gamma_cor, color_order, title_str);

% ------------

% correlate_range.LIM2 = [0 .5];
% correlate_maps_struct = SPAS_make_correlate_maps_struct({}, ...
%     '', {'pa_raw_muA_ds','fa_geoSPAS'}, '', correlate_range, correlate_gamma_cor, ...
%     modulate_prepend, {}, '', modulate_range, modulate_gamma_cor, color_order, title_str);

% correlate_range.LIM1 = [0 1]; correlate_range.LIM2 = [0 1];
% correlate_maps_struct = SPAS_make_correlate_maps_struct({}, ...
%     '', {'muFA','fa_geoSPAS'}, '', correlate_range, correlate_gamma_cor, ...
%     modulate_prepend, {}, '', modulate_range, modulate_gamma_cor, color_order, title_str);

% correlate_range.LIM1 = [0 1]; correlate_range.LIM2 = [0 1];
% correlate_maps_struct = SPAS_make_correlate_maps_struct({}, ...
%     '', {'muFA','OP'}, '', correlate_range, correlate_gamma_cor, ...
%     modulate_prepend, {}, '', modulate_range, modulate_gamma_cor, color_order, title_str);


% correlate_range.LIM1 = [0 1]; correlate_range.LIM2 = [0 1];
% correlate_maps_struct = SPAS_make_correlate_maps_struct({}, ...
%     '', {'OP','fa_geoSPAS'}, '', correlate_range, correlate_gamma_cor, ...
%     modulate_prepend, {}, '', modulate_range, modulate_gamma_cor, color_order, title_str);


% threshold for data zooming (eps = smallest positive value)
opt.data_zoom_thresh = eps; %0.005; % .015 signal projection is larger than threshold

% ------------------
opt.save_bin_rois = 1;

% define bins (relative to full range)
if c == 1 || c == 2
    crossover1 = 0.45;
    crossover2 = 0.2;
elseif c == 3 || c == 4
    crossover1 = 0.45;
    crossover2 = 0.2;
end

opt.bins = {};
opt.bins(end+1).range1 = [0 crossover1];
opt.bins(end).range2 = [crossover2 1];
opt.bins(end).range3 = [];

opt.bins(end+1).range1 = [crossover1 1];
opt.bins(end).range2 = [0 crossover2];
opt.bins(end).range3 = [];



% ------------------

opt.fig.roi_plot_order = 1:numel(data_path_struct.roi_names); % plot ROIs in original input order

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

    % axis limits for scatter plots (remove field or use [] for auto-scaling)
    opt.fig.XLIM = [0.0441913, 1.15639]; 
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

    % generate and save correlation maps for current dataset
    SPAS_make_correlation_maps(root_data_path, merged_names, roi_names, correlate_maps_struct, opt);

end





