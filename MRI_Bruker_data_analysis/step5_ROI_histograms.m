
% Perform ROI statistics using NIfTI maps, such as those from DTI non-linear fits or raw signal difference maps (norm_pa_raw_ds_TDD).
% Generate overlap histograms for experiments specified in the select_subfolders array.
% This script is useful for comparing ROI data within the same session.
% For comparisons across sessions, modifications to the script would be required.
% Saves ROI statistics in 'stats_all.mat' and histogram figures in the 'fig/' subfolder.
% Designed to be used together with step6_ROI_stats.m to generate summary figures and explore ROI statistics from multiple experiments.

clear all
close all

setup_code_path()

data_opt.append_folder_name = '_out';
data_path_struct = setup_data_path(data_opt);

warning ('off','all');


% ------------------------------------------

auto_find_rois = 0; % if 0 read from the data_path_struct
find_roi_name_pattern = 'roi_'; % string pattern(s) for roi files, could be an array {'AAA','BBB',...}

stat_opt.CI = [68, 85, 90, 95]; % 68 [confidence interval(s)]
out_stats_path = '...';

fig_opt.save = 0;
fig_opt.base_fig_name = 'hist';
fig_opt.resolution = '-r200';


% prefixes reflecting different data processing stages
% map_name_prepend = '_s_g_';
% map_name_prepend = '_g_';
% map_name_prepend = '_dn_g_';
map_name_prepend = '_dn_mc_g_';

% Map names for ROI extraction and histogram plotting
map_names = {'dti_s0_geoSPAS','dti_md_geoSPAS','dti_fa_geoSPAS',...
    'pa_raw_TDD_dlogs'};

% map_names = {'dti_s0_geoSPAS','dti_md_geoSPAS','dti_fa_geoSPAS',...
%     'norm_pa_raw_TDD_dlogs'};


title_strings = {'dti_s0_geoSPAS','dti_md_geoSPAS','dti_fa_geoSPAS',...
    'TDD dlogs'};



labels = {'s0','MD [mm^2/ms]','FA','TDD signal log dif.'};


BinWidths   = [100,     0.02,   0.05,   .01];
scales      = [1,       1e-9,   1,        1];
Xmins       = [1e3,     0.3,   0,     -.15];
Xmaxs       = [2.5e3,   1.1,   1,       .2];

sel_maps = [2:4]; %2 3 4

map_names = map_names(sel_maps);
title_strings = title_strings(sel_maps);
labels = labels(sel_maps);
BinWidths = BinWidths(sel_maps);
scales = scales(sel_maps);
Xmins = Xmins(sel_maps);
Xmaxs = Xmaxs(sel_maps);



for n_data_path = 1:numel(data_path_struct)
    root_data_path = data_path_struct(n_data_path).root_data_path;
    exp_folder_name = data_path_struct(n_data_path).exp_folder_name;
    root_data_path = fullfile(root_data_path,exp_folder_name);
    select_subfolders = data_path_struct(n_data_path).select_subfolders;

    % find roi_names
    if auto_find_rois
        roi_names = dir(fullfile(root_data_path));
        roi_names = {roi_names.name};

        if iscell(find_roi_name_pattern)
            ind = contains(roi_names, find_roi_name_pattern{1});
            for n_str = 2:numel(find_roi_name_pattern)
                ind = ind & contains(roi_names, find_roi_name_pattern{n_str});
            end
        else
            ind = contains(roi_names, find_roi_name_pattern);
        end
        roi_names = roi_names(ind);

    else
        roi_names = data_path_struct(n_data_path).roi_names;
    end


    merged_names = SPAS_numbers_to_names(select_subfolders, 'merged_');

    
    % structure s contains statistical values arranged in a 3D array : map x roi x exp-subfolders
    tmp = make_statistics(root_data_path, exp_folder_name, merged_names, roi_names, map_names, map_name_prepend, stat_opt);
    
    % copy fields of s
    for fn = fieldnames(tmp)'
        s(n_data_path).(fn{1}) = tmp.(fn{1});
    end
    s(n_data_path).roi_names = roi_names;
    s(n_data_path).exp_folder_name = exp_folder_name;
    s(n_data_path).select_subfolders = select_subfolders;
    
    stats.map_names = map_names;
    stats.CI = stat_opt.CI;
    stats.s = s;

    % save stats
    mkdir(out_stats_path)
    save(fullfile(out_stats_path,'stats_all.mat'), 'stats')
    
   
    make_histograms(root_data_path, exp_folder_name, merged_names, roi_names, map_names, map_name_prepend, ...
        title_strings, labels, BinWidths, scales, Xmins, Xmaxs, fig_opt)
 
end




% ------------------------- FUNCTIONS -------------------------------------

function make_histograms(root_data_path, exp_folder_name, merged_names, roi_names, map_names, map_name_prepend, ...
    title_strings, labels, BinWidths, scales, Xmins, Xmaxs, fig_opt)


for n_roi = 1:numel(roi_names)
    roi_name = roi_names{n_roi};

    % roi path
    roi_fn = fullfile(root_data_path,roi_name);
    roi = mdm_nii_read(roi_fn);
    ind = find(roi);

    for n_map = 1:numel(map_names)
        map_name = [map_name_prepend map_names{n_map} '.nii.gz'];
      

        % --------  histograms ------------


        title_str = title_strings{n_map};
        label = labels{n_map};
        BinWidth = BinWidths(n_map);
        scale = scales(n_map);
        XLIM = [Xmins(n_map) Xmaxs(n_map)];

        %col = copper(numel(merged_names));
        %col = cividis(numel(merged_names));
        col = viridis(numel(merged_names));


        clear h l;

        fh = figure;
        clf
        hold on
        for n_exp = 1:numel(merged_names)
            merged_name = merged_names{n_exp};

            name = [merged_name map_name];
            nii_fn = fullfile(root_data_path,name);
            nii = mdm_nii_read(nii_fn);
            roi_sig = nii(ind);

            roi_sig = nii(ind)/scale;
            h(n_exp) = histogram(roi_sig);
            h(n_exp).Normalization = 'probability';
            h(n_exp).BinWidth = BinWidth;
            h(n_exp).EdgeAlpha = 0;
            h(n_exp).FaceColor = col(n_exp,:);
            h(n_exp).FaceAlpha = 0.2;

            x = (h(n_exp).BinEdges(1:end-1)+ h(n_exp).BinEdges(2:end))/2;
            y = h(n_exp).Values;

            if length(x) > 2
                x_interp = linspace(min(x),max(x),100);
                y_interp = interp1(x,y,x_interp,'pchip');
            else
                x_interp = x;
                y_interp = y;
            end
            l(n_exp) = plot(x_interp,y_interp,'-','color',h(n_exp).FaceColor,'LineWidth',3);
        end

        legend_str = strrep(merged_names,'_','-');
        lh = legend(l,legend_str);
        lh.Box = 'off';
        title([extractBefore(exp_folder_name,'_') ': ' title_str ' - ' extractBefore(roi_name,'.nii.gz')],'interpreter','none')
        if ~isempty(XLIM)
            xlim(XLIM)
        end
        xlabel(label);

        set(gca,'LineWidth',3,'Yscale','Lin','Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',14)
        pbaspect(gca,[1.62 1 1])

        if fig_opt.save
            fig_path = fullfile(root_data_path, 'fig');
            mkdir(fig_path)
            figName = [extractBefore(exp_folder_name,'_') '_' merged_name extractBefore(map_name,'.nii.gz') '_' ...
                extractBefore(roi_name,'.nii.gz') '_' fig_opt.base_fig_name '.png' ];

            fig_path = fullfile(fig_path, figName);
            display(sprintf('saving ... %s',fig_path))
            print(fh, fig_path, '-dpng', fig_opt.resolution);

        end


    end

end
end


function s = make_statistics(root_data_path, exp_folder_name, merged_names, roi_names, map_names, map_name_prepend, stat_opt)
% Returns structure s: statistics as 3D array (map × ROI × experiment)
% Compute mean, std, and CI bounds per ROI, per map, per experiment.


for n_roi = 1:numel(roi_names)
    roi_name = roi_names{n_roi};

    % roi path
    roi_fn = fullfile(root_data_path,roi_name);
    roi = mdm_nii_read(roi_fn);
    ind = find(roi);

    s.roi_size(n_roi) = numel(ind);

    for n_map = 1:numel(map_names)
        map_name = [map_name_prepend map_names{n_map} '.nii.gz'];

        for n_exp = 1:numel(merged_names)
            merged_name = merged_names{n_exp};
            name = [merged_name map_name];
            nii_fn = fullfile(root_data_path,name);
            nii = mdm_nii_read(nii_fn);
            roi_sig = nii(ind);
            roi_mean = mean(roi_sig);
            roi_std = std(roi_sig);

            [roi_lb, roi_ub] = get_roi_CI_bounds(roi_sig, roi_mean, stat_opt.CI);

            
            % copy to output struct
            s.roi_mean(n_map, n_roi, n_exp) = roi_mean;
            s.roi_std(n_map, n_roi, n_exp) = roi_std;
            s.roi_lb(n_map, n_roi, n_exp,:) = roi_lb;
            s.roi_ub(n_map, n_roi, n_exp,:) = roi_ub;
            

        end
    end
end
end


function [roi_lb, roi_ub] = get_roi_CI_bounds(roi_sig, roi_mean, CI)
% Compute CI bounds by selecting values closest to mean 

tmp = abs(roi_sig - roi_mean);

[~,ind_tmp] = sort(tmp);
for n = 1:length(CI)
    ind_tmp = ind_tmp(1:round(length(ind_tmp) * CI(n)/100));

    % figure(1),clf
    % hold on
    % plot(roi_sig,'.')
    % plot(ind_tmp,roi_sig(ind_tmp),'or')

    roi_lb(n) = min(roi_sig(ind_tmp));
    roi_ub(n) = max(roi_sig(ind_tmp));
end
end
