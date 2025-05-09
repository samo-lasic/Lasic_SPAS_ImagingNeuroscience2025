% Performs ROI statistics on NIfTI maps, e.g., from DTI non-linear fits or raw signal difference maps (norm_pa_raw_ds_TDD).
% Generates overlap histograms for experiments defined in the select_subfolders array.
% Useful for comparing ROI data within the same session or across sessions.


clear all
close all


setup_code_path()

warning ('off','all');

data_opt.append_folder_name = '_out';
data_path_struct = setup_data_path(data_opt);


% ------------------------------------------

% global options
select_map_names = {'md'}; % {} if empty, select all
select_roi_names = {'test1'}; % {} if empty, select all unique for stats or all present for histograms

auto_find_rois = 1; % if 0 read from the data_path_struct
find_roi_name_pattern = 'roi_'; %'ROI_'; % string pattern(s) for roi files, could be an array {'AAA','BBB',...}


% histogram options
hist_opt.plot = 0;
hist_opt.save = 0;
hist_opt.base_fig_name = 'hist';
hist_opt.resolution = '-r200';

% stats options
stat_opt.save = 1;
stat_opt.load = 0;
stat_opt.CI = 68; % 90 95 68 confidence interval
stat_opt.out_path = '...'; % specify output path for saving statistics

stat_opt.save_fig = 0;
stat_opt.resolution = hist_opt.resolution;

show_stats = 1; % show stats per map and roi
show_exp_dif = 0; % show mean difference per map and roi


% map_name_prepend = '_s_g_';
% map_name_prepend = '_g_';
% map_name_prepend = '_dn_g_';
map_name_prepend = '_dn_mc_g_';

% Define input map names
map_names = {...
    'dti_s0_geoSPAS','dti_md_geoSPAS','dti_fa_geoSPAS',...
    'dti_s0_SPAS1','dti_md_SPAS1','dti_fa_SPAS1',...
    'dti_s0_SPAS2','dti_md_SPAS2','dti_fa_SPAS2',...
    'dti_s0_SPAS3','dti_md_SPAS3','dti_fa_SPAS3',...
    'dti_dif_md',...
    'pa_raw_TDD_dlogs'};

title_strings = {...
    'dti_s0_geoSPAS','dti_md_geoSPAS','dti_fa_geoSPAS',...
    'dti_s0_SPAS1','dti_md_SPAS1','dti_fa_SPAS1',...
    'dti_s0_SPAS2','dti_md_SPAS2','dti_fa_SPAS2',...
    'dti_s0_SPAS3','dti_md_SPAS3','dti_fa_SPAS3',...
    'dti_dif_md',...
    'TDD dlogs'};

labels = {...
    's0','MD [mm^2/ms]','FA',...
    's0(SPAS1)','MD(SPAS1) [mm^2/ms]','FA(SPAS1)',...
    's0(SPAS2)','MD(SPAS2) [mm^2/ms]','FA(SPAS2)',...
    's0(SPAS3)','MD(SPAS3) [mm^2/ms]','FA(SPAS3)',...
    'dMD (SPAS3-SPAS1) [mm^2/ms]',...
    'TDD signal log dif.'};

% ----------   s0       MD      FA      dMD      TDD signal log dif.
BinWidths   = [100,     0.02,   0.05,   0.02       .01];
BinWidths = [BinWidths(1:3) BinWidths(1:3) BinWidths(1:3) BinWidths(1:3) BinWidths(4) BinWidths(5)];
scales      = [1,       1e-9,   1,       1e-9,        1];
scales = [scales(1:3) scales(1:3) scales(1:3) scales(1:3) scales(4) scales(5)];
Xmins       = [1e3,     0.3,   0,       -0.2         -.15]; % used for histograms
Xmins = [Xmins(1:3) Xmins(1:3) Xmins(1:3) Xmins(1:3) Xmins(4) Xmins(5)];
Xmaxs      = [2.5e3,   1.1,   1,    0.5,    .2];
Xmaxs = [Xmaxs(1:3) Xmaxs(1:3) Xmaxs(1:3) Xmaxs(1:3) Xmaxs(4)  Xmaxs(5)];


% reduce full set of maps
% sel_maps = [2 3 13];
sel_maps = [2 3 5 6 11 12 13 14]; % includes SPAS1, SPAS3, and dMD maps


map_names = map_names(sel_maps);
title_strings = title_strings(sel_maps);
labels = labels(sel_maps);
BinWidths = BinWidths(sel_maps);
scales = scales(sel_maps);
Xmins = Xmins(sel_maps);
Xmaxs = Xmaxs(sel_maps);


% select maps to plot (for histograms and statistics)
if ~isempty(select_map_names)
    map_names_plot = map_names(contains(map_names, select_map_names));
end


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

    % ------------------- statistics -----------
    if stat_opt.save
        % compute ROI statistics for selected maps and ROIs
        % Structure contains statistical values arranged in a 3D array: map × ROI × experiment subfolders.
        tmp = make_statistics(root_data_path, exp_folder_name, merged_names, roi_names, map_names, map_name_prepend, stat_opt);

        % copy fields of s
        for fieldname = fieldnames(tmp)'
            s(n_data_path).(fieldname{1}) = tmp.(fieldname{1});
        end
        s(n_data_path).roi_names = roi_names;
        s(n_data_path).exp_folder_name = exp_folder_name;
        s(n_data_path).select_subfolders = select_subfolders;
    end

    % ---------------- histograms ---------------
    if hist_opt.plot
        make_histograms(root_data_path, exp_folder_name, merged_names, roi_names, select_roi_names, map_names, select_map_names, map_name_prepend, ...
            title_strings, labels, BinWidths, scales, Xmins, Xmaxs, hist_opt)
    end

end

if stat_opt.save
    stats.map_names = map_names;
    stats.CI = stat_opt.CI;
    stats.s = s;

    % save stats
    mkdir(stat_opt.out_path)
    fn = fullfile(stat_opt.out_path,'stat.mat');
    save(fn, 'stats')
    display(sprintf('stats saved in: %s', fn))
end


if stat_opt.load
    close all

    % load stats
    fn = fullfile(stat_opt.out_path,'stat.mat');
    load(fn, 'stats')
    display(sprintf('stats loaded from: %s', fn))

    treated_names = {'xxx','yyy'}; % specify experiment folder names for exclusion or inclusion

    if (1)
        include_exp_folder_name = {}; % if {} include all
        exclude_exp_folder_name = treated_names; %{}; % if {} include all
        fig_folder_name = 'fig_control';
    end
    if (0)
        include_exp_folder_name = treated_names; % if {} include all
        exclude_exp_folder_name = {}; %{}; % if {} include all
        fig_folder_name = 'fig_selected';
    end
    if(1)
        include_exp_folder_name = {}; % if {} include all
        exclude_exp_folder_name = {}; %{}; % if {} include all
        fig_folder_name = 'fig_all';
    end
    find_global_lim = 1; % if 1 find global limits

    lw1 = 2;
    lw2 = 3;
    fs = 16;



    % --------------------

    if stat_opt.save_fig
        fig_folder = fullfile(stat_opt.out_path,fig_folder_name);
        mkdir(fig_folder)
    end

    roi_names = {};
    for n_exp = 1:numel(stats.s)
        roi_names = [roi_names, stats.s(n_exp).roi_names];
    end

    char_ind = 2;
    for n_roi = 1:numel(roi_names)
        roi_name = roi_names{n_roi};
        roi_name = extractBefore(roi_name,'.nii.gz');
        char_count = strfind(roi_name, '_') + 1;
        if length(char_count) > char_ind-1
            roi_name = roi_name(char_count(char_ind):end);
        else
            roi_name = '';
        end
        roi_names{n_roi} = roi_name;
    end

    roi_names = roi_names(~cellfun('isempty', roi_names));
    [~,uind] = unique(lower(roi_names));
    roi_names = roi_names(uind);

    if ~isempty(select_roi_names)
        roi_names = roi_names(contains(roi_names, select_roi_names));
    end

    % ensure roi_name is present in all datasets
    for n_roi = 1:numel(roi_names)
        roi_ind(n_roi) = all(arrayfun(@(x) any(contains(x.roi_names,roi_names(n_roi))), stats.s));
        %roi_ind(n_roi) = all(arrayfun(@(x) any(contains(lower(x.roi_names),lower(roi_names(n_roi)))), stats.s));
    end
    roi_names = roi_names(roi_ind)


    map_names = stats.map_names;
    if ~isempty(select_map_names)
        map_names = map_names(contains(map_names,select_map_names));
    end

    % get maximum array length for color coding
    max_len_subfolders = max(arrayfun(@(x) length(x.select_subfolders), stats.s));
    col = viridis(max_len_subfolders);



    % ---- show stats per map and roi ----
    if (show_stats)
        for n_map = 1:numel(map_names)
            map_name = map_names{n_map};
            map_ind = find(contains(stats.map_names,map_name));

            fig = [];

            for n_roi = 1:numel(roi_names)
                roi_name = roi_names{n_roi};

                %title_str = sprintf('%s - %s',map_name, roi_name);
                title_str = sprintf('%s (CI = %g %s)', roi_name, stat_opt.CI,'%');

                % ----- new figure ----
                fig(end+1).h = figure;
                hold on

                c = 0;
                xtlab = {''};


                for n_exp = 1:numel(stats.s)

                    % Structure s contains statistical values arranged in a 3D array: map × ROI × experiment subfolders.

                    s = stats.s(n_exp);

                    if isempty(include_exp_folder_name)
                        if contains(s.exp_folder_name, exclude_exp_folder_name)
                            continue
                        end
                    else
                        if ~(contains(s.exp_folder_name, include_exp_folder_name) & ~contains(s.exp_folder_name, exclude_exp_folder_name))
                            continue
                        end
                    end


                    roi_ind = find(contains(s.roi_names,roi_name));
                    if numel(roi_ind) > 1
                        display(sprintf('%s contains multiple equal ROI names (%s)',stats.s(n_exp).exp_folder_name, roi_name))
                        roi_ind = roi_ind(1);
                    end


                    for n_sub = 1:numel(s.select_subfolders)
                        c = c + 1;

                        y1 = [s.roi_lb(map_ind,roi_ind,n_sub)...
                            s.roi_mean(map_ind,roi_ind,n_sub)...
                            s.roi_ub(map_ind,roi_ind,n_sub)];

                        lb(c) = y1(1);
                        ub(c) = y1(3);

                        y(1) = y1(1);
                        y(2:3) = diff(y1);


                        bar(c,y,"stacked",'FaceColor',col(n_sub,:), 'LineWidth',lw1, 'EdgeColor',[1,1,1]);

                        if n_sub == 1
                            xtlab{end+1} = extractBefore(s.exp_folder_name,'_');
                        else
                            xtlab{end+1} = '';
                        end

                    end
                end

                %title_str = sprintf('%s - %s',map_name, roi_name);

                if c > 0
                    fig(end).h.Position(3) = 1200;
                    fig(end).name = sprintf('%s_%s.png', map_name, roi_name);
                    fig(end).lb = min(lb);
                    fig(end).ub = max(ub);

                    ylim([fig(end).lb fig(end).ub] .* [0.5, 1.05])

                    xticks(0:c)
                    xticklabels(xtlab)
                    ylabel(map_name,'interpreter','none')
                    title(title_str,'interpreter','none')
                    set(gca,'LineWidth',lw2,'Yscale','Lin','Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',fs)
                end

            end

            % adjust limits and save figure
            if c > 0
                if find_global_lim
                    LIM = [min([fig.lb]), max([fig.ub])] .* [0.98 1.02];
                    for i = 1:numel(fig)
                        ylim(gca(fig(i).h), LIM);
                    end
                end

                if stat_opt.save_fig
                    for i = 1:numel(fig)
                        fig_path = fullfile(fig_folder, fig(i).name);
                        display(sprintf('saving ... %s',fig_path))
                        print(fig(i).h, fig_path, '-dpng', stat_opt.resolution);
                    end
                end
            end

        end
        close all
        clear fig
    end

    % ---- show difference per map and roi ----
    if (show_exp_dif)

        for n_map = 1:numel(map_names)
            map_name = map_names{n_map};
            map_ind = find(contains(stats.map_names,map_name));

            clear fig

            fig = [];

            for n_roi = 1:numel(roi_names)
                roi_name = roi_names{n_roi};

                title_str = sprintf('%s (mean dif.)', roi_name);

                % ----- new figure ----
                fig(end+1).h = figure;
                hold on

                c = 0;
                xtlab = {''};


                for n_exp = 1:numel(stats.s)

                    % structure s contains statistical values arranged in a 3D array : map x roi x exp-subfolders
                    s = stats.s(n_exp);

                    if isempty(include_exp_folder_name)
                        if contains(s.exp_folder_name, exclude_exp_folder_name)
                            continue
                        end
                    else
                        if ~(contains(s.exp_folder_name, include_exp_folder_name) & ~contains(s.exp_folder_name, exclude_exp_folder_name))
                            continue
                        end
                    end

                    roi_ind = find(contains(s.roi_names,roi_name));
                    if numel(roi_ind) > 1
                        display(sprintf('%s contains multiple equal ROI names (%s)',stats.s(n_exp).exp_folder_name, roi_name))
                        roi_ind = roi_ind(1);
                    end

                    n_subs = numel(s.select_subfolders);

                    if n_subs > 1
                        for n_sub = 1:n_subs-1
                            c = c + 1;

                            y(c) = s.roi_mean(map_ind,roi_ind,n_sub+1) - s.roi_mean(map_ind,roi_ind,n_sub);


                            bar(c,y(c),'FaceColor',col(n_sub,:), 'LineWidth',lw1); %, 'EdgeColor',[1,1,1]);

                            if n_sub == 1
                                xtlab{end+1} = extractBefore(s.exp_folder_name,'_');
                            else
                                xtlab{end+1} = '';
                            end

                        end
                    end
                end

                if c > 0
                    fig(end).h.Position(3) = 1200;
                    fig(end).lb = min(y);
                    fig(end).ub = max(y);
                    fig(end).name = sprintf('%s_%s_dif.png', map_name, roi_name);

                    ylim([fig(end).lb fig(end).ub])

                    xticks(0:c)
                    xticklabels(xtlab)
                    ylabel(map_name,'interpreter','none')
                    title(title_str,'interpreter','none')
                    set(gca,'LineWidth',lw2,'Yscale','Lin','Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',fs)
                end

            end

            % adjust limits and save figure
            if c > 0
                if find_global_lim
                    LIM = [min([fig.lb]), max([fig.ub])]; % .* [0.98 1.02];
                    for i = 1:numel(fig)
                        ylim(gca(fig(i).h), LIM);
                    end
                end

                if stat_opt.save_fig
                    for i = 1:numel(fig)
                        fig_path = fullfile(fig_folder, fig(i).name);
                        display(sprintf('saving ... %s',fig_path))
                        print(fig(i).h, fig_path, '-dpng', stat_opt.resolution);
                    end
                end
            end


        end

        close all
        clear fig
    end


end








% ------------------------------- FUNCTIONS ----------------------------------------
function make_histograms(root_data_path, exp_folder_name, merged_names, roi_names, select_roi_names, map_names, select_map_names, ...
    map_name_prepend, title_strings, labels, BinWidths, scales, Xmins, Xmaxs, hist_opt)

% select roi_names
if ~isempty(select_roi_names)
    roi_names = roi_names(contains(roi_names, select_roi_names));
end


for n_roi = 1:numel(roi_names)
    roi_name = roi_names{n_roi};

    % roi path
    roi_fn = fullfile(root_data_path,roi_name);
    roi = mdm_nii_read(roi_fn);
    ind = find(roi);

    for n_map = 1:numel(map_names)

        if ~contains(map_names{n_map},select_map_names)
            continue
        end
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

        if hist_opt.save
            fig_path = fullfile(root_data_path, 'fig');
            mkdir(fig_path)
            figName = [extractBefore(exp_folder_name,'_') '_' merged_name extractBefore(map_name,'.nii.gz') '_' ...
                extractBefore(roi_name,'.nii.gz') '_' hist_opt.base_fig_name '.png' ];

            fig_path = fullfile(fig_path, figName);
            display(sprintf('saving ... %s',fig_path))
            print(fh, fig_path, '-dpng', hist_opt.resolution);

        end


    end

end
end


function s = make_statistics(root_data_path, exp_folder_name, merged_names, roi_names, map_names, map_name_prepend, stat_opt)
% Structure s contains statistical values arranged in a 3D array: map × ROI × experiment subfolders.

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
            s.roi_lb(n_map, n_roi, n_exp) = roi_lb;
            s.roi_ub(n_map, n_roi, n_exp) = roi_ub;


        end
    end
end
end


function [roi_lb, roi_ub] = get_roi_CI_bounds(roi_sig, roi_mean, CI)
tmp = abs(roi_sig - roi_mean);
[~,ind_tmp] = sort(tmp);

ind_tmp = ind_tmp(1:round(length(ind_tmp) * CI/100));

%     figure(1),clf
%     hold on
%     plot(roi_sig,'.')
%     plot(ind_tmp,roi_sig(ind_tmp),'or')

roi_lb = min(roi_sig(ind_tmp));
roi_ub = max(roi_sig(ind_tmp));


end
