% Creates an Excel table from ROI data (saved in, e.g., stat.mat)

clear all
close all


warning ('off','all');


% ------------------------------------------

% global options

table_case = 1;
% decide if you want to export all maps and all rois and include confidence intervals
% select_map_names and select_roi_names can be empty lists (all) or contain parts of map/roi names

switch table_case
    case 1
        table_name_append = 'md_dif_fa';
        select_map_names = {'md_geoSPAS','dif_md','fa_geoSPAS'}; % {} if empty, select all {'md', 'md_geoSPAS', 'fa', 'fa_geoSPAS' , 'dif_md', 'TDD_dlogs'}
        select_roi_names = {}; %{'test1'}; % {} if empty, select all unique
        export_CI = 1;
        export_roi_size = 1;
    case 2
        table_name_append = 'md_dif';
        select_map_names = {'md_geoSPAS','dif_md'}; % {} if empty, select all {'md', 'md_geoSPAS', 'fa', 'fa_geoSPAS' , 'dif_md', 'TDD_dlogs'}
        select_roi_names = {}; %{'test1'}; % {} if empty, select all unique
        export_CI = 0;
        export_roi_size = 1;
end


% specify output path
stat_opt.out_path = '...'; 


% load stats
fn = fullfile(stat_opt.out_path,'stat.mat');
load(fn, 'stats')
display(sprintf('stats loaded from: %s', fn))

treated_names = {'xxx','yyy'}; % replace placeholder names with dataset names

include_exp_folder_name = {}; % if {} include all
exclude_exp_folder_name = {}; %{}; % if {} include all

if (0)
    include_exp_folder_name = {}; % if {} include all
    exclude_exp_folder_name = treated_names; %{}; % if {} include all
end



% --------------------


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

roi_names = roi_names(~cellfun('isempty', roi_names)); % remove empty ROI names
[~,uind] = unique(lower(roi_names));
roi_names = roi_names(uind); % keep unique ROI names 

if ~isempty(select_roi_names)
    roi_names = roi_names(contains(roi_names, select_roi_names));
end

% ensure roi_name is present in all datasets
for n_roi = 1:numel(roi_names)
    roi_ind(n_roi) = all(arrayfun(@(x) any(contains(x.roi_names,roi_names(n_roi))), stats.s));
end
roi_names = roi_names(roi_ind) % keep only ROI names present in all datasets

% select maps
map_names = stats.map_names;
if ~isempty(select_map_names)
    map_names = map_names(contains(map_names,select_map_names));
end




% ---- stats per map and roi ----

row_ind = 0;
for n_exp = 1:numel(stats.s)


    % structure s contains statistical values arranged in a 3D array : map x roi x exp-subfolders
    s = stats.s(n_exp);

    if isempty(include_exp_folder_name)
        % check exclusion criteria (applied if no inclusion filter provided)
        if contains(s.exp_folder_name, exclude_exp_folder_name) 
            continue
        end
    else
        if ~(contains(s.exp_folder_name, include_exp_folder_name) & ~contains(s.exp_folder_name, exclude_exp_folder_name))
            continue
        end
    end


    for n_sub = 1:numel(s.select_subfolders)
        row_ind = row_ind + 1;
        row_names{row_ind} = sprintf('%s [%d]',extractBefore(s.exp_folder_name,'_'), s.select_subfolders(n_sub));

        col_ind = 0;

        for n_map = 1:numel(map_names)
            map_name = map_names{n_map};
            map_ind = find(contains(stats.map_names,map_name));


            for n_roi = 1:numel(roi_names)
                roi_name = roi_names{n_roi};
                roi_ind = find(contains(s.roi_names,roi_name));

                if numel(roi_ind) > 1
                    display(sprintf('%s contains multiple equal ROI names (%s)',stats.s(n_exp).exp_folder_name, roi_name))
                    roi_ind = roi_ind(1);
                end

                col_ind = col_ind+1;

                if col_ind == 1
                    if row_ind == 1
                        col_names{col_ind} = 'is_treated';
                    end
                    table_data(row_ind,col_ind) = double(contains(s.exp_folder_name, treated_names));
                else

                    % roi - map - val
                    col_names{col_ind} = sprintf('[%s][%s][%s]', roi_name, map_name, 'mean');
                    table_data(row_ind,col_ind) = s.roi_mean(map_ind,roi_ind,n_sub);

                    if (export_CI)
                        col_ind = col_ind+1;
                        col_names{col_ind} = sprintf('[%s][%s][%s]', roi_name, map_name, 'LB');
                        table_data(row_ind,col_ind) = s.roi_lb(map_ind,roi_ind,n_sub);

                        col_ind = col_ind+1;
                        col_names{col_ind} = sprintf('[%s][%s][%s]', roi_name, map_name, 'UB');
                        table_data(row_ind,col_ind) = s.roi_ub(map_ind,roi_ind,n_sub);
                    end

                    if (export_roi_size)
                        col_ind = col_ind+1;
                        col_names{col_ind} = sprintf('[%s][%s]', roi_name, 'size');
                        table_data(row_ind,col_ind) = s.roi_size(roi_ind);
                    end

                end

            end
        end

    end

end

% remove duplicate columns (roi_size)
[~, ind] = unique(col_names,'stable');
col_names = col_names(ind);
table_data = table_data(:,ind);

T = array2table(table_data, 'RowNames', row_names, 'VariableNames', col_names)
T.Properties.DimensionNames{1} = 'experiment';

% export to Excel
% build Excel filename including table name suffix and confidence interval level
excel_fn = 'stat_table';
if numel(table_name_append) > 0
    excel_fn = [excel_fn '_'];
end

if export_CI
    excel_fn = sprintf('%s%s_CI%d.xlsx',excel_fn, table_name_append, stats.CI);
else
    excel_fn = sprintf('%s%s.xlsx',excel_fn, table_name_append);
end

excel_fn = fullfile(stat_opt.out_path, excel_fn);

writetable(T, excel_fn, 'WriteRowNames', true, 'WriteMode', 'overwrite');
