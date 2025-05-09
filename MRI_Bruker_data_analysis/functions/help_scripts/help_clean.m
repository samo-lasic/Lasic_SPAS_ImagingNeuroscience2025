clear all

setup_code_path()

addpath(genpath(fileparts(pwd)));

warning('off','all')

nii_base_str = 'merged_';

% !!! delete files with names containing the strings !!!

contain_strings = {'substring_in_file_name'};
% contain_strings = {'norm_pa_raw_ds_TDD'};


do_delete = 0; % if 0 dry run



root_data_path = fullfile(fileparts(fileparts(fileparts(pwd))),'...');

data_path_struct = {};
data_path_struct(end+1).root_data_path = root_data_path;
data_path_struct(end).exp_folder_name = '...'; %
data_path_struct(end).select_subfolders = []; %  % if empty read all;

%data_path_struct = {};
data_path_struct(end+1).root_data_path = root_data_path;
data_path_struct(end).exp_folder_name = '...'; %
data_path_struct(end).select_subfolders = []; %  % if empty read all;



cnt_all = 0;

for c = 1:numel(data_path_struct)
    p = data_path_struct(c);
    files = dir([fullfile(p.root_data_path,p.exp_folder_name) '/**']);
    files = files(~[files.isdir]);

    cnt = 0;
    for n = 1:length(files)

        if contains(files(n).name, contain_strings)
            cnt = cnt + 1;
            path = fullfile(files(n).folder, files(n).name);
            display(sprintf('%s', path))

            if do_delete
                delete(path);
            end
        end

    end
    display(sprintf('----- %d files deleted in %s -----', cnt, p.exp_folder_name ))

    cnt_all = cnt_all + cnt;


end

display(sprintf('----- %d files deleted in %d folders -----', cnt_all, numel(data_path_struct)))
