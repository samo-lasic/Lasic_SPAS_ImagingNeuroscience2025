% Generates order parameter (OP) maps from FA and µFA maps.
% Saves maps labeled with '_OP_' into experiment folders.

clear all
close all

setup_code_path() % add required folders to MATLAB path

data_opt.append_folder_name = '_out'; % output folder name to append to data path
data_path_struct = setup_data_path(data_opt);


warning ('off','all');

nii_base_str = 'merged_';


% ---------------------------------------------------
% map name prefixes reflecting different data processing stages

%opt.map_name_prepend = '_a_g_';
% opt.map_name_prepend = '_s_a_g_';
opt.map_name_prepend = '_dn_mc_a_g_';
opt.map_name_prepend = '_dn_mc_s_a_g_';

opt.map_names = {'pa_muFA', 'dti_fa_geoSPAS'};
opt.clamp_OP = 1; % clamp OP to range 0-1

% Iterate over all data paths to generate OP maps
for n_data_path = 1:numel(data_path_struct)
    root_data_path = data_path_struct(n_data_path).root_data_path;
    exp_folder_name = data_path_struct(n_data_path).exp_folder_name;
    root_data_path = fullfile(root_data_path,exp_folder_name);
    select_subfolders = data_path_struct(n_data_path).select_subfolders;

    merged_names = SPAS_numbers_to_names(select_subfolders, nii_base_str); % generate list of merged input file names

    SPAS_make_OP_maps(root_data_path, merged_names, opt); % compute and save OP maps
end






