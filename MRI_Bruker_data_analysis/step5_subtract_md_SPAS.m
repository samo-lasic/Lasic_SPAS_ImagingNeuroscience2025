% Subtracts MD maps from SPAS1 and SPAS3.
% Outputs maps labeled 'dti_dif_md'.
% Requires input maps '_dti_md_SPAS1' and '_dti_md_SPAS3'.
% Should be used before running step6_ROI_stats.m if 'dti_dif_md' maps are to be included in ROI statistics or visualization

clear all
close all

setup_code_path()

data_opt.append_folder_name = '_out';
data_path_struct = setup_data_path(data_opt);

warning ('off','all');


% --------------------------
% map name prefixes reflecting different data processing stages

% map_name_prepend = '_s_g_';
% map_name_prepend = '_g_';
% map_name_prepend = '_dn_g_';
map_name_prepend = '_dn_mc_g_';


for n_data_path = 1:numel(data_path_struct)
    root_data_path = data_path_struct(n_data_path).root_data_path;
    exp_folder_name = data_path_struct(n_data_path).exp_folder_name;
    root_data_path = fullfile(root_data_path,exp_folder_name);
    select_subfolders = data_path_struct(n_data_path).select_subfolders;
    merged_names = SPAS_numbers_to_names(select_subfolders, 'merged_');

    for n = 1:numel(merged_names)
        map_name = [merged_names{n} map_name_prepend 'dti_md_SPAS1.nii.gz'];
        nii_fn = fullfile(root_data_path, map_name);
        [nii1,h1] = mdm_nii_read(nii_fn);

        map_name = [merged_names{n} map_name_prepend 'dti_md_SPAS3.nii.gz'];
        nii_fn = fullfile(root_data_path, map_name);
        [nii3,h3] = mdm_nii_read(nii_fn);


        name_out = [merged_names{n} map_name_prepend 'dti_dif_md.nii.gz'];
        nii_fn = fullfile(root_data_path,name_out);
        nii = nii3-nii1;

        mdm_nii_write(nii, nii_fn, h1);
        display(sprintf('%s',nii_fn));

    end

end

