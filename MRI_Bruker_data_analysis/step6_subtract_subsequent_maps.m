% Subtracts maps from subsequent experiments defined in the `select_subfolders` list (e.g., [1 2 3] -> 2-1, 3-2).
% Useful for comparing differences between consecutive scans.
% Supports comparison of non-linear DTI fits with maps such as `norm_pa_raw_ds_TDD`.
% Outputs files like `merged_48-merged_36_s_g_pa_raw_TDD_dlogs` (subtracts `merged_36_s_g_pa_raw_TDD_dlogs` from `merged_48_s_g_pa_raw_TDD_dlogs`).


clear all
close all

setup_code_path()

data_opt.append_folder_name = '_out';
data_path_struct = setup_data_path(data_opt);


warning ('off','all');

map_names = {'dti_md_geoSPAS','dti_fa_geoSPAS','pa_raw_TDD_dlogs'};

% --------------------------
% Choose prefix for input map names (uncomment desired option)
map_name_prepend = '_s_g_';
map_name_prepend = '_g_';
% map_name_prepend = '_dn_g_';
% map_name_prepend = '_dn_mc_g_';



for n_data_path = 1:numel(data_path_struct)
    root_data_path = data_path_struct(n_data_path).root_data_path;
    exp_folder_name = data_path_struct(n_data_path).exp_folder_name;
    root_data_path = fullfile(root_data_path,exp_folder_name);
    select_subfolders = data_path_struct(n_data_path).select_subfolders;
    merged_names = SPAS_numbers_to_names(select_subfolders, 'merged_');

    for n_map = 1:numel(map_names)
        map_name = [map_name_prepend map_names{n_map} '.nii.gz'];

        for n = 2:numel(merged_names)
            merged_name1 = merged_names{n-1};
            merged_name2 = merged_names{n};

            nii_fn = fullfile(root_data_path, [merged_name1 map_name]);
            [nii1,h1] = mdm_nii_read(nii_fn);

            nii_fn = fullfile(root_data_path, [merged_name2 map_name]);
            [nii2,h2] = mdm_nii_read(nii_fn);

            name = [merged_name2 '-' merged_name1 map_name];
            nii_fn = fullfile(root_data_path,name);
            nii = nii2-nii1;

            mdm_nii_write(nii, nii_fn, h1);
            display(sprintf('%s',nii_fn));

        end

    end

end

