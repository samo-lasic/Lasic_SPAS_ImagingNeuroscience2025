% Read Bruker files, extract waveforms, create NIfTI files and XPS structure for analysis.
% Uses functions from the Multidimensional Diffusion MRI repository: https://github.com/markus-nilsson/md-dmri

clear all
close all

setup_code_path()

warning ('off','all'); % Suppress all warnings during conversion

data_path_struct = setup_data_path();

% ----------------------------------------------------------------------------------
% opt.root_out_path = ...  % If not defined, outputs to folders appended with "_out"

opt.outputSubFolders = 1;
opt.save_info = 1; % save basic experiment info to file
opt.rotate_image = 1; % rotate/permute phase and read pixel matrix
opt.save_waveforms = 1; % save the extracted wfm figure and g(t)
opt.show_waveform = 1; % save png
opt.include_dummy_wfm = 0; % expects waveform name dummy_..., e.g. dummy_STE1_21000_5040_20_1051_252_230302_20230307_163542_xps
opt.hide_fig = 0;
% opt.pdata = 2; set alternative pdata used for analysis; if this field is omitted pdata/1 is used

% Loop through each dataset and convert Bruker files
for c = 1:numel(data_path_struct)

    root_data_path = data_path_struct(c).root_data_path;
    exp_folder_name = data_path_struct(c).exp_folder_name;
    select_subfolders = data_path_struct(c).select_subfolders;
    display(sprintf('converting: %s', fullfile(root_data_path,exp_folder_name)))

    % Perform conversion and extraction for current dataset
    convert_FWF_Bruker(root_data_path, exp_folder_name, select_subfolders, opt)

end




