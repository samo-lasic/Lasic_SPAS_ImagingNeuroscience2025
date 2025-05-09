function data_path_struct = setup_data_path(opt)
% Set up data paths for Bruker MRI analysis
% Returns a structure with root data path, experiment folder name, and ROI names.
% NOTE: Specify valid paths in 'default_root_data_path' and 'exp_folder_name'.


    path_case = 2; % 1 - server, 2 - local
    
    use_roi_names = 1; % if 0: set empty fields for roi_names
    
    switch path_case
        case 1
            default_root_data_path = '...';
        case 2
            default_root_data_path = '...';
    
    end
    
    % default inputs
    if nargin < 1
        data_path_struct = {};
        root_data_path = default_root_data_path;
        append_folder_name = '';
    else
    
        if isfield(opt,'data_path_struct')
            data_path_struct = opt.data_path_struct;
        else
            data_path_struct = {};
        end
    
        if isfield(opt,'root_data_path')
            root_data_path = opt.root_data_path;
        else
            root_data_path = default_root_data_path;
        end
    
        if isfield(opt,'append_folder_name')
            append_folder_name = opt.append_folder_name;
        else
            append_folder_name = '';
        end
    end
    
    % ------------------------------------------------------------------
    
    % Example entries for multiple datasets (duplicated block defines second dataset).
    % NOTE: Copy and edit this block to add more datasets as needed.

    data_path_struct = {};
    data_path_struct(end+1).root_data_path = root_data_path;
    exp_folder_name = 'xxx/zzz';  % if empty read all;
    data_path_struct(end).exp_folder_name = strcat(exp_folder_name, append_folder_name);
    data_path_struct(end).select_subfolders = [];  % if empty read all;
    data_path_struct(end).roi_names = {'roi_xxx.nii.gz', 'roi_yyy.nii.gz', ...};
    
    % data_path_struct = {};
    data_path_struct(end+1).root_data_path = root_data_path;
    exp_folder_name = 'xxx/zzz';  % if empty read all;
    data_path_struct(end).exp_folder_name = strcat(exp_folder_name, append_folder_name);
    data_path_struct(end).select_subfolders = [];  % if empty read all;
    data_path_struct(end).roi_names = {'roi_xxx.nii.gz', 'roi_yyy.nii.gz', ...};
    
    
    
    % ------------------------------------------------------------------
    
    % make roi_names empty (for now, some steps expect this field)
    if ~use_roi_names
        for n = 1:numel(data_path_struct)
            data_path_struct(n).roi_names = {};
        end
    end
    
    
    
    
    