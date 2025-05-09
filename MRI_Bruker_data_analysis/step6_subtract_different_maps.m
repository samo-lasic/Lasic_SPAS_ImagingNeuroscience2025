% Subtract different maps like muFA from STE and geoSPAS LTE and muFA from STE and SPAS3 LTE.


clear all
close all

setup_code_path()

data_opt.append_folder_name = '_out';
data_path_struct = setup_data_path(data_opt);


warning ('off','all');

nii_base_str = 'merged_';



% ---------------------------------------------------
% example: 
% merged_xxx_a_g_pa_muFA.nii.gz
% merged_xxx_a_g_dti_fa_geoSPAS.nii.gz

% map name prefixes reflecting different data processing stages

%opt.map_name_prepend = '_a_g_';
% opt.map_name_prepend = '_s_a_g_';
% opt.map_name_prepend = '_mc_s_a_g_';
% opt.map_name_prepend = '_mc_a_g_';
%opt.map_name_prepend = '_dn_a_g_';
opt.map_name_prepend = '_dn_mc_a_g_';
opt.map_name_prepend = '_dn_mc_s_a_g_pa_';


% Define strings for input and output: '1','2','3' represent SPAS indices, 'g' represents geoSPAS.

% map_name = 'muFA';
map_name = 'muFA2';

c = 5;

if c == 1
    wfm_str = {'1','g'};
elseif c == 2
    wfm_str = {'g','3'};
elseif c == 3
    wfm_str = {'1','3'};
elseif c == 4
    wfm_str = {'1','2'};
elseif c == 5
    wfm_str = {'2','3'};
end

opt.out_name = '';
for n = 1:length(wfm_str)
    tmp_str = wfm_str{n};
    opt.out_name = strcat(opt.out_name,tmp_str);
    if strcmp(tmp_str,'g')
        opt.map_names{n} = map_name;
    else
        opt.map_names{n} = strcat('SPAS', tmp_str, '_', map_name);
    end
end
opt.out_name = strcat('SPAS_dif','_',opt.out_name,'_',map_name);


opt.lim_range = [0 1]; % optional limits (omit this field for no limits)

for n_data_path = 1:numel(data_path_struct)
    root_data_path = data_path_struct(n_data_path).root_data_path;
    exp_folder_name = data_path_struct(n_data_path).exp_folder_name;
    root_data_path = fullfile(root_data_path,exp_folder_name);
    select_subfolders = data_path_struct(n_data_path).select_subfolders;

    merged_names = SPAS_numbers_to_names(select_subfolders, nii_base_str);

    SPAS_make_difference_maps(root_data_path, merged_names, opt);
end

% ------------------------- FUNCTIONS -------------------------------------

function out_nii_fn = SPAS_make_difference_maps(root_data_path, merged_names, opt)
% Compute difference maps from two input maps (e.g., muFA maps from different sources)

if ~isfield(opt,'lim_range')
    opt.lim_range = Inf * [-1 1];
end


for n_data = 1:numel(merged_names)
    merged_name = merged_names{n_data};

    % check that both maps exist
    names = strcat(merged_name, opt.map_name_prepend, opt.map_names, '.nii.gz');
    nii_fn = fullfile(root_data_path,names);
    ind = isfile(nii_fn);

    N = length(ind);

    if N ~= 2
        display('expecting 2 files')
        continue
    end

    for n = 1:N
        if ind(n) == 0
            display(sprintf('missing: %s', nii_fn{n}))
        end
    end

    if prod(ind) == 0
        continue
    end


    for n = 1:N
        display(sprintf('%s',nii_fn{n}))
        [map(n).I, map(n).h] = mdm_nii_read(nii_fn{n});
        map(n).I = double(map(n).I);
    end

    dif_map = map(1).I - map(2).I;

    display(sprintf('dif. map has %.1f %% nans', 100*numel(find(isnan(dif_map)))/numel(dif_map)))
    dif_map(isnan(dif_map)) = 0;

    lim = opt.lim_range(1);
    display(sprintf('dif. map has %.1f %% values < %g', 100*numel(find(dif_map < lim))/numel(dif_map), lim))
    dif_map(dif_map < lim) = lim;

    lim = opt.lim_range(2);
    display(sprintf('dif. map has %.1f %% values > %g', 100*numel(find(dif_map > lim))/numel(dif_map), lim))
    dif_map(dif_map > lim) = lim;

    out_nii_fn = strcat(merged_name, opt.map_name_prepend, opt.out_name, '.nii.gz');
    out_nii_fn = fullfile(root_data_path,out_nii_fn);
    display(sprintf('saving: %s',out_nii_fn))

    mdm_nii_write(single(dif_map), out_nii_fn, map(1).h);


end

end



