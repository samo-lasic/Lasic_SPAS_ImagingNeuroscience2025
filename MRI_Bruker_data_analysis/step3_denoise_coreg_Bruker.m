% De-noise and co-register Bruker MRI datasets.
% Applies Marchenko-Pastur denoising and extrapolation-based motion/eddy-current correction.

clear all
close all

setup_code_path()

%root_data_path = ...
data_opt.append_folder_name = '_out';
data_path_struct = setup_data_path(data_opt);


% --------- Set processing options ----------------------
nii_base_str = 'merged_';

do.denoise = 1;
do.coreg = 1;
do.coreg_single = 1; % if 1, co-register each dataset separately, else co-register multiple datasets together
use_denoised_for_coreg = 1;%[0 1];
do.clean_up = 1;


% -------------------------------

im_cor_array = [1];  % Select co-registration parameter setting (can run multiple cases)

for n_cor = 1:length(im_cor_array)
    cor_case = im_cor_array(n_cor);

    switch cor_case

        case 1

            % de-noise
            window_size = 12;

            % eddy/motion correction
            p = elastix_p_affine(200); % iterations
            p.Scales = power(10, [6 6 5  5 5 5  5 5 6  0.5 0.5 5 ]);
            % cost for rotations, sheer, scale, translation (high value not allowed)

        case 2

            % de-noise
            window_size = 6;

            % eddy/motion correction
            p = elastix_p_affine(200); % iterations
            p.Scales = power(10, [6 6 5  5 5 5  5 5 6  0.5 0.5 5 ]);
            % cost for rotations, sheer, scale, translation (high value not allowed)

        case 3

            % de-noise
            window_size = 3; 
            % eddy/motion correction
            p = elastix_p_affine(200); % iterations
            p.Scales = power(10, [6 6 5  5 5 5  5 5 6  0.5 0.5 5 ]);
            % cost for rotations, sheer, scale, translation (high value not allowed)

    end

    % common co-registration options
    opt = mdm_opt();
    opt.mdm.mec.do_rotate_bvec = 0;
    opt.mdm.mec.do_cleanup = 1; % delete extrapolation reference

    opt = mio_opt(opt);
    opt.mio.coreg.assume_las = 0; % if this is 1, error: MC: Rotate gradients only for LAS data (fix deviated) in mdm_coreg / mdm_mec_b0

    % Loop through each dataset for denoising and co-registration
    for n_data_path = 1:numel(data_path_struct)
        root_data_path = data_path_struct(n_data_path).root_data_path;
        exp_folder_name = data_path_struct(n_data_path).exp_folder_name;
        root_data_path = fullfile(root_data_path,exp_folder_name);
        select_subfolders = data_path_struct(n_data_path).select_subfolders;

        % de-noise
        if do.denoise
            for n_fn = 1:numel(select_subfolders)
                denoise1(root_data_path, select_subfolders(n_fn),nii_base_str, window_size)
            end
        end

        % eddy/motion correction
        if do.coreg

            for n_use_dn = 1:length(use_denoised_for_coreg)
                coreg1(root_data_path, nii_base_str, select_subfolders, p, opt, do, use_denoised_for_coreg(n_use_dn))
            end

        end

    end

end

function denoise1(root_data_path, select_subfolders, nii_base_str, window_size)
merged_name = SPAS_numbers_to_names(select_subfolders, nii_base_str);
ind = ~contains(merged_name,'.nii.gz');
merged_name = strcat(merged_name(ind),'.nii.gz');
merged_name = merged_name{1};

% Load NIfTI and XPS for denoising
s.nii_fn = fullfile(root_data_path, merged_name);
s.xps = mdm_xps_load(mdm_fn_nii2xps(s.nii_fn));

window = window_size * [1 1];
[nii, h] = mdm_nii_read(s.nii_fn);
[nii_dn,S2,P] = denoise(nii,window); %,mask);

% clamp negative values
nii_dn(nii_dn<0) = 0;


s.nii_fn = fullfile(root_data_path,[extractBefore(merged_name,'.nii.gz') '_dn.nii.gz']);

mdm_nii_write(nii_dn,s.nii_fn,h);

mdm_xps_save(s.xps,mdm_fn_nii2xps(s.nii_fn))
end

function coreg1(root_data_path, nii_base_str, select_subfolders, p, opt, do, use_dn)

str_array = SPAS_numbers_to_names(select_subfolders, nii_base_str);
ind = ~contains(str_array,'.nii.gz');
str_array(ind) = strcat(str_array(ind),'.nii.gz');


if use_dn
    str_append = '_dn.nii.gz';
else
    str_append = '.nii.gz';
end

nii_fn_array_in = fullfile(root_data_path,strrep(str_array,'.nii.gz',str_append));

if ~do.coreg_single
    % make temporary merged nifti for co-registration
    out_nii_fn = fullfile(root_data_path,'cor',['cor' str_append]);
    out_nii_fn = mdm_nii_merge(nii_fn_array_in, out_nii_fn);

    % collect xps info
    max_s_ind = 0;
    n_vol_split = [];
    xps_without_cells_array = {};
    xps_with_cells_array = {};

    for n_fn = 1:numel(nii_fn_array_in)
        xps_tmp = mdm_xps_load(mdm_fn_nii2xps(nii_fn_array_in{n_fn}));
        xps_tmp.s_ind = xps_tmp.s_ind + max_s_ind;
        max_s_ind = max(max_s_ind + xps_tmp.s_ind);
        n_vol_split(n_fn) = xps_tmp.n;

        [xps_without_cells_array{n_fn}, xps_with_cells_array{n_fn}] = struct_split(xps_tmp);
    end
    xps = mdm_xps_merge(xps_without_cells_array);

    % save temporary xps
    mdm_xps_save(xps, mdm_fn_nii2xps(out_nii_fn));

    nii_fn_array = {out_nii_fn};
else
    nii_fn_array = nii_fn_array_in;
    n_vol_split = [];
end

cor_o_path = fullfile(root_data_path,'cor');

for n_fn = 1:numel(nii_fn_array)

    s.nii_fn = nii_fn_array{n_fn};
    s.xps = mdm_xps_load(mdm_fn_nii2xps(s.nii_fn));
    xps_src = s.xps;

    % Write the elastix parameter file
    if do.coreg_single
        p_fn = fullfile(root_data_path,sprintf('elastix_par_%d.txt', select_subfolders(n_fn)));
        p_fn = elastix_p_write(p, p_fn);

        % reduce xps for coreg
        [s.xps, ~] = struct_split(s.xps);

    else

        p_fn = fullfile(root_data_path,'elastix_par.txt');
        p_fn = elastix_p_write(p, p_fn);
    end


    % Run an extrapolation-based registration
    s_registered = mdm_s_mec(s, p_fn, cor_o_path, opt);

    % clamp negative values
    [nii_reg, h_reg] = mdm_nii_read(s_registered.nii_fn);
    append_name = '';
    if isempty(append_name)
        nii_reg_fn = s_registered.nii_fn;
    else
        nii_reg_fn = strrep(s_registered.nii_fn,'.nii.gz', ['_' append_name '.nii.gz']);
    end
    nii_reg(nii_reg<0) = 0;

    mdm_nii_write(nii_reg, nii_reg_fn, h_reg);
    mdm_xps_save(xps_src, mdm_fn_nii2xps(nii_reg_fn));


    % If co-registration was performed on merged dataset, split merged result back into individual files.
    % Otherwise, move individual registered files to output location.

    if ~do.coreg_single
        % split temporary merged nifti and save to root

        ind_split = [0 cumsum(n_vol_split)];
        for n = 1:length(nii_fn_array_in)
            I1 = nii_reg(:,:,:,1+ind_split(n):ind_split(n+1));
            nii_reg_fn_new = strrep(nii_fn_array_in{n},'.nii.gz','_mc.nii.gz');
            nii_reg_fn_new = mdm_nii_write(I1, nii_reg_fn_new, h_reg);

            % copy xps
            copyfile(mdm_fn_nii2xps(nii_fn_array_in{n}), mdm_fn_nii2xps(nii_reg_fn_new))
        end

    else
        % move files
        [~,fn] = fileparts(extractBefore(nii_reg_fn,'.nii.gz'));
        nii_reg_fn_new = fullfile(fileparts(cor_o_path),strcat(fn,'.nii.gz'));
        copyfile(nii_reg_fn, nii_reg_fn_new)
        copyfile(mdm_fn_nii2xps(nii_reg_fn), mdm_fn_nii2xps(nii_reg_fn_new))

    end

    if do.clean_up
        rmdir(cor_o_path,'s')
    end
end
end