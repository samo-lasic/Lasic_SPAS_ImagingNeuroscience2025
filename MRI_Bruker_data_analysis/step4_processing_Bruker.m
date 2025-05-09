% Main script for analyzing Bruker data after merging
% Performs various analysis cases (1–11) on merged NIfTI files including smoothing, masking, averaging, fitting, and map generation.
% Inputs: merged NIfTI files, ROI masks; Outputs: processed NIfTI files, parameter maps

% Analysis cases:
% 1 - smoothing
% 2 - mask
% 3 - average repetitions (subsample b, volume threshold)
% 4 - signal attenuation from ROI
% 5 - add geoSPAS (geometric average waveforms)
% 6 - fit decays and normalize data (output also gamma fit results and fit decays, this step could take a very long time but it should not be necessary)
% 7 - fit DT from gamma fit ADCs
% 8 - directional average
% 9 - make muFA SPAS map (use direction average data, can be normalized or not)
% 10 - DTI (nonlinear)
% 11 - raw signal differences (SPAS, geoSPAS & STE), i.e. signal difference @ max b from SPAS1-SPAS3 and geoSPAS-STE, save also s0

% Inputs:
% Composite names of nifti files that are inputs for different analysis cases (see make_fn)
% raw file name can be appended by:
% _dn - de-noise, _mc - motion corrected, _s - smooth,  _a - average repetitions, g - geometric average, _norm - normalized, _fit - fitted, _pa - powder average
% a combination of these is possible in different order, e.g. _s_a_pa
% use_nii has fields like use_nii.s = 1, use_nii.a = 1, use_nii.pa = 1, use_nii.g = 0 ...
% zero_use() sets all fields to 0 so raw data names are used (without appending)

% Fields of "use" options determining which data is used as input for different analysis cases:
% use.dn = 1; % de-noised
% use.mc = 1; % eddy/motion corrected
% use.s = 1; % smooth
% use.a = 0; % set to 1 if averaging repetitions is needed
% use.g = 1; % geometric average
% use.norm = 1;
% use.fit = 0;
% use.pa = 0; % powder average


clear all
close all

% setup code paths
restoredefaultpath
setup_code_path()

% setup data paths
% root_data_path = ...  % Optional: override root data path here
data_opt.append_folder_name = '_out';
data_path_struct = setup_data_path(data_opt);

warning ('off','all');

nii_base_str = 'merged_';


% ----------------------
do_merged = 1; % if not go into sub-folders and do separate waveforms
select_wfm_names = {'SPAS1','SPAS2','SPAS3','tLTE','STE1','STE2'}; % only if do_merged = 0;
% select_wfm_names = select_wfm_names([1 3]);

% -----------------------

do_thresh_vol = 0; % for DTI

% mask-threshold volumes-geoSPAS-DTI for which the masked signal sum is too low compared to the non-masked signal sum (relative to total signal) - (sum_s-sum_not_s)./tot_sig
thresh_vol = 2.5; % standard deviations


use_nii = zero_use();

% using de-noised & eddy/motion corrected
if (1)
    use_nii.dn = 1; use_nii.mc = 1;
    anal_case = [4]; opt.normalize = 1; % just check ROI
end

% using de-noised & eddy/motion corrected & smooth & average rep. & geometric average & powder average
if (1)
    use_nii.dn = 1; use_nii.mc = 1; use_nii.s = 1; use_nii.a = 1; use_nii.g = 1; use_nii.pa = 1;
    % define mask
    use_mask = zero_use; use_mask.dn = 1; use_mask.mc = 1; use_mask.s = 1;
    anal_case = [9];
    % anal_case = [1 2 3 5 8 9 10 11]; % fixed or in vivo rat - full SPAS with smoothing
    % anal_case = [2 3 5 8 9 10 11]; % fixed or in vivo rat - full SPAS no smoothing
end

% --- smooth size ---
% Smoothing per slice with 2D Gaussian smoothing kernel
smooth_sz = .5; % in vivo and ex vivo rat

% Ensure use_nii and use_mask are defined; if missing, initialize to indicate raw data (no processing flags enabled)
if ~exist('use_nii')
    use_nii = zero_use; % init all processing flags to zero
end

if ~exist('use_mask')
    use_mask = zero_use; % init all processing flags to zero
end


dryRUN = 0; % testing

for n_data_path = 1:numel(data_path_struct)
    root_data_path = data_path_struct(n_data_path).root_data_path;
    exp_folder_name = data_path_struct(n_data_path).exp_folder_name;
    root_data_path = fullfile(root_data_path,exp_folder_name);
    select_subfolders = data_path_struct(n_data_path).select_subfolders;

    if isfield(data_path_struct(n_data_path),'roi_names')
        roi_names = data_path_struct(n_data_path).roi_names;
        ind = ~contains(roi_names,'.nii.gz');
        roi_names(ind) = strcat(roi_names(ind),'.nii.gz');
    end

    if do_merged
        merged_names = SPAS_numbers_to_names(select_subfolders, nii_base_str);
    else
        merged_names = {};
        for n = 1:numel(select_subfolders)
            nii_paths = dir(fullfile(root_data_path,num2str(select_subfolders(n)),'**','*.nii.gz'));
            nii_paths = nii_paths(contains({nii_paths.name},select_wfm_names)); % select waveforms
            merged_names = [merged_names fullfile(extractAfter({nii_paths.folder},root_data_path),extractBefore({nii_paths.name},'.nii.gz'))];
        end
    end

    for n_names = 1:numel(merged_names)
        merged_name = merged_names{n_names};
        clear nii_fn
        clear mask_nii_fn

        for n_case = 1:length(anal_case)
            a_case = anal_case(n_case);

            switch a_case
                case 1 % 1 - smoothing

                    opt.smooth_sz = smooth_sz;

                    display(sprintf('smoothing with sz %.1f ...', opt.smooth_sz))

                    % nifti path
                    % nii_fn = fullfile(root_data_path,[merged_name '.nii.gz']);

                    use_here = use_nii;

                    % construct filename by appending suffixes based on processing steps
                    nii_fn = make_fn(merged_name, use_here, '.nii.gz');
                    nii_fn = fullfile(root_data_path, nii_fn);

                    if dryRUN
                        nii_fn = append_nii_fn(nii_fn, 's');
                    else
                        nii_fn = smooth_slices(nii_fn, opt);
                    end

                    display(sprintf('%s',nii_fn))

                case 2 % 2 - mask
                    display('mask ...')
                    % nifti path
                    if ~exist('nii_fn')
                        use_here = use_nii;
                        %use_here.a = 1;
                        use_here.g = 0;
                        use_here.pa = 0;
                        use_here.norm = 0;
                        use_here.fit = 0;
                        nii_fn = make_fn(merged_name, use_here, '.nii.gz');
                        nii_fn = fullfile(root_data_path, nii_fn);
                    end

                    display(sprintf('%s',nii_fn))

                    opt.mask.thresh = 0.1; opt.mask.Ncontours = 10;
                    opt.density_std = 3 * 1; opt.kurtosis_pow = 0.5;
                    opt.show_mask = 0;

                    if dryRUN
                        mask_nii_fn = append_nii_fn(nii_fn, 'mask');
                    else
                        mask_nii_fn = SPAS_mask(nii_fn, opt);
                    end

                    display(sprintf('%s',mask_nii_fn))

                    % --------------------
                    % threshold volumes for which the masked signal sum is too low compared to the non-masked signal sum (relative to total signal) - (sum_s-sum_not_s)./tot_sig

                    if (~dryRUN)
                        if ~exist('thresh_vol')
                            thresh_vol = 100;
                        end
                        do_plot = 0;

                        mask_nii_fn = append_nii_fn(nii_fn, 'mask');

                        good_volumes = SPAS_get_good_volumes(nii_fn, mask_nii_fn, thresh_vol, do_plot);
                        xps = mdm_xps_load(mdm_fn_nii2xps(nii_fn));
                        xps.good_vol = good_volumes;
                        mdm_xps_save(xps, mdm_xps_fn_from_nii_fn(nii_fn));
                    end

                case 3 % 3 - average repetitions
                    display('average repetitions ...')
                    % nifti path
                    if ~exist('nii_fn')
                        use_here = use_nii;
                        use_here.a = 0;
                        use_here.g = 0;
                        use_here.pa = 0;
                        use_here.norm = 0;
                        use_here.fit = 0;
                        nii_fn = make_fn(merged_name, use_here, '.nii.gz');
                        nii_fn = fullfile(root_data_path, nii_fn);
                    end
                    display(sprintf('%s',nii_fn))

                    if dryRUN
                        nii_fn = append_nii_fn(nii_fn, 'a');
                    else
                        nii_fn = SPAS_average_repetitions(nii_fn);
                    end

                    display(sprintf('%s',nii_fn))


                case 4 % 4 - check signal merged
                    if ~exist('opt')
                        opt.normalize = 1;
                    end

                    opt.normalize_to_first_b = 0; % if 0 and  opt.normalize == 1 normalize with gamma fit
                    opt.gamma_fit = 1; % for non mono exponential decays

                    opt.data_is_normalized = 0; % if 1 then guess s0 = 1 else estimate from data
                    % opt.b_low_lim = 600 *1e6; % for initial slope estimation, if 0 then use first two decay points
                    opt.b_low_lim = 1000 *1e6; % for initial slope estimation, if 0 then use first two decay points
                    opt.weighted = 0;

                    opt.check_sig_dif = 0;

                    opt.XLim = []; % apply only if not empty
                    opt.YLim = []; % apply only if not empty
                    % opt.fig.XLim = [0 1300]; % apply only if not empty
                    % opt.YLim = [4 100]; % [0.3 1.02] % apply only if not empty
                    % opt.check_sig_dif = 0; % TDD contrast

                    % opt.YLim = [0.032 1]; % rat in vivo
                    % opt.YLim = [0.112 1]; % rat ex vivo

                    opt.change_order = 0; % optionally change the order of SPAS waveforms
                    opt.show_geoSPAS = 1; % optionally remove geoSPAS
                    opt.show_separate_wfm_sig = 0; % optionally show signals for each waveform separately

                    opt.show_labels = 0;
                    opt.showTitle = 0;
                    opt.TitleFontSize = 8;
                    opt.fig.save = 0;
                    opt.fig.resolution = '-r300';


                    % asign line styles, markers and colors
                    opt.styles = {};
                    opt.styles(end+1).wfm_name = 'SPAS1';
                    opt.styles(end).line = '-'; opt.styles(end).lw = 2; opt.styles(end).col = [1 0 0]; opt.styles(end).marker = '.'; opt.styles(end).ms = 28;

                    opt.styles(end+1).wfm_name = 'SPAS2';
                    opt.styles(end).line = '-'; opt.styles(end).lw = 4; opt.styles(end).col = [0 .8 0]; opt.styles(end).marker = '.'; opt.styles(end).ms = 36;

                    opt.styles(end+1).wfm_name = 'SPAS3';
                    opt.styles(end).line = '-'; opt.styles(end).lw = 2; opt.styles(end).col = [0 .3 1]; opt.styles(end).marker = '.'; opt.styles(end).ms = 28;

                    opt.styles(end+1).wfm_name = 'STE1';
                    opt.styles(end).line = '-'; opt.styles(end).lw = 2; opt.styles(end).col = [0 0 0]; opt.styles(end).marker = 'o'; opt.styles(end).ms = 8;

                    opt.styles(end+1).wfm_name = 'STE2';
                    opt.styles(end).line = '-'; opt.styles(end).lw = 2; opt.styles(end).col = [.3 .3 .3]; opt.styles(end).marker = 'o'; opt.styles(end).ms = 8;

                    opt.styles(end+1).wfm_name = 'tLTE';
                    opt.styles(end).line = '-'; opt.styles(end).lw = 2; opt.styles(end).col = [1 .8 0]; opt.styles(end).marker = 'o'; opt.styles(end).ms = 8;

                    opt.styles(end+1).wfm_name = 'geoSPAS';
                    opt.styles(end).line = ':'; opt.styles(end).lw = 3; opt.styles(end).col = [.5 .5 .5]; opt.styles(end).marker = '.'; opt.styles(end).ms = 0.1;

                    % additional waveforms
%                     opt.styles(end+1) = opt.styles(end);
%                     opt.styles(end).wfm_name = '...';

                    opt.select_wfm_names = select_wfm_names;


                    display('check signal merged ...')
                    % nifti path
                    if ~exist('nii_fn')
                        use_here = use_nii;
                        %use_here.a = 1;
                        use_here.g = 0;
                        use_here.pa = 0;
                        use_here.norm = 0;
                        use_here.fit = 0;
                        nii_fn = make_fn(merged_name, use_here, '.nii.gz');
                        nii_fn = fullfile(root_data_path, nii_fn);
                    end
                    display(sprintf('%s',nii_fn))


                    % roi
                    roi_fn = fullfile(root_data_path,roi_names);

                    if isfile(roi_fn)
                        if ~dryRUN
                            SPAS_signal_check_merged(n_names, nii_fn, roi_fn, opt)
                        end
                    end



                case 5 % 5 - add geoSPAS (geometric average waveforms)
                    display('add geoSPAS ...')

                    % nifti path
                    if ~exist('nii_fn')
                        use_here = use_nii;
                        %use_here.a = 1;
                        use_here.g = 0;
                        use_here.pa = 0;
                        use_here.fit = 0;
                        nii_fn = make_fn(merged_name, use_here, '.nii.gz');
                        nii_fn = fullfile(root_data_path, nii_fn);
                    end
                    display(sprintf('%s',nii_fn))

                    if dryRUN
                        nii_fn = append_nii_fn(nii_fn, 'g');
                    else
                        nii_fn = SPAS_geo_average(nii_fn);
                    end

                    display(sprintf('%s',nii_fn))

                case 6 % 6 - fit decays and normalize data (output also gamma fit results and fit decays)
                    display('fit decays and normalize data ...')

                    % nifti path
                    if ~exist('nii_fn')
                        use_here = use_nii;
                        %use_here.a = 1;
                        use_here.g = 1;
                        use_here.pa = 0;
                        use_here.norm = 0;
                        use_here.fit = 0;
                        nii_fn = make_fn(merged_name, use_here, '.nii.gz');
                        nii_fn = fullfile(root_data_path, nii_fn);
                    end
                    display(sprintf('%s',nii_fn))

                    if ~exist('mask_nii_fn')
                        use_here = use_mask;
                        %use_here.a = 1;
                        use_here.g = 0;
                        use_here.pa = 0;
                        use_here.norm = 0;
                        use_here.fit = 0;
                        mask_nii_fn = make_fn(merged_name, use_here, '_mask.nii.gz');
                        mask_nii_fn = fullfile(root_data_path, mask_nii_fn);

                    end

                    % gamma fit & normalization parameters
                    opt.SPAS_fit_decays_normalize.weighted = 0;
                    %opt.SPAS_fit_decays_normalize.NRMSE_thresh = .035;

                    % b-range limit for initial slope
                    opt.SPAS_fit_decays_normalize.b_low_lim = 0*3e9; % if 0 use first two b-values

                    opt.SPAS_fit_decays_normalize.save.gamma_fit_parameters = 1;
                    opt.SPAS_fit_decays_normalize.save.gamma_fit_signal = 1; % normalized fits (not original data)
                    opt.SPAS_fit_decays_normalize.save.s0_map = 0;
                    opt.SPAS_fit_decays_normalize.save.adc_map = 0;
                    opt.SPAS_fit_decays_normalize.save.mu2_map = 0;

                    if dryRUN
                        nii_fn = append_nii_fn(nii_fn, 'norm');
                    else
                        out = SPAS_fit_decays_normalize(nii_fn,mask_nii_fn,opt.SPAS_fit_decays_normalize);
                        nii_fn = out.norm_nii_fn;

                        if isempty(out)
                            display(sprintf('%s','missing data'))
                        else
                            f = fields(out);
                            for n = 1:numel(f)
                                display(sprintf('%s',out.(f{n})))
                            end
                        end
                    end

                case 7 % 7 - fit DT from gamma fit ADCs
                    display('DTI from gamma ...')

                    % nifti path
                    if ~exist('nii_fn')
                        use_here = use_nii;
                        %use_here.a = 1;
                        use_here.norm = 1;
                        use_here.g = 1;
                        use_here.pa = 0;
                        nii_fn = make_fn(merged_name, use_here, '.nii.gz');
                        nii_fn = fullfile(root_data_path, nii_fn);
                    end
                    display(sprintf('%s',nii_fn))

                    if contains(nii_fn,'_fit.nii.gz')
                        fit_fn = strrep(nii_fn,'_fit.nii.gz','_par.mat');
                    else
                        fit_fn = strrep(nii_fn,'.nii.gz','_par.mat');
                    end

                    if ~dryRUN
                        out = SPAS_dti_from_gamma(nii_fn,fit_fn,'geoSPAS');

                        % outputs:
                        % out.md_nii_fn
                        % out.fa_nii_fn
                        % out.vd_nii_fn
                        % out.ap_dADC_nii_fn
                        % out.ap_rdADC_nii_fn

                        if isempty(out)
                            display(sprintf('%s','missing data'))
                        else
                            f = fields(out);
                            for n = 1:numel(f)
                                display(sprintf('%s',out.(f{n})))
                            end
                        end
                    end

                case 8 % 8 - directional average
                    display('direction averaging ...')

                    % nifti path
                    if ~exist('nii_fn')
                        use_here = use_nii;
                        %use_here.a = 1;
                        use_here.g = 1;
                        use_here.pa = 0;

                        nii_fn = make_fn(merged_name, use_here, '.nii.gz');
                        nii_fn = fullfile(root_data_path, nii_fn);

                        display(sprintf('%s',nii_fn))
                        if dryRUN
                            nii_fn = append_nii_fn(nii_fn, 'pa');
                        else
                            nii_fn = SPAS_powder_average(nii_fn);
                        end

                    else

                        display(sprintf('%s',nii_fn))
                        if dryRUN
                            nii_fn = append_nii_fn(nii_fn, 'pa');
                        else
                            nii_fn = SPAS_powder_average(nii_fn);
                        end


                        nii_fn = strrep(nii_fn,'norm.nii.gz','norm_fit.nii.gz');
                        display(sprintf('%s',nii_fn))
                    end


                case 9 %9 - make muFA - SPAS maps (use direction average data, can be normalized or not)
                    display('muFA SPAS maps ...')
                    % fit muFA
                    % this can be done on non- or normalized data

                    opt.SPAS_muFA.data_is_normalized = 0; % if 1 then guess s0 = 1 else estimate from data
                    opt.SPAS_muFA.b_low_lim = 600 *1e6; % for initial slope estimation, if 0 then use first two decay points
                    opt.SPAS_muFA.weighted = 0;
                    opt.SPAS_muFA.n_inits = 5; % number of initial guess iterations (this field can be omitted for a single guess)
                    opt.SPAS_muFA.thresh = 1; % 1-all pass, .8 - only 20% best fits pass, .85 with NRMSE gives similar threshold as .7 with RSS
                    opt.SPAS_muFA.thresh_with_NRMSE = 1;

                    % clamp maps if array is not empty
                    opt.SPAS_muFA.s0_lim = []; %[] or e.g. [0 1];
                    opt.SPAS_muFA.md_lim = []; %[] or e.g. [0 1];
                    opt.SPAS_muFA.muFA_lim = [0 1]; %[] or e.g. [0 1];
                    opt.SPAS_muFA.muFA2_lim = [0 1]; %[] or e.g. [0 1];

                    % save selected
                    opt.SPAS_muFA.save.muFA = 0;
                    opt.SPAS_muFA.save.muFA2 = 1;
                    opt.SPAS_muFA.save.s0 = 0;
                    opt.SPAS_muFA.save.md = 0;
                    opt.SPAS_muFA.save.NRMSE = 0;
                    opt.SPAS_muFA.save.RSS = 0;
                    opt.SPAS_muFA.save.mask = 0;

                    % select LTE to be used in fitting (if no option use geoSPAS)
                    % opt.SPAS_muFA.LTE_wfm_name = {'SPAS1'};
                    % opt.SPAS_muFA.LTE_wfm_name = {'SPAS2'};
                    % opt.SPAS_muFA.LTE_wfm_name = {'SPAS3'};
                    % opt.SPAS_muFA.LTE_wfm_name = {'tLTE'};

                    if ~exist('mask_nii_fn')
                        use_here = use_mask;
                        use_here.g = 0;
                        use_here.pa = 0;
                        use_here.norm = 0;
                        use_here.fit = 0;
                        mask_name = make_fn(merged_name, use_here, '_mask.nii.gz');
                        mask_nii_fn = fullfile(root_data_path, mask_name);
                    end


                    % nifti path
                    if ~exist('nii_fn')

                        use_here = use_nii;
                        %use_here.a = 1;
                        %use_here.g = 1;
                        use_here.pa = 1;
                        nii_fn = make_fn(merged_name, use_here, '.nii.gz');
                        nii_fn = fullfile(root_data_path, nii_fn);

                        display(sprintf('%s',nii_fn))

                        if ~dryRUN
                            out = SPAS_muFA(nii_fn, mask_nii_fn, opt.SPAS_muFA);
                        end

                    else

                        display(sprintf('%s',nii_fn))
                        if ~dryRUN
                            out = SPAS_muFA(nii_fn, mask_nii_fn, opt.SPAS_muFA);
                        end

                        nii_fn = strrep(nii_fn,'norm_pa.nii.gz','norm_fit_pa.nii.gz');
                        display(sprintf('%s',nii_fn))



                    end

                    if ~dryRUN
                        if isempty(out)
                            display(sprintf('%s','missing data'))
                        else
                            for n_out = 1:numel(out)
                                f = fields(out(n_out));
                                for n = 1:numel(f)
                                    display(sprintf('%s',out(n_out).(f{n})))
                                end
                            end
                        end
                    end




                case 10 % 10 - DTI (nonlinear)
                    display('nonlinear DTI ...')

                    % nifti path

                    if ~exist('nii_fn')
                        use_here = use_nii;
                    else
                        use_here = use_from_fn(nii_fn); % pars filename suffixes to determine applied processing steps
                    end

                    use_here.norm = 0;
                    use_here.g = 1;
                    use_here.pa = 0;
                    use_here.fit = 0;

                    nii_fn = make_fn(merged_name, use_here, '.nii.gz');
                    nii_fn = fullfile(root_data_path, nii_fn);

                    display(sprintf('%s',nii_fn))

                    if ~exist('mask_nii_fn')
                        use_here = use_mask;
                        mask_nii_fn = make_fn(merged_name, use_here, '_mask.nii.gz');
                        mask_nii_fn = fullfile(root_data_path, mask_nii_fn);
                    end

                    display(sprintf('%s',mask_nii_fn))

                    opt.wfm_name = 'geoSPAS';
                    % opt.wfm_name = {'SPAS1','SPAS3','geoSPAS'};
                    % opt.wfm_name = {'SPAS1','SPAS3'};
                    %opt.wfm_name = 'tLTE';

                    % options for dti fit
                    opt.weightedFit = 0;
                    % b-range limit for initial slope
                    opt.b_lim = 1e9; % if 0 use first two b-values

                    opt.do_thresh_vol = do_thresh_vol; % threshold volumes for which the masked signal sum is too low compared to the non-masked signal sum (relative to total signal) - (sum_s-sum_not_s)./tot_sig
                    % can skip this field or set it to 0 to avoid filtering volumes

                    opt.save_dti_fit_parameters = 1;

                    opt.save_s0_map = 1;
                    opt.save_md_map = 1;
                    opt.save_fa_map = 1;
                    opt.save_v_map = 1; % main eigen vector


                    if ~dryRUN
                        SPAS_nls_dti(nii_fn, mask_nii_fn, opt);
                    end


                case 11 % 11 - SPAS maps, i.e. signal difference @ max b from SPAS1-SPAS3 and geoSPAS-STE, save also s0
                    display('SPAS maps, signal differences @ max b ...')
                    % make subtraction maps from max. b signals - SPAS AND/OR geoSPAS-STE
                    % data may be normalized prior to this, but in that case the s0 extraction is useless

                    opt.SPAS_sig_dif_to_nifti.normalize_case = 2; % 0-no normalization, 1-divide by smin, 2-divide by smax

                    opt.SPAS_sig_dif_to_nifti.save.TDD_smin = 1;
                    opt.SPAS_sig_dif_to_nifti.save.TDD_smax = 1;
                    opt.SPAS_sig_dif_to_nifti.save.TDD_ds = 1;
                    opt.SPAS_sig_dif_to_nifti.save.TDD_dlogs = 1;

                    opt.SPAS_sig_dif_to_nifti.save.muA_smin = 1;
                    opt.SPAS_sig_dif_to_nifti.save.muA_smax = 1;
                    opt.SPAS_sig_dif_to_nifti.save.muA_ds = 1;
                    opt.SPAS_sig_dif_to_nifti.save.muA_dlogs = 1;

                    if ~exist('mask_nii_fn')

                        use_here = use_mask;
                        use_here.s = use_nii.s;
                        use_here.a = 0;
                        use_here.g = 0;
                        use_here.pa = 0;
                        use_here.norm = 0;
                        use_here.fit = 0;
                        mask_name = make_fn(merged_name, use_here, '_mask.nii.gz');
                        mask_nii_fn = fullfile(root_data_path, mask_name);
                    end

                    if ~exist('nii_fn')
                        use_here = use_nii;
                    else
                        use_here = use_from_fn(nii_fn);
                    end
                    use_here.g = 1;
                    use_here.pa = 1;

                    nii_fn = make_fn(merged_name, use_here, '.nii.gz');
                    nii_fn = fullfile(root_data_path, nii_fn);
                    display(sprintf('%s',nii_fn))

                    if ~dryRUN
                        out = SPAS_sig_dif_to_nifti(nii_fn, mask_nii_fn);

                        if isempty(out)
                            display(sprintf('%s','missing data'))
                        else
                            f = fields(out);
                            for n = 1:numel(f)
                                display(sprintf('%s',out.(f{n})))
                            end
                        end
                    end

            end
        end
    end
end


% -----------------  FUNCTIONS ------------------


function fn = make_fn(fn,use,fn_appendix)
if use.dn
    fn = [fn '_dn'];
end
if use.mc
    fn = [fn '_mc'];
end
if use.s
    fn = [fn '_s'];
end
if use.a
    fn = [fn '_a'];
end
if use.g
    fn = [fn '_g'];
end
if use.norm
    fn = [fn '_norm'];
end
if use.fit
    fn = [fn '_fit'];
end
if use.pa
    fn = [fn '_pa'];
end
fn = [fn fn_appendix];
end

function use = use_from_fn(fn)
[~,fn] = fileparts(fn);
if contains(fn,'_dn_') use.dn = 1; else use.dn = 0; end
if contains(fn,'_mc_') use.mc = 1; else use.mc = 0; end
if contains(fn,'_s_') use.s = 1; else use.s = 0; end
if contains(fn,'_a_') use.a = 1; else use.a = 0; end
if contains(fn,'_g_') use.g = 1; else use.g = 0; end
if contains(fn,'_norm_') use.norm = 1; else use.norm = 0; end
if contains(fn,'_fit_') use.fit = 1; else use.fit = 0; end
if contains(fn,'_pa_') use.pa = 1; else use.pa = 0; end
end

function use = zero_use()
use.dn = 0;
use.mc = 0;
use.s = 0;
use.a = 0;
use.g = 0;
use.norm = 0;
use.fit = 0;
use.pa = 0;
end
