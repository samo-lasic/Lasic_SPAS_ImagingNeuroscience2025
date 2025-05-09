% Helper script to color-code spectral variance (Vw) or centroid frequency,
% and plot the SPAS axis (from low-pass filtered b-tensor) alongside spectral moment axes.


% can be used to compare results from different settings
if (1) % start fresh
    clear all
    close all
    count_figs = 0; % init fig counter
else
    count_figs = 1000; % init fig counter
end

restoredefaultpath
addpath(genpath(fullfile(pwd,'functions')));

% path to waveforms 
waveform_dir = fullfile('..', 'waveforms', 'M1_SPAS_2x21ms_rot');

% waveform names
waveform_names{1} = 'STE1_21000_5040_20_1051_252_220829_rot';


fig_dir = fullfile(waveform_dir,'waveform_figs')


% ------   what to do  ---------
% Select plotting options

show_LP_color = 1; % If true: color by low-pass filtered b-tensor (SPAS); if false: color by spectral variance (Vw) or centroid frequency
show_Vw_axis = 0;  % Only used if show_LP_color == 0; if true: plot spectral variance axis; if false: plot centroid frequency axis

% define colors
color_XYZ(1,:) = 1.0*[1 0 0];%XX
color_XYZ(2,:) = 0.8*[0 1 0];%YY
color_XYZ(3,:) = 1.0*[0 0 1];%ZZ

color_SPAS(1,:) = 0.2*[1 1 1];%XX
color_SPAS(2,:) = 0.5*[1 1 1];%YY
color_SPAS(3,:) = 0.8*[1 1 1];%ZZ

% swap colors (RGB -> BGR)
if (1)
    tmp = color_SPAS;
    color_SPAS = color_XYZ;
    color_XYZ = tmp;
end

color_tLTE = [1 .8 0];
color_optLTE = [1 .5 0];


do.save_figs = 0;
fig_resolution = '-r300';

% some general parameters
do.show_SQ_tensors = 0; % superquadric glyphs, else peanuts
do.showBshape = 1;
do.show_trajectory = 1;
do.show_SPAS_axes = 1;
do.plot_tuning = 0; % tunedLTE, tuning contours [single D(w)]: tuning_parameters, plot_tuning_parameters
do.plot_color_weighted_spectra = 1; % b-tensor, SPAS, tuning contours [multiple D(w)]: col_spectra_parameters


% ------------ some general plotting parameters ------------------------------------------------------------
f0 = 210;

az = -130; el = 40;
lightPos1 = [1 -0.5 .5]; lightPos2 = 0*[0 1 0.5];
lightPos1 = [.5 .5 .5]; lightPos2 = 0*[0 1 0.5];

% ------------------------------------------------------------------------------------
% check which waveforms (and analysis) exist
d = dir(waveform_dir);
d = d(find(~strcmp({d.name},'.') & ~strcmp({d.name},'..') & ~strcmp({d.name},'.DS_Store')));
existing_waveform_names = {d(find(~contains({d.name},'_info'))).name};
d_info_files = d(find(contains({d.name},'_info')));

s = '{';
for n = 1:numel(existing_waveform_names)
    existing_waveform_names{n} = extractBefore(existing_waveform_names{n},'.mat');
    s = [s sprintf('''%s'',', existing_waveform_names{n})];
end
disp('list of all waveform names:')
disp([s(1:end-1) '}'])

% -------------------------- WAVEFORMS -------------------------------------------------------------


% Ensure waveform names do not include file extensions
[~, waveform_names, ~] = fileparts(waveform_names);
waveform_names = cellstr(waveform_names);

% find all waveforms exactly matching waveform_names
waveform_names_tmp = {};
for n = 1:numel(waveform_names)
    waveform_names_tmp(end+1) = existing_waveform_names(find(strcmp(existing_waveform_names,waveform_names{n})));
end
waveform_names = waveform_names_tmp;


%-------------------- PARAMETERS --------------------------------------


% ---------------- PLOT TENSORS -------------
tensor_view_parameters.showSQ = 1;
tensor_view_parameters.az = az;
tensor_view_parameters.el = el;
if ~isempty(lightPos1)
    tensor_view_parameters.lightPos1 = lightPos1;
end
if ~isempty(lightPos2)
    tensor_view_parameters.lightPos2 = lightPos2;
end


% ---------------- PLOT COLOR CODED DIRECTIONAL SPECTRA PROJECTIONS -------------

col_spectra_parameters = get_col_spectra_parameters(tensor_view_parameters);

col_spectra_parameters.show_LP_color = show_LP_color; % else color code Vw or centroid frequency
col_spectra_parameters.show_Vw_axis = show_Vw_axis; % else show centroid frequency axis besides LP SPAS axis

col_spectra_parameters.ax_color = color_XYZ;
col_spectra_parameters.SPAS_color = color_SPAS;

col_spectra_parameters.tuning.colormap_showDw = 0; % show D(w) in the legend

col_spectra_parameters.showSQ = do.show_SQ_tensors;
col_spectra_parameters.show_trajectory = do.show_trajectory;

col_spectra_parameters.showBshape = do.showBshape;
col_spectra_parameters.showSPAS = do.show_SPAS_axes;
col_spectra_parameters.showTuning = do.plot_tuning;

col_spectra_parameters.thresh = 0.998; % if 0 then use tensor_evolution_parameters.thresh_PSf
col_spectra_parameters.f0 = f0; %if >0 overrides thresh, if 0 then obtain from col_spectra_parameters.thresh


col_spectra_parameters.az = az;
col_spectra_parameters.el = el;
col_spectra_parameters.lightPos1 = lightPos1;
col_spectra_parameters.lightPos2 = lightPos2;



% -----------------------------------------------------------------------

for count = 1:numel(waveform_names)

    waveform_name = waveform_names{count}

    if ~exist('waveform_full_names')
        waveform_full_names = waveform_names;
    end
    if length(waveform_full_names) ~= length(waveform_names)
        waveform_full_names = waveform_names;
    end

    waveform_full_name = waveform_full_names{count};

    load(fullfile(waveform_dir,[waveform_name '_info.mat']),'wfm')

    gx = wfm.g(:,1);
    gy = wfm.g(:,2);
    gz = wfm.g(:,3);

    qx = wfm.q(:,1);
    qy = wfm.q(:,2);
    qz = wfm.q(:,3);

    bt = wfm.bt;

    PS_full = wfm.PS_full;
    PS = wfm.PS;
    cPS = wfm.cPS;
    trPS = wfm.trPS;


    f_full = wfm.f_full;
    f = wfm.f;
    dt = wfm.TE/length(wfm.g);

    load UDSRTriN1000
    u = [UDSR.x UDSR.y UDSR.z];

    if do.save_figs
        mkdir(fullfile(fig_dir, waveform_full_name))
    end


    g = [gx gy gz];
    q = [qx qy qz];



% --------------------------- color weighted spectra -------------------------------
    if do.plot_color_weighted_spectra
        if col_spectra_parameters.thresh == 0
            col_spectra_parameters.thresh = tensor_evolution_parameters.thresh_PSf;
        end

        if isfield(wfm,'LF_ang')
            col_spectra_parameters.LF_ang = wfm.LF_ang;
        end

        count_figs = count_figs+1;

        bt_plot = squeeze(bt(end,:,:));

        [hplot, hlegend, ~] = plot_color_weighted_SPAS_moments(count_figs, count_figs+1, f, PS, trPS, UDSR, g, q, bt_plot, dt, col_spectra_parameters);
        hplot.Color = 'White';
        hlegend.Color = 'White';

        count_figs = count_figs+1;

        if isfield(col_spectra_parameters.tuning,'contour')
            if col_spectra_parameters.tuning.contour.restricted == col_spectra_parameters.colormap.restricted
                if col_spectra_parameters.colormap.restricted
                    str = '_res';
                else
                    str = '_flow';
                end
            else
                str ='';
            end
        end

        if do.save_figs
            figName = [waveform_full_name  '_RGB_Power_colormap' str ]; % _' num2str(col_spectra_parameters.az,3) num2str(col_spectra_parameters.el,3)];
            figName = [figName '.png'];
            print(hlegend,fullfile(fig_dir, waveform_full_name, figName),'-dpng',fig_resolution);
        end

        if do.save_figs
            figName = [waveform_full_name  '_RGB_Power' str]; % _' num2str(col_spectra_parameters.az,3) num2str(col_spectra_parameters.el,3)];
            figName = [figName '.png'];
            print(hplot,fullfile(fig_dir, waveform_full_name, figName),'-dpng',fig_resolution);
        end


    end


end



















