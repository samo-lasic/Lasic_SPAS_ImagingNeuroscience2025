% Main script for analyzing and plotting waveforms in multiple steps.
% Analyze and plot b-tensor waveforms, 
% Extract and plot the Spectral Principal Axis System (SPAS) 
% Extract and plot the SPAS LTE projections, 
% Extract and plot the tuned LTE projections.

% step 1: perform FFT of source b-tensor waveforms
% step 2: show g(t), q(t), spectra, cumulative power trace, b-tensor, extract SPAS and tuned LTE
% Spectral anisotropy is color-coded with RGB using 3 frequency bands defined in one of three ways: 
% (1) encoding power (spectral trace) split into bands with equal power, 
% (2) normalized diffusion spectrum D(w)/D0​ with equal changes within each band, 
% (3) arbitrary frequency bands.
% SPAS is derived as the eigensystem of the low-frequency filtered b-tensor (threshold: 1/3 total power).
% Tuned LTE is extracted for a specified restriction size.
% Step 3: Perform FFT on SPAS and tuned LTE waveforms.
% Step 4: Display SPAS, tuned waveforms, and their power spectra.

% ----------------------------------------------------------

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
% step 1: Analyze original waveform (.mat).
% Computes q(t), spectra, tensor metrics, exchange sensitivity, Maxwell index.
% Extracts tuned LTE and SPAS-LTE waveforms. Saves *_info.mat and .mat files. No plotting.

% step 2: Plot SPAS + tuning results from *_info.mat.
% Shows g(t), q(t), spectra, SPAS projections, tuning contours, and RGB maps. Saves figures.

% step 3: Analyze SPAS-LTE waveforms (from step 1).
% Performs same analysis as step 1 on tunedLTE_res, SPAS1–3, etc. Saves *_info.mat. No plotting.

% step 4: Plot SPAS-LTE results from *_info.mat.
% Shows g(t), q(t), and spectra only. No tuning or SPAS projections. Saves figures.

analyze_select = [1 1;  % step 1: analyze original waveform + extract tuned/SPAS-LTE
                  0 1;  % step 2: plot SPAS projections + tuning results
                  1 2;  % step 3: analyze SPAS-LTE waveforms 
                  0 2]; % step 4: plot SPAS-LTE results
% -------------
step = 1; % 1–4: steps 2 and 4 do plotting only
% -------------

do.analyze = analyze_select(step,1); % if 0 then load saved, if 1 then calculate power spectra
do.plotting = ~do.analyze; % plot or load data for tuning

do.save_figs = 1;
fig_resolution = '-r300';%'-r600'; '-r300'; '-r72';

% some general parameters
do.show_SQ_tensors = 0; % superquadric glyphs, else peanuts 
do.showBshape = 1;
do.show_trajectory = 1;

do.showSPAS = 1;
do.make_SPAS_grey = 1;
do.show_SPAS_axes = 1;
do.show_SPAS_scatter_axes = 0;

do.tuning = 1;
tuning_R = 2.5*1e-6; %0.7*13/2*1e-6;

do.shade_flow = 0; % mark flow encoding
do.Dw_is_restricted = 1; % for colormap and tuning
do.check_background_cross_terms = 0;
do.plot_color_weighted_exchange = 1;

do.clear_info = 0; % applies only if do.analyze = 1, clear all old analysis, i.e. ..._info files

% define colors
color_XYZ(1,:) = 1.0*[1 0 0];%XX
color_XYZ(2,:) = 0.8*[0 1 0];%YY
color_XYZ(3,:) = 1.0*[0 0 1];%ZZ

color_SPAS(1,:) = 0.2*[1 1 1];%XX
color_SPAS(2,:) = 0.5*[1 1 1];%YY
color_SPAS(3,:) = 0.8*[1 1 1];%ZZ

% b-value for showing scaled waveforms
b_used = 3806 * 1e6;

% swap colors
if (1)
    tmp = color_SPAS;
    color_SPAS = color_XYZ;
    color_XYZ = tmp;
end

%1 - STE or PTE (bdelta<1, tuning, SPAS etc.), 2 - LTE only, 3 - custom
do_select = analyze_select(step,2);

switch do_select
    case 1
        do.extract_tuned_LTE_3D = 1; %tuning_parameters
        do.plot_tuning = 1; % tunedLTE, tuning contours [single D(w)]: tuning_parameters, plot_tuning_parameters
        do.save_SPAS_g = 1; %SPAS
        do.plot_waveforms = 1; % g(t), q(t), s(w): scatter_plot_parameters.wfm
        do.plot_full_spectra = 1; % scatter_plot_parameters.fullPS
        do.plot_color_weighted_spectra = 1; % b-tensor, SPAS, tuning contours [multiple D(w)]: col_spectra_parameters 
        do.plot_spectral_PAS_projections = 1; %  spectral projections along SPAS: scatter_plot_parameters.SPAS, col_spectra_parameters.spectralPAS


    case 2
        % -------------------------- JUST DO LTE -------------------------------------------------------------
        do.extract_tuned_LTE_3D = 0; %tuning_parameters
        do.plot_tuning = 0; % tunedLTE, tuning contours [single D(w)]: tuning_parameters, plot_tuning_parameters
        do.save_SPAS_g = 0; %SPAS
        do.plot_waveforms = 1; % g(t), q(t), s(w): scatter_plot_parameters.wfm
        do.plot_full_spectra = 1; % scatter_plot_parameters.fullPS
        do.plot_color_weighted_spectra = 0; % b-tensor, SPAS, tuning contours [multiple D(w)]: col_spectra_parameters 
        do.plot_spectral_PAS_projections = 0; %  spectral projections along SPAS: scatter_plot_parameters.SPAS, col_spectra_parameters.spectralPAS

    case 3
        % ---------------------------- custom  ------------------------------------------------
        do.extract_tuned_LTE_3D = 0; %tuning_parameters
        do.plot_tuning = 1; % tunedLTE, tuning contours [single D(w)]: tuning_parameters, plot_tuning_parameters
        do.save_SPAS_g = 0; %SPAS
        do.plot_waveforms = 1; % g(t), q(t), s(w): scatter_plot_parameters.wfm
        do.plot_full_spectra = 1; % scatter_plot_parameters.fullPS
        do.plot_color_weighted_spectra = 0; % b-tensor, SPAS, tuning contours [multiple D(w)]: col_spectra_parameters 
        do.plot_spectral_PAS_projections = 0; %  spectral projections along SPAS: scatter_plot_parameters.SPAS, col_spectra_parameters.spectralPAS

        do.save_SPAS_g = 0; %SPAS
        do.plot_spectral_PAS_projections = 1; %  scatter_plot_parameters.SPAS, col_spectra_parameters.spectralPAS
        do.plot_color_weighted_spectra = 1; % col_spectra_parameters

end

if do.tuning == 0 % ensure no tuning
    do.extract_tuned_LTE_3D = 0;
    do.plot_tuning = 0;
end

if do.showSPAS == 0 % ensure no SPAS
    do.plot_spectral_PAS_projections = 0;
end

% ------------ some general plotting parameters ------------------------------------------------------------

f0 = 210;
equal_power_bins = 1; % if 0 colormap from D(w), adjust with Dw_multiplication_factor

% multiplication factor for size/correlation time
Dw_multiplication_factor =  1;

az = -130; el = 40;
lightPos1 = [1 -0.5 .5]; lightPos2 = 0*[0 1 0.5];
lightPos1 = [.5 .5 .5]; lightPos2 = 0*[0 1 0.5];

%tuning_colormap_limits = []; %if [] no limits, else limit relative to traceD
if do.Dw_is_restricted
    tuning_colormap_limits = [0.252752 1.5625]; 
else
    tuning_colormap_limits = [0.832808 1.21938];
end


scatter_plot_parameters.wfm.plotLegends = 0;
scatter_plot_parameters.wfm.show_axes = 1;
scatter_plot_parameters.wfm.show_Yticks = 1;
scatter_plot_parameters.wfm.show_labels = 0;
scatter_plot_parameters.wfm.show_cum_trace = 1;
scatter_plot_parameters.wfm.normalizeXaxis = 0; % show time from 0 to 1
scatter_plot_parameters.wfm.units_g = 1; %1 - T/m, 2 - mT/m
scatter_plot_parameters.wfm.units_q = 3; %1 - 1/m, 2 - 1/mm, 3 - 1/um

scatter_plot_parameters.wfm.spectrum_hide_Yticks = 1;

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

% ensure no extension
[~, waveform_names, ~] = fileparts(waveform_names);
waveform_names = cellstr(waveform_names);

if (do_select == 2)
    base_waveform_name = waveform_names;
    clear waveform_names
    waveform_names = {};
    for n = 1:numel(base_waveform_name)
        waveform_names{end+1} = sprintf('%s_%s',base_waveform_name{n},'optTunedLTE_res');
        waveform_names{end+1} = sprintf('%s_%s',base_waveform_name{n},'tunedLTE_res');
        waveform_names{end+1} = sprintf('%s_%s',base_waveform_name{n},'SPAS1');
        waveform_names{end+1} = sprintf('%s_%s',base_waveform_name{n},'SPAS2');
        waveform_names{end+1} = sprintf('%s_%s',base_waveform_name{n},'SPAS3');
    end
end

% find all waveforms matching waveform_names
waveform_names_tmp = {};
for n = 1:numel(waveform_names)
    waveform_names_tmp(end+1) = existing_waveform_names(find(strcmp(existing_waveform_names,waveform_names{n})));
end
waveform_names = waveform_names_tmp;


% waveform info - interpolation
% if dt > 0 &  encoding_time = 0 then we set encoding_time = samples*dt (no interpolation)
% if dt = 0 &  encoding_time > 0 then we set dt = encoding_time/samples (no interpolation)
% if dt = 0 &  encoding_time = 0 then we interpolate samples = encoding_time/dt (do interpolation)

dt = 2e-05; % Bruker
encoding_time = 0; % if 0 then set automatically
set_encoding_power = 0; % if 0 then use original from waveform else set new b


% FT info - for analysis
N0 = 1e7; % 0-padding


%--------------------  PARAMETERS --------------------------------------

% ---------------- SCATTER PLOT PARAMETERS -------------
scatter_plot_parameters = get_scatter_plot_parameters_wfm(scatter_plot_parameters);

% ---------------- SCATTER PLOT PARAMETERS: SPECTRA -------------
scatter_plot_parameters = get_scatter_plot_parameters_fullPS(scatter_plot_parameters);

%scatter_plot_parameters.fullPS.showPatch = 0;
scatter_plot_parameters.fullPS.thresh = 0.998; %if 0 then use tensor_evolution_parameters.thresh_PSf
scatter_plot_parameters.fullPS.f0 = f0; %if >0 overrides thresh, if 0 then obtain from col_spectra_parameters.thresh


% ---------------- SCATTER PLOT PARAMETERS: spectral PAS (SPAS) -------------
% search: do.plot_spectral_PAS_projections
scatter_plot_parameters.SPAS = scatter_plot_parameters.wfm;
scatter_plot_parameters.SPAS.show_axes = do.show_SPAS_scatter_axes;

scatter_plot_parameters.SPAS.frq_ticks = []; % for frequency axis
scatter_plot_parameters.SPAS.PS_labels = {'LF power','MF power','HF power'};

scatter_plot_parameters.wfm.color = color_XYZ;
scatter_plot_parameters.SPAS.color = color_SPAS;

% ---------------- PLOT TENSORS -------------
% tensor_view_parameters = get_tensor_evolution_parameters(scatter_plot_parameters);
tensor_view_parameters.showSQ = 1;
tensor_view_parameters.az = az;
tensor_view_parameters.el = el;
if ~isempty(lightPos1)
    tensor_view_parameters.lightPos1 = lightPos1;
end
if ~isempty(lightPos2)
    tensor_view_parameters.lightPos2 = lightPos2;
end


% ---------------- TUNING PARAMETERS -------------
% --- extract tunedLTE ----
tuning_parameters.save_tunedLTE_wfm = 1;

tuning_parameters.contour.restricted = do.Dw_is_restricted; % if 0 then mimic incoherent flow

tuning_parameters.contour.D0 = 2e-9;
if tuning_parameters.contour.restricted
    % restricted
    R = Dw_multiplication_factor*tuning_R;
    tuning_parameters.contour.D2Rm4 = tuning_parameters.contour.D0^2./R.^4;
else
    % incoherent flow
    tau = Dw_multiplication_factor^2*0.5*1e-2;
    tuning_parameters.contour.D2Rm4 = 1./tau.^2;
end

tuning_parameters.thresh = 0.999; % this needs to be higher than for color weighted spectra
tuning_parameters.contour.limits = 0.1*[1 1];

% ----- PLOT TUNING -------
% search: do.plot_tuning

plot_tuning_parameters = get_plot_tuning_parameters(tuning_parameters);

plot_tuning_parameters.plot.showSQ = do.show_SQ_tensors;
plot_tuning_parameters.plot.show_trajectory = do.show_trajectory;
plot_tuning_parameters.plot.showBshape = do.showBshape;
plot_tuning_parameters.plot.colormap_limits = tuning_colormap_limits;

plot_tuning_parameters.plot.az = az;
plot_tuning_parameters.plot.el = el;
if ~isempty(lightPos1)
    plot_tuning_parameters.plot.lightPos1 = lightPos1;
end
if ~isempty(lightPos2)
    plot_tuning_parameters.plot.lightPos2 = lightPos2;
end

plot_tuning_parameters.plot.show_optTuned = 1;
plot_tuning_parameters.plot.tunedLTE_color = [1 .5 0];
plot_tuning_parameters.plot.optLTE_color = [1 .8 0];
plot_tuning_parameters.plot.LTE_MarkerSize = 32;


% ---------------- PLOT COLOR CODED DIRECTIONAL SPECTRA PROJECTIONS -------------

col_spectra_parameters = get_col_spectra_parameters(tensor_view_parameters,plot_tuning_parameters);
col_spectra_parameters.ax_color = scatter_plot_parameters.wfm.color;
col_spectra_parameters.SPAS_color = scatter_plot_parameters.SPAS.color;

plot_tuning_parameters.plot.ax_color = scatter_plot_parameters.wfm.color;

col_spectra_parameters.showSQ = do.show_SQ_tensors;
col_spectra_parameters.show_trajectory = do.show_trajectory;
col_spectra_parameters.tuning.colormap_showDw = 1; % show D(w) in the legend
col_spectra_parameters.showBshape = do.showBshape;
col_spectra_parameters.showSPAS = do.show_SPAS_axes;
col_spectra_parameters.showTuning = do.plot_tuning;

col_spectra_parameters.thresh = 0.998; % if 0 then use tensor_evolution_parameters.thresh_PSf
col_spectra_parameters.f0 = f0; %if >0 overrides thresh, if 0 then obtain from col_spectra_parameters.thresh

% ----------- D(w) for color map ----------------------------
col_spectra_parameters.colormap.restricted = tuning_parameters.contour.restricted; % if 0 then mimic incoherent flow
col_spectra_parameters.colormap.D0 = tuning_parameters.contour.D0;
col_spectra_parameters.colormap.D2Rm4 = (1-equal_power_bins)*tuning_parameters.contour.D2Rm4; % if zero use equal power bins

% ----------- FLOW ----------------------------
col_spectra_parameters.flow.do_shade = do.shade_flow;
col_spectra_parameters.flow.shade.band = 0.5; % lower value for wider band
col_spectra_parameters.flow.shade.amp = 1;
col_spectra_parameters.flow.shade.p = 18;
col_spectra_parameters.flow.shade.color = 0.0*[1 1 0];



% ----------- D(w) for contours ----------------------------
col_spectra_parameters.tuning.contour.restricted = tuning_parameters.contour.restricted; % if 0 then mimic incoherent flow
col_spectra_parameters.tuning.contour.D0 = tuning_parameters.contour.D0;

if col_spectra_parameters.tuning.contour.restricted
    % restricted
    Rs = Dw_multiplication_factor*tuning_R*[1 2];
    col_spectra_parameters.tuning.contour.D2Rm4s = col_spectra_parameters.tuning.contour.D0^2./Rs.^4;
else
    % incoherent flow
    taus = Dw_multiplication_factor^2*0.5*[10 1 0.2]*1e-2; 
    col_spectra_parameters.tuning.contour.D2Rm4s = 1./taus.^2;
end

col_spectra_parameters.az = plot_tuning_parameters.plot.az;
col_spectra_parameters.el = plot_tuning_parameters.plot.el;
col_spectra_parameters.lightPos1 = plot_tuning_parameters.plot.lightPos1;
col_spectra_parameters.lightPos2 = plot_tuning_parameters.plot.lightPos2;



% -----------------------------------------------------------------------
gmr = 26.75e7;

if do.analyze

    if do.clear_info
        for n = 1:numel(d_info_files)
            delete(fullfile(d_info_files(n).folder,d_info_files(n).name))
        end
    end

    for count = 1:numel(waveform_names)
        waveform_name = waveform_names{count};

        load(fullfile(waveform_dir,[waveform_name '.mat']))

        Nwave = size(g,1);

        if (encoding_time == 0) & (dt > 0)
            TE = Nwave*dt; %duration of the waveform without 0-padding
            t = linspace(0,TE,Nwave)';
        elseif (encoding_time > 0) & (dt == 0)
            dt = encoding_time/Nwave;
            TE = encoding_time; %duration of the waveform without 0-padding
            t = linspace(0,TE,Nwave)';
        else % do interpolation
            Nwave_new = round(encoding_time/dt);
            TE = Nwave_new*dt; %duration of the waveform without 0-padding
            tmp = linspace(0,TE,Nwave)';
            t = linspace(0,TE,Nwave_new)';
            g = interp1(tmp,g,t);
            Nwave = Nwave_new;
        end

        if exist('xps') == 1
            wfm.xps = xps;
        end

        q = gmr*cumsum(g)*dt;
        [V, L] = eig(q'*q*dt);
        b = trace(L); % source b

        if set_encoding_power ~= 0
            wfm.b = set_encoding_power;
            g = g*sqrt(set_encoding_power/b);
            q = gmr*cumsum(g)*dt;
            [V, L] = eig(q'*q*dt);

        else
            wfm.b = b;
        end

        if wfm.b ~= 0
            % ------  encoding properties -----------------------
            slew_rate = diff(g)/dt;

            %Sjölund J, Szczepankiewicz F, Nilsson M, Topgaard D, Westin C-F, Knutsson H.
            %Constrained optimization of gradient waveforms for generalized diffusion encoding. J Magn Reson 2015; 261: 157–168.

            % Szczepankiewicz F, Westin CF, Nilsson M. Gradient waveform design for tensor-valued encoding in diffusion MRI. J Neurosci Methods 2020: 109007.
            
            wfm.heat = max(sum(g.^2)*dt); % also b/TE^2

            wfm.slew_rate_max = max(abs(slew_rate(:)));
            wfm.gradient_max = max(abs(g(:)));
            wfm.slew_x_gradient_max = wfm.slew_rate_max*wfm.gradient_max; % Power

            wfm.energy_max = wfm.gradient_max^2;

            wfm.heat_factor = wfm.heat/(wfm.gradient_max^2*TE); % unitless
            wfm.b_factor = wfm.b/(gmr^2*wfm.gradient_max^2*TE^3); % unitless (1/12 for constant gradient)
            wfm.b_per_T_per_maxG = wfm.b/(wfm.gradient_max*TE);
            wfm.TE = TE;


            if (1) 
                % -----  Concomitant gradients -----------------------

                % Szczepankiewicz F, Westin C-F, Nilsson M. Maxwell-compensated design of asymmetric gradient waveforms for tensor-valued diffusion encoding. Magn Reson Med 2019; 00: 1–14.
                % "optimization should constrain the Maxwell index to be less than ~1000 (mT/m)2ms."
                % "To create additional headroom, we constrain all Maxwell‐compensated waveforms to have m ≤ 100 (mT/m)2ms."
                
                % Lasič S, Szczepankiewicz F, Dall’Armellina E, Das A, Kelly C, Plein S, et al. Motion compensated b-tensor encoding for in vivo cardiac diffusion-weighted imaging. NMR Biomed 2019: 1–16.

                % add h180(t) next to g(t) in the mat file to be loaded
                if ~exist('h180')
                    h180 = -sign(t-TE/2); % dephasing direction (flips after 180)
                end
               
                % the actual gradient - not the "effective" one
                glab = g.*repmat(h180,1,3);

                gij = vec_outer_prod(glab);
                h180_mat = size(repmat(h180,1,3,3));


                Maxwell_ij = cumsum(gij.*h180_mat)*dt; % proportional to Filip's M

                M = squeeze(sum(gij.*h180_mat*dt));
                wfm.Maxwell_index =  sqrt(trace(M*M));
                wfm.norm_Maxwell_index = wfm.Maxwell_index/wfm.gradient_max^2;

            end


            if (1) 

                % -----  Exchange  -----------------------
                % just for XYZ waveforms and a mean projection on the sphere, the exchange weighting is probably a rank 4 tensor
                
                % Ning L, Nilsson M, Lasič S, Westin C-F, Rathi Y. Cumulant expansions for measuring water exchange using diffusion MRI. J Chem Phys 2018; 148
                % Szczepankiewicz F, Westin CF, Nilsson M. Gradient waveform design for tensor-valued encoding in diffusion MRI. J Neurosci Methods 2020: 109007.
                % q4(t) = integral_0^tau q^2(x)q^2(x+t) dx
                % Gamma = 2/b^2*integral_0^tau q4(t)*t*dt

                q = gmr*cumsum(g)*dt;
                
                [~, L] = eig(q'*q*dt);
                b = trace(L);

                wfm.exchange_Gamma = exchange_gamma_tensorial(q,b,dt,[1 0 0; 0 1 0; 0 0 1]);

                load UDSRTriN1000
                u = [UDSR.x UDSR.y UDSR.z];
                exchange_Gamma_u = exchange_gamma_tensorial(q,b,dt,u);
                wfm.exchange_Gamma_min = min(exchange_Gamma_u);
                wfm.exchange_Gamma_max = max(exchange_Gamma_u);
                wfm.exchange_Gamma_mean = mean(exchange_Gamma_u);

                display(sprintf('exchange: mean = %g, range = [%g-%g] ms', wfm.exchange_Gamma_mean*1e3, wfm.exchange_Gamma_min*1e3,wfm.exchange_Gamma_max*1e3))

            end

            clear h180

            [wfm.bt, wfm.bT, wfm.f_full, wfm.f, wfm.PS_full, wfm.PS, wfm.cPS, wfm.trPS, wfm.g, wfm.q] = ...
                wfm_power_spectra_ensure_b(g, b, N0, TE);

            T = tensor_shape(wfm.bT);

            wfm.b_delta = T.Delta;
            wfm.b_eta = T.Eta;
            wfm.b_FA = T.FA;
            wfm.V = T.V; % after normalization to b = 1
            wfm.l = T.l; % after normalization to b = 1
            wfm.trace = T.m; % after normalization to b = 1

            save(fullfile(waveform_dir,[waveform_name '_info.mat']),'wfm')
            clear g wfm
        end
    end
end

if do.plotting

    if (do_select == 2)
        % ----- show SPAS waveforms in one plot

        load(fullfile(waveform_dir,[waveform_names{contains(waveform_names,'SPAS1')} '_info.mat']),'wfm')
        
        gSPAS = wfm.g;
        qSPAS = wfm.q;

        ind = find(wfm.f_full>=-f0 & wfm.f_full<=f0);
        f_full = wfm.f_full(ind);
        trPS_SPAS(:,1) = wfm.PS_full(ind,1,1) + wfm.PS_full(ind,2,2) + wfm.PS_full(ind,3,3);

        load(fullfile(waveform_dir,[waveform_names{contains(waveform_names,'SPAS2')} '_info.mat']),'wfm')
        ind = find(wfm.f_full>=-f0 & wfm.f_full<=f0);

        gSPAS = gSPAS + wfm.g;
        qSPAS = qSPAS + wfm.q;

        trPS_SPAS(:,2) = wfm.PS_full(ind,1,1) + wfm.PS_full(ind,2,2) + wfm.PS_full(ind,3,3);

        load(fullfile(waveform_dir,[waveform_names{contains(waveform_names,'SPAS3')} '_info.mat']),'wfm')

        gSPAS = gSPAS + wfm.g;
        qSPAS = qSPAS + wfm.q;

        ind = find(wfm.f_full>=-f0 & wfm.f_full<=f0);
        trPS_SPAS(:,3) = wfm.PS_full(ind,1,1) + wfm.PS_full(ind,2,2) + wfm.PS_full(ind,3,3);

        % scale to gradients used
        scale_b = sqrt(b_used./sum(qSPAS(:,1).^2)/dt);

        % figure
        count_figs = count_figs+1;
        fig_h = figure(count_figs);
        fig_h.Color = 'white';

        % --------------------------- g(t)  -----------------------------

        ax1 = axes('Position',[scatter_plot_parameters.wfm.left_margin ...
            scatter_plot_parameters.wfm.bottom_margin+2*scatter_plot_parameters.wfm.height+2*scatter_plot_parameters.wfm.yspacing1+scatter_plot_parameters.wfm.yspacing2 ...
            scatter_plot_parameters.wfm.width scatter_plot_parameters.wfm.height]);
        hold on

        X = linspace(0,1,length(gSPAS))';
        if ~scatter_plot_parameters.wfm.normalizeXaxis
            X = X*dt*length(gSPAS)*1e3;
        end
        Y = scale_b * scatter_plot_parameters.wfm.scale_g*gSPAS;

        plot_scatter_axis(X, Y, [], scatter_plot_parameters.wfm.g_label, scatter_plot_parameters.SPAS.color', '-', scatter_plot_parameters.wfm)

        % --------------------------- q(t)  -----------------------------

        ax2 = axes('Position',[scatter_plot_parameters.wfm.left_margin scatter_plot_parameters.wfm.bottom_margin+scatter_plot_parameters.wfm.height+scatter_plot_parameters.wfm.yspacing1+scatter_plot_parameters.wfm.yspacing2 scatter_plot_parameters.wfm.width scatter_plot_parameters.wfm.height]);
        hold on

        Y = scale_b * scatter_plot_parameters.wfm.scale_q*qSPAS;
        plot_scatter_axis(X, Y, scatter_plot_parameters.wfm.t_label, scatter_plot_parameters.wfm.q_label, scatter_plot_parameters.SPAS.color', '-', ...
            scatter_plot_parameters.wfm)

        ax3 = axes('Position',[scatter_plot_parameters.wfm.left_margin scatter_plot_parameters.wfm.bottom_margin+0*scatter_plot_parameters.wfm.yspacing1 scatter_plot_parameters.wfm.width scatter_plot_parameters.wfm.height]);
        hold on

        % --------------------------- power spectra  -----------------------------

        plot_scatter_axis(f_full,trPS_SPAS, scatter_plot_parameters.wfm.f_label,scatter_plot_parameters.wfm.trPS_label, scatter_plot_parameters.SPAS.color', '-', ...
            scatter_plot_parameters.wfm)


        if scatter_plot_parameters.wfm.spectrum_hide_Yticks
            yticks([])
        end


        q = scale_b * qSPAS(:,1);
        exchange_Gamma(1) = exchange_gamma_projection([q 0*q 0*q ],sum(q.^2) * dt,dt,[1 0 0]);
        q = scale_b * qSPAS(:,2);
        exchange_Gamma(2) = exchange_gamma_projection([q 0*q 0*q ],sum(q.^2) * dt,dt,[1 0 0]);
        q = scale_b * qSPAS(:,3);
        exchange_Gamma(3) = exchange_gamma_projection([q 0*q 0*q ],sum(q.^2) * dt,dt,[1 0 0]);

        display(sprintf('Gamma SPAS1: %.2f ms ', exchange_Gamma(1)*1e3))
        display(sprintf('Gamma SPAS2: %.2f ms ', exchange_Gamma(2)*1e3))
        display(sprintf('Gamma SPAS3: %.2f ms ', exchange_Gamma(3)*1e3))

        g = scale_b * gSPAS(:,1);
        display(sprintf('max(g)/slew(g) SPAS1: %.3f / %.3f  ', max(abs(g)), max(abs(diff(g)/dt))))
        g = scale_b * gSPAS(:,2);
        display(sprintf('max(g)/slew(g) SPAS2: %.3f / %.3f  ', max(abs(g)), max(abs(diff(g)/dt))))
        g = scale_b * gSPAS(:,3);
        display(sprintf('max(g)/slew(g) SPAS3: %.3f / %.3f  ', max(abs(g)), max(abs(diff(g)/dt))))


        if do.save_figs
            waveform_name = extractBefore(waveform_names{contains(waveform_names,'SPAS1')},'_SPAS1');
            figName = [waveform_name '_SPAS_wfms.png' ]; 
            print(fig_h,fullfile(fig_dir, waveform_name, figName),'-dpng',fig_resolution);
        end

        % show tunedLTE and optTunedLTE waveforms in one plot
        load(fullfile(waveform_dir,[waveform_names{contains(waveform_names,'tunedLTE')} '_info.mat']),'wfm')
        g_tLTE = wfm.g(:,find(sum(wfm.g)));
        q_tLTE = wfm.q(:,find(sum(wfm.g)));

        ind = find(wfm.f_full>=-f0 & wfm.f_full<=f0);

        f_full = wfm.f_full(ind);
        trPS_tLTE = wfm.PS_full(ind,1,1) + wfm.PS_full(ind,2,2) + wfm.PS_full(ind,3,3);

        g_Tuned(:,1) = g_tLTE;
        q_Tuned(:,1) = q_tLTE;
        trPS_Tuned(:,1) = trPS_tLTE;

        col(:,1) = plot_tuning_parameters.plot.tunedLTE_color;

        % scale to gradients used
        b_used = 3806 * 1e6;
        scale_b = sqrt(b_used./sum(q_Tuned(:,1).^2)/dt);


        q = scale_b * q_Tuned(:,1);
        exchange_Gamma = exchange_gamma_projection([q 0*q 0*q ],sum(q.^2) * dt,dt,[1 0 0]);
        display(sprintf('Gamma tLTE: %.2f ms ', exchange_Gamma*1e3))

        g = scale_b * g_Tuned(:,1);
        display(sprintf('max(g)/slew(g) tLTE: %.3f / %.3f  ', max(abs(g)), max(abs(diff(g)/dt))))


        if plot_tuning_parameters.plot.show_optTuned
            load(fullfile(waveform_dir,[waveform_names{contains(waveform_names,'optTunedLTE')} '_info.mat']),'wfm')
            g_optLTE = wfm.g(:,find(sum(wfm.g)));
            q_optLTE = wfm.q(:,find(sum(wfm.g)));

            ind = find(wfm.f_full>=-f0 & wfm.f_full<=f0);

            f_full = wfm.f_full(ind);
            trPS_optLTE = wfm.PS_full(ind,1,1) + wfm.PS_full(ind,2,2) + wfm.PS_full(ind,3,3);

            g_Tuned(:,2) = g_optLTE;
            q_Tuned(:,2) = q_optLTE;
            trPS_Tuned(:,2) = trPS_optLTE;

            col(:,2) = plot_tuning_parameters.plot.optLTE_color;

            q = scale_b * q_Tuned(:,2);
            exchange_Gamma = exchange_gamma_projection([q 0*q 0*q ],sum(q.^2) * dt,dt,[1 0 0]);
            display(sprintf('Gamma optLTE: %.2f ms ', exchange_Gamma*1e3))

            g = scale_b * g_Tuned(:,2);
            display(sprintf('max(g)/slew(g) tLTE: %.3f / %.3f  ', max(abs(g)), max(abs(diff(g)/dt))))

        end


        % figure
        count_figs = count_figs+1;
        fig_h = figure(count_figs);
        fig_h.Color = 'white';

        % --------------------------- g(t)  -----------------------------

        ax1 = axes('Position',[scatter_plot_parameters.wfm.left_margin ...
            scatter_plot_parameters.wfm.bottom_margin+2*scatter_plot_parameters.wfm.height+2*scatter_plot_parameters.wfm.yspacing1+scatter_plot_parameters.wfm.yspacing2 ...
            scatter_plot_parameters.wfm.width scatter_plot_parameters.wfm.height]);
        hold on

        X = linspace(0,1,length(g_Tuned))';
        if ~scatter_plot_parameters.wfm.normalizeXaxis
            X = X*dt*length(g_Tuned)*1e3;
        end
        Y = scale_b * scatter_plot_parameters.wfm.scale_g*g_Tuned;

        plot_scatter_axis(X, Y, [], scatter_plot_parameters.wfm.g_label, col, '-', scatter_plot_parameters.wfm)

        % --------------------------- q(t)  -----------------------------

        ax2 = axes('Position',[scatter_plot_parameters.wfm.left_margin scatter_plot_parameters.wfm.bottom_margin+scatter_plot_parameters.wfm.height+scatter_plot_parameters.wfm.yspacing1+scatter_plot_parameters.wfm.yspacing2 scatter_plot_parameters.wfm.width scatter_plot_parameters.wfm.height]);
        hold on

        Y = scale_b * scatter_plot_parameters.wfm.scale_q*q_Tuned;
        plot_scatter_axis(X, Y, scatter_plot_parameters.wfm.t_label, scatter_plot_parameters.wfm.q_label, col, '-', ...
            scatter_plot_parameters.wfm)

        ax3 = axes('Position',[scatter_plot_parameters.wfm.left_margin scatter_plot_parameters.wfm.bottom_margin+0*scatter_plot_parameters.wfm.yspacing1 scatter_plot_parameters.wfm.width scatter_plot_parameters.wfm.height]);
        hold on

        % --------------------------- power spectra  -----------------------------
        plot_scatter_axis(f_full, trPS_Tuned, scatter_plot_parameters.wfm.f_label,scatter_plot_parameters.wfm.trPS_label, col, '-', ...
            scatter_plot_parameters.wfm)


        if scatter_plot_parameters.wfm.spectrum_hide_Yticks
            yticks([])
        end

        if do.save_figs
            waveform_name = extractBefore(waveform_names{contains(waveform_names,'tunedLTE')},'_tunedLTE');
            figName = [ waveform_name '_tuned_wfms.png' ]; 
            print(fig_h,fullfile(fig_dir, waveform_name, figName),'-dpng',fig_resolution);
        end


    end


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


        % --------------------------- plot just time-domain waveforms (g and q)  -----------------------------
        if do.plot_waveforms

            % scale to gradients used
            scale_b = sqrt(b_used./(trace(wfm.q' * wfm.q * dt)));

            % --------------------------- g(t)  -----------------------------
            count_figs = count_figs+1;
            fig_h = figure(count_figs);
            fig_h.Color = 'white';

            clf
            ax1 = axes('Position',[scatter_plot_parameters.wfm.left_margin ...
                scatter_plot_parameters.wfm.bottom_margin+2*scatter_plot_parameters.wfm.height+2*scatter_plot_parameters.wfm.yspacing1+scatter_plot_parameters.wfm.yspacing2 ...
                scatter_plot_parameters.wfm.width scatter_plot_parameters.wfm.height]);
            hold on

            X = linspace(0,1,length(gx))';
            if ~scatter_plot_parameters.wfm.normalizeXaxis
                X = X*dt*length(gx)*1e3;
            end
            Y = scale_b * scatter_plot_parameters.wfm.scale_g*[gx gy gz];

            nonzero_ind = sum(abs(Y))/sum(abs(Y(:))) > 0.01;% find zero elements
            Y = Y(:,nonzero_ind);
            col = scatter_plot_parameters.wfm.color(nonzero_ind,:)';

            legendText1{:,3} = scatter_plot_parameters.wfm.legendText(nonzero_ind);
            legendText1{:,2} = scatter_plot_parameters.wfm.legendText(nonzero_ind);
            legendText1{:,1} = 'PS';

            plot_scatter_axis(X, Y, [], scatter_plot_parameters.wfm.g_label, col, '-', scatter_plot_parameters.wfm)

            % --------------------------- q(t)  -----------------------------

            ax2 = axes('Position',[scatter_plot_parameters.wfm.left_margin scatter_plot_parameters.wfm.bottom_margin+scatter_plot_parameters.wfm.height+scatter_plot_parameters.wfm.yspacing1+scatter_plot_parameters.wfm.yspacing2 scatter_plot_parameters.wfm.width scatter_plot_parameters.wfm.height]);
            hold on

            Y = scale_b * scatter_plot_parameters.wfm.scale_q*[qx qy qz];

            nonzero_ind = sum(abs(Y))/sum(abs(Y(:))) > 0.01;% find zero elements
            Y = Y(:,nonzero_ind);
            col = scatter_plot_parameters.wfm.color(nonzero_ind,:)';


            plot_scatter_axis(X, Y, scatter_plot_parameters.wfm.t_label, scatter_plot_parameters.wfm.q_label, col, '-', ...
                scatter_plot_parameters.wfm)

            ax3 = axes('Position',[scatter_plot_parameters.wfm.left_margin scatter_plot_parameters.wfm.bottom_margin+0*scatter_plot_parameters.wfm.yspacing1 scatter_plot_parameters.wfm.width scatter_plot_parameters.wfm.height]);

            hold on

            % --------------------------- spectral trace (mean power spectrum)  -----------------------------

            if scatter_plot_parameters.fullPS.thresh == 0
                scatter_plot_parameters.fullPS.thresh = tensor_evolution_parameters.thresh_PSf;
            end

            if scatter_plot_parameters.fullPS.f0 == 0
                f0 = f(max(find(trPS<scatter_plot_parameters.fullPS.thresh)));
            else
                f0 = scatter_plot_parameters.fullPS.f0;
            end

            ind = find(f_full>=-f0 & f_full<=f0);

            X = f_full(ind);
            Y = PS_full(ind,1,1) + PS_full(ind,2,2) + PS_full(ind,3,3);
            col = 0*[1 1 1]';

            plot_scatter_axis(X,Y,scatter_plot_parameters.wfm.f_label,scatter_plot_parameters.wfm.trPS_label, col, '-', ...
                scatter_plot_parameters.wfm)

            if scatter_plot_parameters.wfm.show_cum_trace
                maxY = max(Y);
                Y = trPS(f<f0);
                Y = [flipud(Y); Y(2:end)];
                X = f(f<f0);
                X = [-fliplr(X) X(2:end)];
                hold on
                plot(X,Y/max(Y)*maxY,':k','LineWidth',scatter_plot_parameters.wfm.lwIn)
            end

            if scatter_plot_parameters.wfm.spectrum_hide_Yticks
                yticks([])
            end


            if do.save_figs
                figName = [waveform_full_name '_waveforms.png' ]; 
                print(fig_h,fullfile(fig_dir, waveform_full_name, figName),'-dpng',fig_resolution);
            end

            % ------   show legend separately -------
            if scatter_plot_parameters.wfm.plotLegends

                if scatter_plot_parameters.wfm.showPatch
                    legendType = 'patch';
                else
                    legendType = 'line';
                end
                h = fig_h;


                count_figs = count_figs+1;
                fig_h = figure(count_figs);
                fig_h.Color = 'white';
                clf
                h = copyobj(h.Children,fig_h);
                for n = 1:numel(h)
                    h(n).XLim(1) = 0; % to avoid clashing with legend
                    h(n).XLim(2) = h(n).XLim(2)*2; % to avoid clashing with legend

                    h_obj = findobj(h(n),'Type',legendType);
                    lh = legend(flipud(h_obj),legendText1{:,n});
                    lh.Location = 'best';
                    lh.Box = 'off';
                end

                if do.save_figs
                    figName = [waveform_full_name  '_waveform_legend.png']; 
                    print(fig_h,fullfile(fig_dir, waveform_full_name, figName),'-dpng',fig_resolution);
                end
            end

        end


        % --------------------------- full cross correlation spectra  -----------------------------
        if do.plot_full_spectra
            if scatter_plot_parameters.fullPS.thresh == 0
                scatter_plot_parameters.fullPS.thresh = tensor_evolution_parameters.fullPS.thresh_PSf;
            end

            if scatter_plot_parameters.fullPS.f0 == 0
                f0 = f(max(find(trPS<scatter_plot_parameters.fullPS.thresh)));
            else
                f0 = scatter_plot_parameters.fullPS.f0;
            end
            scatter_plot_parameters.fullPS_new = scatter_plot_parameters.fullPS;
            scatter_plot_parameters.fullPS_new.f0 = f0;


            ind = find(f_full>=-f0 & f_full<=f0);

            X = f_full(ind);

            Y = [PS_full(ind,1,1) PS_full(ind,2,2) PS_full(ind,3,3) ...
                PS_full(ind,1,2) PS_full(ind,2,3) PS_full(ind,1,3)];

            realY = real(Y);
            imagY = imag(Y);

            real_ind = sum(abs(realY))/sum(abs(realY(:))) > 0.01;% find zero elements
            imag_ind = sum(abs(imagY))/sum(abs(imagY(:))) > 0.01;% find zero elements

            Y = [realY(:,real_ind) imagY(:,imag_ind)];

            for n = 1:6
                real_LineStyle{n} = scatter_plot_parameters.fullPS.real_LineStyle;
                imag_LineStyle{n} = scatter_plot_parameters.fullPS.imag_LineStyle;
            end

            LineStyles = [real_LineStyle(real_ind) imag_LineStyle(imag_ind)];
            col = [scatter_plot_parameters.fullPS.colors(:,real_ind) scatter_plot_parameters.fullPS.colors(:,imag_ind)];
            legendText1 = {'','','','','',''};

            count_figs = count_figs+1;
            hplotFullPSf = figure(count_figs);
            hplotFullPSf.Color = 'white';

            axes('Position',[scatter_plot_parameters.fullPS.left_margin scatter_plot_parameters.fullPS.bottom_margin ...
                scatter_plot_parameters.fullPS.width scatter_plot_parameters.fullPS.height])

            hold on

            % plot full PS
            plot_scatter_axis(X,Y,scatter_plot_parameters.fullPS.f_label,scatter_plot_parameters.fullPS.PS_label, col, '-', scatter_plot_parameters.fullPS)

            axis tight

            if ~isempty(scatter_plot_parameters.fullPS.frq_ticks)
                xticks(scatter_plot_parameters.fullPS.frq_ticks)
            end

            set(gca,'LineWidth',scatter_plot_parameters.fullPS.lwOut,'Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',...
                scatter_plot_parameters.fullPS.fs)
            if ~scatter_plot_parameters.fullPS.show_axes
                axis off
            elseif ~scatter_plot_parameters.fullPS.show_Yticks
                set(gca,'ytick',[])
            end

            if scatter_plot_parameters.fullPS.showPatch
                legendType = 'patch';
            else
                legendType = 'line';
            end

            h = findobj(gca,'Type',legendType);

            if (scatter_plot_parameters.fullPS.plotLegends)
                % show legend separately
                h = gca;

                count_figs = count_figs+1;
                hplotfullPSf_legend = figure(count_figs);
                hplotfullPSf_legend.Color = 'white';
                clf
                h = copyobj(h,hplotfullPSf_legend);
                h.XLim = h.XLim*2; % avoid clashing with legend
                hplotfullPSf_legend.Visible = 'on';
                if (0)
                    h = findobj(gca); % make them invisible
                    for count = 1:numel(h)
                        if strcmp(h(count).Type,'line') || strcmp(h(count).Type,'patch')
                            h(count).Visible = 'off';
                        end
                    end
                end
                h = findobj(gca,'Type',legendType);

                lh = legend(flipud(h),legendText1);
                lh.Location = 'best';
                legend boxoff
                lh.ItemTokenSize(1) = lh.ItemTokenSize(1)*scatter_plot_parameters.fullPS.legend_width*scatter_plot_parameters.fullPS.lwOut;

            end


            if do.save_figs
                figName = [waveform_full_name  '_fullPS.png' ];
                print(hplotFullPSf,fullfile(fig_dir, waveform_full_name, figName),'-dpng',fig_resolution);

                if exist('hplotfullPSf_legend')
                    figName = [waveform_full_name  '_fullPS_legend.png' ];
                    print(hplotfullPSf_legend,fullfile(fig_dir, waveform_full_name, figName),'-dpng',fig_resolution);
                end
            end


        end

        % --------------------------- power spectra  -----------------------------

        g = [gx gy gz]; 
        q = [qx qy qz];

        % --------------------------- plot full spectra along the spectral PAS (low frequency band)  -----------------------------
        if do.plot_spectral_PAS_projections


            % --------------    find spectral PAS
            spectral_PAS = find_SPAS(f,PS,trPS, col_spectra_parameters.SPAS.LF_power_thresh);

            % select the relevant range
            if scatter_plot_parameters.fullPS.thresh == 0
                scatter_plot_parameters.fullPS.thresh = tensor_evolution_parameters.thresh_PSf;
            end

            if scatter_plot_parameters.fullPS.f0 == 0
                f0 = f(max(find(trPS<scatter_plot_parameters.fullPS.thresh)));
            else
                f0 = scatter_plot_parameters.fullPS.f0;
            end


            ind = find(abs(f_full) <= f0);
            spectral_PAS.PS_full = zeros(length(ind),3);


            for m = 1:3 % from LF to HF band
                u = spectral_PAS.filteredB_shape.V(:,m);
                outer_prod = u*u';
                % spectral projections
                spectral_PAS.PS_full(:,m) = real(sum(sum(wfm.PS_full(ind,:,:) .* permute(repmat(outer_prod,1,1,length(wfm.f_full(ind))),[3 1 2]),2),3));
            end


            % --------------------------- spectral PAS - PS projections along LF, MF and HF axes -----------------------------
            count_figs = count_figs+1;
            fig_h = figure(count_figs);
            fig_h.Color = 'white';

            X = f_full(ind);

            clf
            ax1 = axes('Position',[scatter_plot_parameters.SPAS.left_margin ...
                scatter_plot_parameters.SPAS.bottom_margin+2*scatter_plot_parameters.SPAS.height+...
                2*scatter_plot_parameters.SPAS.yspacing1+0*scatter_plot_parameters.SPAS.yspacing2 ...
                scatter_plot_parameters.SPAS.width scatter_plot_parameters.SPAS.height]);
            hold on

            Y = real(spectral_PAS.PS_full(:,1));

            % plot SPAS - 1
            plot_scatter_axis(X,Y,[],[], scatter_plot_parameters.SPAS.color(:,1), '-', scatter_plot_parameters.SPAS)

            ax2 = axes('Position',[scatter_plot_parameters.SPAS.left_margin ...
                scatter_plot_parameters.SPAS.bottom_margin+scatter_plot_parameters.SPAS.height+...
                scatter_plot_parameters.SPAS.yspacing1+0*scatter_plot_parameters.SPAS.yspacing2 ...
                scatter_plot_parameters.SPAS.width scatter_plot_parameters.SPAS.height]);
            hold on

            Y = real(spectral_PAS.PS_full(:,2));

            % plot SPAS - 2
            plot_scatter_axis(X,Y,[],scatter_plot_parameters.SPAS.PS_label, scatter_plot_parameters.SPAS.color(:,2), '-', scatter_plot_parameters.SPAS)


            ax3 = axes('Position',[scatter_plot_parameters.SPAS.left_margin ...
                scatter_plot_parameters.SPAS.bottom_margin+0*scatter_plot_parameters.SPAS.yspacing1 ...
                scatter_plot_parameters.SPAS.width scatter_plot_parameters.SPAS.height]);
            hold on

            Y = real(spectral_PAS.PS_full(:,3));

            % plot SPAS - 3
            plot_scatter_axis(X,Y,scatter_plot_parameters.SPAS.f_label,[], scatter_plot_parameters.SPAS.color(:,3), '-', scatter_plot_parameters.SPAS)


            if do.save_figs
                figName = [waveform_full_name '_spectralPAS.png' ];
                print(fig_h,fullfile(fig_dir, waveform_full_name, figName),'-dpng',fig_resolution);
            end
 
            % show legend separately 
            if scatter_plot_parameters.SPAS.plotLegends

                legendText1{:,3} = 'LF band PS';
                legendText1{:,2} = 'MF band PS';
                legendText1{:,1} = 'HF band PS';

                h = fig_h;
                count_figs = count_figs+1;
                fig_h = figure(count_figs);
                fig_h.Color = 'white';
                clf
                h = copyobj(h.Children,fig_h);

                if scatter_plot_parameters.SPAS.showPatch
                    legendType = 'patch';
                else
                    legendType = 'line';
                end

                for n = 1:numel(h)
                    h(n).XLim(1) = 0; % avoid clashing with legend
                    h(n).XLim(2) = h(n).XLim(2)*2; % avoid clashing with legend

                    h_obj = findobj(h(n),'Type',legendType);
                    lh = legend(flipud(h_obj),legendText1{:,n});
                    lh.Location = 'best';
                    lh.Box = 'off';
                end

                if do.save_figs
                    figName = [waveform_full_name  '_spectralPAS_legend.png' ];
                    print(fig_h,fullfile(fig_dir, waveform_full_name, figName),'-dpng',fig_resolution);
                end
            end

        end


        %-------------------------- show tuning ---------------------
        if do.plot_tuning

            count_figs = count_figs+1;

            [hplot, hlegend, res] = plot_tuning(count_figs, count_figs+1, f,PS, trPS, ...
                UDSR, g, q, dt, N0, plot_tuning_parameters);
            count_figs = count_figs+1;


            if plot_tuning_parameters.contour.restricted
                str = '_res';
            else
                str = '_flow';
            end


            if do.save_figs
                figName = [waveform_full_name  '_tuning_legend' str];
                figName = [figName '.png'];
                print(hlegend,fullfile(fig_dir, waveform_full_name, figName),'-dpng',fig_resolution);

                figName = [waveform_full_name  '_tuning' str];
                figName = [figName '.png'];
                print(hplot,fullfile(fig_dir, waveform_full_name, figName),'-dpng',fig_resolution);
            end

        end

        
        % --------------------------- color weighted EXCHANGE sensitivity -------------------------------
        if do.plot_color_weighted_exchange

            if col_spectra_parameters.thresh == 0
                col_spectra_parameters.thresh = tensor_evolution_parameters.thresh_PSf;
            end

            if isfield(wfm,'LF_ang')
                col_spectra_parameters.LF_ang = wfm.LF_ang;
            end

            count_figs = count_figs+1;

            bt_plot = squeeze(bt(end,:,:));
        
            [hplot, hlegend] = plot_color_weighted_exchange(count_figs, count_figs+1, f, PS, trPS, q, bt_plot, dt, wfm.b, UDSR, col_spectra_parameters); %col_exchange_parameters);
            count_figs = count_figs+1;


            if do.save_figs
                figName = [waveform_full_name  '_RGB_Exchange_colormap.png']; 
                print(hlegend,fullfile(fig_dir, waveform_full_name, figName),'-dpng',fig_resolution);
            end

            if do.save_figs
                figName = [waveform_full_name  '_RGB_Exchange.png']; 
                print(hplot,fullfile(fig_dir, waveform_full_name, figName),'-dpng',fig_resolution);
            end


        end


        % --------------------------- colored tensor (color weighted spectra) -------------------------------
        if do.plot_color_weighted_spectra
            if col_spectra_parameters.thresh == 0
                col_spectra_parameters.thresh = tensor_evolution_parameters.thresh_PSf;
            end

            if isfield(wfm,'LF_ang')
                col_spectra_parameters.LF_ang = wfm.LF_ang;
            end

            count_figs = count_figs+1;

            bt_plot = squeeze(bt(end,:,:));
            
            [hplot, hlegend, col_spectra_parameters_new] = plot_color_weighted_spectra_extended(count_figs, count_figs+1, f, PS, trPS, UDSR, g, q, bt_plot, dt, col_spectra_parameters);

            count_figs = count_figs+1;


            if col_spectra_parameters.tuning.contour.restricted == col_spectra_parameters.colormap.restricted
                if col_spectra_parameters.colormap.restricted
                    str = '_res';
                else
                    str = '_flow';
                end
            else
                str ='';
            end

            if do.save_figs
                figName = [waveform_full_name  '_RGB_Power_colormap' str];
                figName = [figName '.png'];
                print(hlegend,fullfile(fig_dir, waveform_full_name, figName),'-dpng',fig_resolution);
            
                figName = [waveform_full_name  '_RGB_Power' str];
                figName = [figName '.png'];
                print(hplot,fullfile(fig_dir, waveform_full_name, figName),'-dpng',fig_resolution);
            end


        end


        if do.save_figs

            info_path = fullfile(fig_dir, waveform_full_name,'fig_info.mat');
            if isfile(info_path) % load pre existing fields
                load(fullfile(fig_dir, waveform_full_name,'fig_info.mat'),'fig_info');
            end


            if do.plot_full_spectra
                fig_info.scatter_plot_parameters.fullPS = scatter_plot_parameters.fullPS_new;
            end

            if do.plot_color_weighted_spectra
                fig_info.col_spectra_parameters = col_spectra_parameters_new;
            end

            % save waveform parameters to print with Tex
            % remove some heavy data
            fig_info.wfm = wfm;
            fig_info.wfm = rmfield(fig_info.wfm,'bt');
            fig_info.wfm = rmfield(fig_info.wfm,'f_full');
            fig_info.wfm = rmfield(fig_info.wfm,'f');
            fig_info.wfm = rmfield(fig_info.wfm,'PS_full');
            fig_info.wfm = rmfield(fig_info.wfm,'PS');
            fig_info.wfm = rmfield(fig_info.wfm,'cPS');
            fig_info.wfm = rmfield(fig_info.wfm,'trPS');
            fig_info.wfm = rmfield(fig_info.wfm,'g');
            fig_info.wfm = rmfield(fig_info.wfm,'q');

            save(fullfile(fig_dir, waveform_full_name,'fig_info.mat'),'fig_info');

        end

        if do.extract_tuned_LTE_3D && wfm.b_delta<0.99

            count_figs = count_figs+1;

            % a temporary fix for PTE
            if wfm.b_delta == -0.5
                res = extract_tuned_LTE_2D(f,PS, trPS, UDSR, g, q, dt, N0, tuning_parameters);
            else
                res = extract_tuned_LTE_3D(f,PS, trPS, UDSR, g, q, dt, N0, tuning_parameters);
            end

            X = linspace(0,dt*length(gx),length(gx))*1e3;

            count_figs = count_figs+1;
            fh = figure(count_figs);clf, hold on, plot(X,g,'--');, plot(X,res.tunedLTE.g,'k-'), xlim([0 1.05*max(X)])
            fh.Color = 'white';
            title(sprintf('tuned LTE: tune, ||Du/trD-1|| = %g, amp = %g, slew = %g', ...
                res.tunedLTE.tuning, res.tunedLTE.maxG, res.tunedLTE.maxSlew))


            count_figs = count_figs+1;
            fh = figure(count_figs);, clf, hold on, plot(X,g,'--');, plot(X,res.optLTE.g,'k-'), xlim([0 1.05*max(X)])
            fh.Color = 'white';
            title(sprintf('optimized LTE (%s): tune, ||Du/trD-1|| = %g, amp = %g, slew = %g', ...
                res.opt, res.optLTE.tuning, res.optLTE.maxG, res.optLTE.maxSlew))


            if tuning_parameters.contour.restricted
                str = '_res';
            else
                str = '_flow';
            end


            if tuning_parameters.save_tunedLTE_wfm

                mkdir(waveform_dir)
                clear g
                g(:,3) = res.tunedLTE.g;
                save(fullfile(waveform_dir,[waveform_name '_tunedLTE' str]),'g')
            end

            if tuning_parameters.save_tunedLTE_wfm
                mkdir(waveform_dir)
                clear g
                g(:,3) = res.optLTE.g;
                save(fullfile(waveform_dir,[waveform_name '_optTunedLTE' str]),'g')
            end

        end

        if do.save_SPAS_g && wfm.b_delta<0.99
            % --------------    find spectral PAS
            spectral_PAS = find_SPAS(wfm.f,wfm.PS,wfm.trPS,col_spectra_parameters.SPAS.LF_power_thresh);

            gSPAS = wfm.g*spectral_PAS.filteredB_shape.V;

            mkdir(waveform_dir)
            for n = 1:3
                str = sprintf('_SPAS%d', n);
                g = zeros(size(gSPAS));
                g(:,n) = gSPAS(:,n);

                save(fullfile(waveform_dir,[waveform_name str]),'g')
            end


        end

        if do.check_background_cross_terms

            % find 180, assuming it is in the middle of zeros
            ind180 = find(gx == 0 & gy == 0 & gz == 0);
            tmp = find(diff(ind180) == 1);
            ind180 = round((ind180(min(tmp)) + ind180(max(tmp)))/2);
            N = length(gx);
            TE = dt*N;
            t = linspace(0,TE,N)*1e3;
            t180 = t(ind180);
            h180 = -sign(t-t180); % dephasing direction (flips after 180)
            H180 = gmr*cumsum(h180)*dt;

            q = [qx qy qz];
            c = cumsum(q.*repmat(H180',1,3));
            qmax = max(q(:));
            cmax = max(c(:));

            count_figs = count_figs+1;
            fh = figure(count_figs);
            clf
            fh.Color = 'white';
            hold on
            plot(t,qx/qmax,'-','Color',0.6*[1 0 0],'LineWidth',2)
            plot(t,qy/qmax,'-','Color',0.6*[0 0.8 0],'LineWidth',2)
            plot(t,qz/qmax,'-','Color',0.6*[0 0 1],'LineWidth',2)

            plot(t,c(:,1)/cmax,'-','Color',[1 0 0],'LineWidth',6)
            plot(t,c(:,2)/cmax,'-','Color',[0 0.8 0],'LineWidth',6)
            plot(t,c(:,3)/cmax,'-','Color',[0 0 1],'LineWidth',6)

            set(gca,'LineWidth',4,'Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',18)
            ylabel('normalize q and c (thin and thick lines)')
            xlabel('time [ms]')
        end
    end
end


