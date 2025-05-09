
function col_spectra_parameters = get_col_spectra_parameters(tensor_view_parameters,plot_tuning_parameters);

col_spectra_parameters.showBshape = 1;
col_spectra_parameters.showSQ = 0;
col_spectra_parameters.SQoffset = 0.18;
col_spectra_parameters.showSPAS = 0;

col_spectra_parameters.thresh = 0.98; %if 0 then use tensor_evolution_parameters.thresh_PSf
col_spectra_parameters.f0 = 0; %if >0 overrides thresh, if 0 then obtain from col_spectra_parameters.thresh
col_spectra_parameters.binning = [];% 2 fractions of entire frq. range, e.g. [1/3 1/3] for even binning, if [] then even binning according to power


% bins 1
% ----------- D(w) for color map ----------------------------
col_spectra_parameters.colormap.restricted = 1; % if 0 then mimic incoherent flow
col_spectra_parameters.colormap.D0 = 2e-9;
R = 5e-6;
col_spectra_parameters.colormap.D2Rm4 = col_spectra_parameters.colormap.D0^2/R^4; % %D^2/R^4, if 0 then bin according to encoding power (lower value for encoding power at lower frequencies)


% ----------- FLOW ----------------------------

col_spectra_parameters.flow.do_shade = 1;
col_spectra_parameters.flow.shade.band = 0.5; % lower value for wider band
col_spectra_parameters.flow.shade.amp = 1;
col_spectra_parameters.flow.shade.p = 18;
col_spectra_parameters.flow.shade.color = 0.0*[1 1 0];

col_spectra_parameters.SPAS.LF_power_thresh = 1/3; % threshold for power in the LF band, i.e. LF b-tensor, for laculating spectralPAS

col_spectra_parameters.relative2mean = 1; % if 0 e.g. f^2-<f^2>, if 1 e.g. (f^2-<f^2>)/<f^2>
col_spectra_parameters.p = 2; % f^p
col_spectra_parameters.frq_interval_forSA = 0.6;  %if 0 find optimum

col_spectra_parameters.show_trajectory = 1;
col_spectra_parameters.trajectory_scale = 1.15;
col_spectra_parameters.trajectory_lw = 4;
col_spectra_parameters.ax_style = [20,0.5]; % number of segments, ratio length/spacing of segments (0-1: 0.5 is even)
col_spectra_parameters.ax_width = 2;

col_spectra_parameters.path_brightness = 1; % for showing scale
col_spectra_parameters.edge_brightness = 0.5;
col_spectra_parameters.edge_alpha = 0;
col_spectra_parameters.face_alpha = 0.88;
col_spectra_parameters.ambient = 0.95;
col_spectra_parameters.specular = 0.3;

col_spectra_parameters.az = tensor_view_parameters.az; col_spectra_parameters.el = tensor_view_parameters.el; col_spectra_parameters.lightPos1 = tensor_view_parameters.lightPos1; col_spectra_parameters.lightPos2 = tensor_view_parameters.lightPos2;

col_spectra_parameters.CameraViewAngle = 6;  % controls zoom

if nargin > 1
    % ----------  tuning - compare Dw --------------------------
    col_spectra_parameters.tuning = plot_tuning_parameters;

    col_spectra_parameters.showTuning = 0;
    col_spectra_parameters.tuning.contour.restricted = 1; % if 0 then mimic incoherent flow


    % ----------- D(w) for contours ----------------------------
    col_spectra_parameters.tuning.contour.D0 = col_spectra_parameters.colormap.D0;

    if col_spectra_parameters.tuning.contour.restricted
        % restricted
        Rs = [7 13 39]*1e-6;
        col_spectra_parameters.tuning.contour.D2Rm4s = col_spectra_parameters.tuning.contour.D0^2./Rs.^4;
    else
        % incoherent flow
        taus = [10 1 0.1]*1e-3;
        col_spectra_parameters.tuning.contour.D2Rm4s = 1./taus.^2;
    end

    col_spectra_parameters.tuning.contour.limits = 0.1*[1 1]; % plot contours, if [0 1] then full range
    col_spectra_parameters.tuning.contour.Nlevels = 1;% plot contours if 0 no contour lines
    %col_spectra_parameters.tuning.contour.style = {':','-','--','-.'};
    col_spectra_parameters.tuning.contour.style = {'-','--',':'};
end

end
