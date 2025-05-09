% Generates figures showing:
% - tuning landscapes across different substrate sizes, visualizing tuning proximity
%   (1 - |1 - D_u/MD|) as a function of polar and azimuthal angles
% - angular deviations between the tuned projections for R = 2.5 µm and other sizes,
%   assessing tuning robustness across restriction sizes
% Utilizes fminsearch for optimization
% Saves output to the specified 'fig' folder.

clear all
close all

restoredefaultpath
addpath(genpath(fullfile(pwd,'functions')));

% path to waveform files
waveform_dir = fullfile('..', 'waveforms', 'M1_SPAS_2x21ms_rot');

% waveform names
waveform_names{1} = 'STE1_21000_5040_20_1051_252_220829_rot';

fig_dir = fullfile(waveform_dir,'waveform_figs')

% figure options 
save_fig = 0;
lw = 4;
fs = 14;


% Generate tuning contours for different sizes
Rs_grid = [2 3 7 20 40] * 1e-6;   % sphere radii for tuning contour plots 

tuning_plot_lim = [.9, 1];


R = 2.5 * 1e-6;                   % reference substrate radius for angular comparison
Rs = [.05:.05:20] * 1e-6;         % sphere radii for optimizing tuned projection (via fminsearch) used in angle vs size plot 

PS_thresh = .999;
D0 = 2e-9;

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

    PS_full = wfm.PS_full;
    PS = wfm.PS;
    cPS = wfm.cPS;
    trPS = wfm.trPS;


    f_full = wfm.f_full;
    f = wfm.f;
    dt = wfm.TE/length(wfm.g);

    load UDSRTriN1000
    u = [UDSR.x UDSR.y UDSR.z];
    ODF = UDSR;

    bt = real(squeeze(sum(wfm.PS_full)));
    bt = bt/trace(bt);


    g = [gx gy gz];
    q = [qx qy qz];


    % ----- tuning with fminsearch(fun,x0,options) -----

    indR = find(Rs == R);             % index of reference radius in Rs

    % limit frequency range
    ind = find(trPS<PS_thresh);
    f = f(ind);
    PS = PS(ind,:,:);

    % interpolate to speed up
    [f, PS] = interpolate_power_spectra(f,PS, 1000);
    tracePS = PS(:,1,1)+PS(:,2,2)+PS(:,3,3);

    % normalize to trace
    trace_norm = sum(tracePS);
    tracePS = tracePS/trace_norm;
    PS = real(PS/trace_norm);


    % calculate tuning contours
    num_columns = 1;  % Number of columns
    plot_height = 0.14;
    plot_width = 0.35;
    vertical_spacing = 0.02;
    horizontal_spacing = 0.05;

    Ngrid = 200;
    theta = linspace(0, pi, Ngrid);
    phi = linspace(0, 2 * pi, Ngrid);
    num_plots = length(Rs_grid);

    tuning = zeros(num_plots, Ngrid, Ngrid);

    % Compute tuning proximity (1 - |1 - D_u/MD|) over polar and azimuthal angles for each sphere radius
    for nr = 1:num_plots
        normDw = DwSpherical(2 * pi * f, Rs_grid(nr), D0, 0, 50)' / D0;  % spherical restriction
        normMD = sum(tracePS .* normDw);      

        for n_theta = 1:length(theta)
            for n_phi = 1:length(phi)
                PSij = tensor_projection_ang(PS, theta(n_theta), phi(n_phi));
                %tuning(n_theta, n_phi) = 1 - abs(3 * sum(PSij .* normDw) - normMD);
                tuning(nr, n_theta, n_phi) = 1 - abs(1 - 3 * sum(PSij .* normDw) / normMD);
            end
        end
    end

    tuning_lim = [min(tuning(:)) 1];
    sprintf('min tuning = %.3f', tuning_lim(1))


     % plot tuning contours
    fh = figure;
    clf;
    fh.Position = [440 428 600 750];  % figure size
    fh.Color = 'white';

    Xticks = [0 90 180 270 360];
    Yticks = [0 90 180];

    col = linspace(1, 0, num_plots);
    alpha = linspace(1, .2, num_plots);
    col = [col; 0 * col; 1 - col];

    rows = ceil(num_plots / num_columns);  % rows

    X = phi * 180 / pi;
    Y = theta * 180 / pi;

    for nr = 1:num_plots

        % Determine column and row indices
        col_num = mod(nr-1, num_columns) + 1;  % Column number
        row_num = ceil(nr / num_columns);  % Row number

        % Manual positioning of subplots
        left = 0.1 + (col_num-1) * (plot_width + horizontal_spacing);  % based on column
        bottom = 1 - (row_num * (plot_height + vertical_spacing));  % based on row

        % axes with specific position
        ax = axes('Position', [left, bottom, plot_width, plot_height]);

        surf(X, Y, squeeze(tuning(nr,:,:)), 'EdgeColor', 'None', 'FaceAlpha', 1);
        xticks(Xticks);
        yticks(Yticks);

        daspect([1 1 1]);
        view([0 90]);
        clim(tuning_plot_lim);
        axis tight;

        % tick sizes and line width
        set(gca, 'LineWidth', 1.5, 'TickDir', 'Out', 'TickLength', .01*[1 1], 'FontSize', fs);

        % grid lines
        hold on;
        for n = 2:length(Xticks)-1
            plot3(Xticks(n)*[1 1], [Yticks(1) Yticks(end)], [1 1], 'w:', 'LineWidth', 1);
        end
        for n = 2:length(Yticks)-1
            plot3([Xticks(1) Xticks(end)], Yticks(n)*[1 1], [1 1], 'w:', 'LineWidth', 1);
        end
        hold off;

        % remove axes labels
        if row_num < rows
            set(gca, 'XTickLabel', []);
        end
        if col_num > 1
            set(gca, 'YTickLabel', []);
        end
    end

    colormap('copper');

    if save_fig
        print(fh, fullfile(fig_dir, 'tuning_landscapes_rev1.png'), '-dpng', '-r300');
    end


    
    % make the colorbar separately
    if (1)
        Ticks = 100 * [tuning_plot_lim(1):diff(tuning_plot_lim)/5:tuning_plot_lim(2)];
        Ticks = round(Ticks);

        fh = figure;
        clf;
        fh.Position = [440 428 600 750];  % figure size
        fh.Color = 'white';

        subplot('Position', [0.15, 0.1, plot_width, plot_height]);
        ax = surf(X, Y, squeeze(tuning(1,:,:)) * 100, 'EdgeColor', 'None', 'FaceAlpha', 1);
        view([0 90])
        axis tight
        colormap('copper')
        clim([min(Ticks), max(Ticks)])
 
        xticks([]);
        yticks([]);

        % Colorbar position
        %cb = colorbar('Position', [.91 0.1 0.03 0.8]); 
        cb = colorbar('northoutside'); % colorbar on top  
        cb.Position = [0.15, plot_height+0.11, plot_width, 0.025];
        set(cb, 'Ticks', Ticks, 'TickLabels', num2str(Ticks'), 'TickDirection', 'in', 'TickLength', [0.01 0.01], 'FontSize', fs)
       

        if save_fig
            print(fh, fullfile(fig_dir, 'tuning_colorbar_rev.png'), '-dpng', '-r300');
        end
    end

    

    % Compute angular deviation between tuned projections for each sphere radius and the reference radius 

    global optForTuning
    optForTuning.PS = PS;

    % find angles for the selected size
    optForTuning.normDw = DwSpherical(2*pi*f,Rs(indR),D0,0,50)'/D0; % spherical restriction
    optForTuning.normMD = sum(tracePS.*optForTuning.normDw);
    theta_phi0 = fminsearch(@tuningSphere,[0, 0]);
    u0 = [sin(theta_phi0(1)) * cos(theta_phi0(2)), sin(theta_phi0(1)) * sin(theta_phi0(2)), cos(theta_phi0(1))];
    cosA = [];

    for nR = 1:length(Rs)
        optForTuning.normDw = DwSpherical(2*pi*f,Rs(nR),D0,0,50)'/D0; % spherical restriction
        optForTuning.normMD = sum(tracePS.*optForTuning.normDw);
        theta_phi = fminsearch(@tuningSphere, theta_phi0);
        cosA(end+1) = u0 * [sin(theta_phi(1)) * cos(theta_phi(2)), sin(theta_phi(1)) * sin(theta_phi(2)), cos(theta_phi(1))]';
    end


    % plot angular distance
    A = acos(cosA)*180/pi;
    YLIM = [.99 1.01];

    fh = figure;,clf
    set(gcf,'color','white')
    hold on
    plot(Rs*1e6, cosA,'k-','LineWidth',2*lw)
    plot(R*1e6*[1 1], YLIM .*[min(cosA), max(cosA)], 'k:','LineWidth',lw)
    set(gca,'LineWidth',2*lw,'Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',2*fs)
    ylim(YLIM .*[min(cosA), max(cosA)])

    fh = figure;,clf
    set(gcf,'color','white')
    hold on
    plot(Rs*1e6, A,'k-','LineWidth',2*lw)
    plot(R*1e6*[1 1], YLIM .*[min(A), max(A)], 'k:','LineWidth',lw)
    set(gca,'LineWidth',2*lw,'Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',2*fs)
    ylim(YLIM .*[min(A), max(A)])

    color1 = [0 0 0];
    color2 = [.9 .6 0];

    fh = figure;,clf
    set(fh,'defaultAxesColorOrder',[color1; color2]);
    set(gcf,'color','white')
    hold on
    yyaxis left
    plot(Rs*1e6, cosA,'k-','LineWidth',lw*2,'Color',color1)
    plot(R*1e6*[1 1], YLIM .*[min(cosA), max(cosA)], 'k:','LineWidth',lw)
    ylim(YLIM .*[min(cosA), max(cosA)])
    yyaxis right
    plot(Rs*1e6, A,'k-','LineWidth',lw*2,'color',color2)
    set(gca,'LineWidth',2*lw,'Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',2*fs)
    ylim(YLIM .*[min(A), max(A)])

    if save_fig
        print(fh, fullfile(fig_dir, 'tuning_landscapes_rev2.png'), '-dpng', '-r300');
    end


end

%--------------------------------- FUNCTIONS ------------------------------

function tuning = tuningSphere(theta_phi)
% tuningSphere - objective function for fminsearch to optimize projection angles
% Inputs: polar and azimuthal angles in radians
% Output: deviation from target tuning (D_u/MD = 1)

global optForTuning
PSij = tensor_projection_ang(optForTuning.PS,theta_phi(1),theta_phi(2));
tuning = abs(3*sum(PSij.*optForTuning.normDw)/optForTuning.normMD - 1);
end

