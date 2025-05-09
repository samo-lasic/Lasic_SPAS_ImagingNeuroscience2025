
function scatter_plot_parameters = get_scatter_plot_parameters_wfm(scatter_plot_parameters)

scatter_plot_parameters.wfm.showPatch = 1; % what is shown in the legend, lines or patches?  (only non-zero components are shown)


if scatter_plot_parameters.wfm.normalizeXaxis
    scatter_plot_parameters.wfm.t_label = 'time [x 1/\tau]';
else
    scatter_plot_parameters.wfm.t_label = 'time [ms]';
end

switch scatter_plot_parameters.wfm.units_g
    case 1
        scatter_plot_parameters.wfm.scale_g = 1;
        scatter_plot_parameters.wfm.g_label = '\bfg\rm [T/m]';
    case 2
        scatter_plot_parameters.wfm.scale_g = 1e3;
        scatter_plot_parameters.wfm.g_label = '\bfg\rm [mT/m]';
end

switch scatter_plot_parameters.wfm.units_q
    case 1
        scatter_plot_parameters.wfm.scale_q = 1;
        scatter_plot_parameters.wfm.q_label = '\bfq\rm [1/m]';
    case 2
        scatter_plot_parameters.wfm.scale_q = 1e-3;
        scatter_plot_parameters.wfm.q_label = '\bfq\rm [1/mm]';
    case 3
        scatter_plot_parameters.wfm.scale_q = 1e-6;
        scatter_plot_parameters.wfm.q_label = '\bfq\rm [1/\mum]';
end

scatter_plot_parameters.wfm.trPS_label = 'trace(\bfs\rm)';
scatter_plot_parameters.wfm.f_label = 'frequency [Hz]';
scatter_plot_parameters.wfm.PS_label = 'power spectra';

scatter_plot_parameters.wfm.color(:,1) = 1.0*[1 0 0];%XX
scatter_plot_parameters.wfm.color(:,2) = 0.8*[0 1 0];%YY
scatter_plot_parameters.wfm.color(:,3) = 1.0*[0 0 1];%ZZ
scatter_plot_parameters.wfm.legendText = {'X','Y','Z'};

scatter_plot_parameters.wfm.Npatch = 1; % avery Np-th point is used for pathing (sometimes can speed up things)

scatter_plot_parameters.wfm.Ninterpolate = 1e5; %1e5; % if > 0 interpolate all scatter plots
scatter_plot_parameters.wfm.lwIn = 4;
scatter_plot_parameters.wfm.lwOut = 4;
scatter_plot_parameters.wfm.fs = 18;
scatter_plot_parameters.wfm.patch_face_alpha = 0.5;

scatter_plot_parameters.wfm.aspectRatio = 1.618*4;%1.618;


if scatter_plot_parameters.wfm.show_axes
    scatter_plot_parameters.wfm.left_margin = 0.09;
    scatter_plot_parameters.wfm.bottom_margin = 0.1; %0.12;
    scatter_plot_parameters.wfm.width = 0.88;
    scatter_plot_parameters.wfm.height = 0.19;
    scatter_plot_parameters.wfm.yspacing1 = 0.11; % 0.07; % top
    scatter_plot_parameters.wfm.yspacing2 = 0.07; % 0.07; % bottom

    if scatter_plot_parameters.wfm.show_labels
        scatter_plot_parameters.wfm.left_margin = 0.14;
        scatter_plot_parameters.wfm.bottom_margin = 0.15; %0.12;
        scatter_plot_parameters.wfm.width = 0.82;
        scatter_plot_parameters.wfm.height = 0.18;
        scatter_plot_parameters.wfm.yspacing1 = 0.11;% 0.07; % top
        scatter_plot_parameters.wfm.yspacing2 = 0.07; %0.07; % bottom
    end

else
    scatter_plot_parameters.wfm.left_margin = 0.025;
    scatter_plot_parameters.wfm.bottom_margin = 0.05;
    scatter_plot_parameters.wfm.width = 0.95;
    scatter_plot_parameters.wfm.height = 0.25;
    scatter_plot_parameters.wfm.yspacing1 = 0.05;
    scatter_plot_parameters.wfm.yspacing2 = 0.00;
end

% -----------------------------

scatter_plot_parameters.wfm.remove_axis_g = 0;
scatter_plot_parameters.wfm.remove_axis_q = 0;
scatter_plot_parameters.wfm.spaceY_g = 0;
scatter_plot_parameters.wfm.spaceX_g = 0;
scatter_plot_parameters.wfm.spaceY_q = 0;
scatter_plot_parameters.wfm.spaceX_q = 0;

if scatter_plot_parameters.wfm.show_axes
    scatter_plot_parameters.wfm.scale_size_g = 1.14;
    scatter_plot_parameters.wfm.left_margin_g = -0.82;
    scatter_plot_parameters.wfm.bottom_margin_g = -0.85;

    scatter_plot_parameters.wfm.scale_size_q = scatter_plot_parameters.wfm.scale_size_g;
    scatter_plot_parameters.wfm.left_margin_q = scatter_plot_parameters.wfm.left_margin_g;
    scatter_plot_parameters.wfm.bottom_margin_q = -1.18;

    if scatter_plot_parameters.wfm.show_labels
        scatter_plot_parameters.wfm.scale_size_g = 1.1;
        scatter_plot_parameters.wfm.left_margin_g = -0.72;
        scatter_plot_parameters.wfm.bottom_margin_g = -0.767;

        scatter_plot_parameters.wfm.scale_size_q = scatter_plot_parameters.wfm.scale_size_g;
        scatter_plot_parameters.wfm.left_margin_q = scatter_plot_parameters.wfm.left_margin_g;
        scatter_plot_parameters.wfm.bottom_margin_q = -1.07;
    end

else
    scatter_plot_parameters.wfm.scale_size_g = 1.25;
    scatter_plot_parameters.wfm.left_margin_g = -0.95;
    scatter_plot_parameters.wfm.bottom_margin_g = -0.95;

    scatter_plot_parameters.wfm.scale_size_q = scatter_plot_parameters.wfm.scale_size_g;
    scatter_plot_parameters.wfm.left_margin_q = scatter_plot_parameters.wfm.left_margin_g;
    scatter_plot_parameters.wfm.bottom_margin_q = -1.25;
end

end
