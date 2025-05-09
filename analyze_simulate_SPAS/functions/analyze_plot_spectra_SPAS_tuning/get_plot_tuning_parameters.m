function plot_tuning_parameters = get_plot_tuning_parameters(tuning_parameters);

% ----  plot contor lines for one Dw
plot_tuning_parameters = tuning_parameters;

plot_tuning_parameters.contour.limits = [0.1 0.8]; % plot contours, if [0 1] then full range
plot_tuning_parameters.contour.Nlevels = 5;% plot contours if 0 no contour lines
plot_tuning_parameters.contour.log = 0;% plot contours logarithmic scale
plot_tuning_parameters.contour.Nres = 500;% plot contours, resample 300, 500

% ---- plot band for the lowest contour line
plot_tuning_parameters.contour.do_shade = 1;
plot_tuning_parameters.contour.shade.band = 5e-2; % if 0 no band
plot_tuning_parameters.contour.shade.amp = 0.45; % 0-1 amplitude
plot_tuning_parameters.contour.shade.p = 2; % plot band
plot_tuning_parameters.contour.shade.color = [0 1 0];

plot_tuning_parameters.plot.colormap_limits = [];
plot_tuning_parameters.plot.patch_brightness = 1; %0.8; for showing scale

plot_tuning_parameters.plot.show_trajectory = 1;
plot_tuning_parameters.plot.trajectory_scale = 1.15;
plot_tuning_parameters.plot.trajectory_lw = 4;
plot_tuning_parameters.plot.ax_style = [20,0.5]; % number of segments, ratio length/spacing of segments (0-1: 0.5 is even)
plot_tuning_parameters.plot.ax_width = 2;
plot_tuning_parameters.plot.ax_color = eye(3);

plot_tuning_parameters.plot.edge_brightness = 0.5; %0.8;
plot_tuning_parameters.plot.edge_alpha = 0;%0.5
plot_tuning_parameters.plot.face_alpha = 0.88;
plot_tuning_parameters.plot.ambient = 0.95;
plot_tuning_parameters.plot.specular = 0.3;

plot_tuning_parameters.plot.az = -30;
plot_tuning_parameters.plot.el = 20;
plot_tuning_parameters.plot.lightPos1 = 1*[1 -0.5 .5];
plot_tuning_parameters.plot.lightPos2 = 0*[0 1 0.5];

plot_tuning_parameters.plot.CameraViewAngle = 6; % controls zoom

plot_tuning_parameters.plot.showBshape = 1;
plot_tuning_parameters.plot.showSQ = 0;
plot_tuning_parameters.plot.SQoffset = 0.18;

plot_tuning_parameters.plot.show_optTuned = 1;
plot_tuning_parameters.plot.tunedLTE_color = [1 .8 0];
plot_tuning_parameters.plot.optLTE_color = [1 .5 0];
plot_tuning_parameters.plot.LTE_MarkerSize = 32;


end

