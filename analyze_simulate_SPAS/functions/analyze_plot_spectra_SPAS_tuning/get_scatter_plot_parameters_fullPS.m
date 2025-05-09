function scatter_plot_parameters = get_scatter_plot_parameters_fullPS(scatter_plot_parameters);

scatter_plot_parameters.fullPS = scatter_plot_parameters.wfm;
%scatter_plot_parameters.fullPS.legendType = 'line'; %'patch' 'line'; % other
scatter_plot_parameters.fullPS.lwIn = 6;

scatter_plot_parameters.fullPS.legendText = {'XX','YY','ZZ','XY','YZ','XZ'};

scatter_plot_parameters.fullPS.colors(:,1) = 1.0*[1 0 0];%XX
scatter_plot_parameters.fullPS.colors(:,2) = 0.8*[0 1 0];%YY
scatter_plot_parameters.fullPS.colors(:,3) = 1.0*[0 0 1];%ZZ
scatter_plot_parameters.fullPS.colors(:,4) = 0.8*[1 1 0];%XY
scatter_plot_parameters.fullPS.colors(:,5) = 0.8*[0 1 1];%YZ
scatter_plot_parameters.fullPS.colors(:,6) = 0.8*[1 0 1];%XZ

scatter_plot_parameters.fullPS.showPatch = 0;
scatter_plot_parameters.fullPS.thresh = 0.995; %if 0 then use tensor_evolution_parameters.thresh_PSf
scatter_plot_parameters.fullPS.f0 = 0; %if >0 overrides thresh, if 0 then obtain from col_spectra_parameters.thresh

scatter_plot_parameters.fullPS.real_LineStyle = '-';
scatter_plot_parameters.fullPS.imag_LineStyle = ':';
scatter_plot_parameters.fullPS.legend_width = 0.18; 
scatter_plot_parameters.fullPS.frq_ticks = [];  % for frequency axis

if scatter_plot_parameters.fullPS.show_axes
    scatter_plot_parameters.fullPS.left_margin = 0.09;
    scatter_plot_parameters.fullPS.width = scatter_plot_parameters.wfm.width; %1-2*scatter_plot_parameters.fullPS.left_margin;
    scatter_plot_parameters.fullPS.height = scatter_plot_parameters.fullPS.width/1.618;
    scatter_plot_parameters.fullPS.bottom_margin = 0.15;
else
    scatter_plot_parameters.fullPS.left_margin = 0.02;
    scatter_plot_parameters.fullPS.width = 1-2*scatter_plot_parameters.fullPS.left_margin;
    scatter_plot_parameters.fullPS.height = scatter_plot_parameters.fullPS.width/1.618;
    scatter_plot_parameters.fullPS.bottom_margin = 0.05;
end
end
