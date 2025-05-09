function fig = SPAS_tiling_fig_opt(c)

if nargin < 1
    c = 1;
end

fig.num = 1;

fig.show_rois = 0; % if not present, don't show roi
fig.roi_lw = 2;

fig.display_ratio = 1.2e-4;

fig.showTitle = 1;
fig.title_size = 0.1;
fig.title_fs = 12;



%LPS (right to left)
if ~isfield(fig,'reverseX')
    fig.reverseX = 1;
end
if ~isfield(fig,'reverseY')
    fig.reverseY = 0;
end

fig.save = 0;
fig.resolution = '-r300'; % '-r200'; '-r300';


% used for signal contrast maps
fig.show_roi_scatter = 1; % if not present, don't show roi
opt.fig.show_colormap = 1;

fig.colormap_space = .15; %.12;
fig.x_space = 0.01; % 0.02
fig.y_space = 0;


switch c
    case 1

        % -----------------     tiling     ----------------------

        fig.size_mm = [80 36]; % 6 slices
        fig.n_cols = 3;
        fig.LimitSliceRange = []; % e.g. [], [2 10] if empty show all slice else limit min and max slice index

        %     fig.size_mm = [132 36]; % 10 slices
        %     fig.n_cols = 5;
        %     fig.LimitSliceRange = []; % e.g. [], [2 10] if empty show all slice else limit min and max slice index

        %     fig.size_mm = [80 60]; % 10 slices
        %     fig.n_cols = 3;
        %     fig.LimitSliceRange = [2 10]; % e.g. [], [2 10] if empty show all slice else limit min and max slice index

    case 2
        fig.size_mm = [80 14]; % 4 slices
        fig.n_cols = 4;
        fig.LimitSliceRange = [2 5]; % e.g. [], [2 10] if empty show all slice else limit min and max slice index

    case 3
        fig.size_mm = [60 35]; %[63 36]; % 4 slices
        fig.n_cols = 2;
        fig.LimitSliceRange = [2 5]; % e.g. [], [2 10] if empty show all slice else limit min and max slice index

    case 4
        fig.size_mm = [60 35]; %[63 36]; % 4 slices
        fig.n_cols = 2;
        fig.LimitSliceRange = [1 4]; % e.g. [], [2 10] if empty show all slice else limit min and max slice index

    case 5
        fig.size_mm = [70 35];  % 6 slices
        fig.n_cols = 3;
        fig.LimitSliceRange = [1 6]; % e.g. [], [2 10] if empty show all slice else limit min and max slice index

    case 6
        fig.size_mm = [120 22]; %[63 36]; % 4 slices

        fig.LimitSliceRange = [2 5]; % e.g. [], [2 10] if empty show all slice else limit min and max slice index
        fig.n_cols = 1+diff(fig.LimitSliceRange);

        fig.colormap_space = .08; %.12;
        fig.x_space = 0.01; % 0.02
        fig.y_space = 0;

    case 7
        fig.size_mm = [140 22]; % 6 slices

        fig.LimitSliceRange = [1 6]; % e.g. [], [2 10] if empty show all slice else limit min and max slice index
        fig.n_cols = 1+diff(fig.LimitSliceRange);

        fig.colormap_space = .06; %.12;
        fig.x_space = 0.00; % 0.02
        fig.y_space = 0;

end




