function [fig, fh] = SPAS_map_fig_opt(n_map, size_nii, zoom_x, zoom_y, fig)

if isempty(fig.LimitSliceRange)
    fig.n_rows = ceil(size_nii(3)/fig.n_cols);
else
    fig.n_rows = ceil(diff(fig.LimitSliceRange)/fig.n_cols);
end

fig.sub_fig_width = (1 - fig.colormap_space) / fig.n_cols - fig.x_space;
fig.sub_fig_height = (1 - fig.title_size )/fig.n_rows - fig.y_space;

if ~fig.showTitle
    fig.title_size = 0;
end

fig.dimX = length(zoom_x);
fig.dimY = length(zoom_y);

fig.im_ratio = (fig.dimY*size_nii(2))/(fig.dimX*size_nii(1));


if n_map > 0
    fh = figure(fig.num + n_map - 1);
else
    fh = figure;
end

clf
set(gcf, 'InvertHardCopy', 'off')
if isfield(fig,'color')
    fh.Color = fig.color;
else
    fh.Color = 'Black';
end

set(fh, 'Units', 'inches', ...
    'PaperPosition',[0 0 fig.size_mm / 25.4 ], 'PaperPositionMode', 'auto');

screen_size = get(0,'screensize');
fh.Position(3:4) = fig.display_ratio * fig.size_mm .* screen_size(4);

