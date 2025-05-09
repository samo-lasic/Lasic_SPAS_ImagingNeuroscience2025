function [fig, fh] = fig_opt(n_fig, fig)

if n_fig > 0
    fh = figure(n_fig);
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

