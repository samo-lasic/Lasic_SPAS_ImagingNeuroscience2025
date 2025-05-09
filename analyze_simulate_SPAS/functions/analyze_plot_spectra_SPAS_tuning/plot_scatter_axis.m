
function plot_scatter_axis(X,Y,x_label,y_label,col,line_styles,parameters)
% ----   patches between 0 and Y   -------

for n = 1:size(Y,2)
    if numel(line_styles) == 1
        line_style = char(line_styles);
    else
        line_style = line_styles{n};
    end

    if (parameters.showPatch)
        [px,py] = patch_XY(X,0*Y(:,n),Y(:,n),parameters.Npatch);
        p = patch(px,py,'k');
        p.FaceColor = col(:,n);
        p.FaceAlpha = parameters.patch_face_alpha;
        p.EdgeColor = 'none';
        p.LineStyle = line_style;

    end

end


for n = 1:size(Y,2)
    if numel(line_styles) == 1
        line_style = char(line_styles);
    else
        line_style = line_styles{n};
    end

    [X1,Y1] = interpolate_plot(X,Y(:,n),parameters.Ninterpolate);
    plot(X1,Y1,line_style,'Color',col(:,n),'LineWidth',parameters.lwIn,'MarkerSize',10)
end


axis tight
set(gca,'LineWidth',parameters.lwOut,'Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',parameters.fs)
if ~parameters.show_axes
    axis off
elseif ~parameters.show_Yticks
    set(gca,'ytick',[])
end

if parameters.show_labels
    xlabel(x_label);
    ylabel(y_label);
end
end

