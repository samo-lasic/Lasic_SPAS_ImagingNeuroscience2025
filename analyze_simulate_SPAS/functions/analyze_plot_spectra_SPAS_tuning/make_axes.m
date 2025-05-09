
function make_axes(parameters, spectral_PAS)
% help plotting

hold on
if (0)
    L = 1.5;
    make_cylindrical_axes(L,parameters.ax_width*L/200,20,parameters.ax_style(1),parameters.ax_style(2),parameters.ax_color)
else

    L = 1.5;

    trail = L * 0.9*linspace(-1,1,500);
    l = 2*max(trail)/round(3+2*25/parameters.ax_width);
    trail = mod(floor(trail/l),2) .* trail;
    trail(trail == 0) = nan;

    ms = 12 * parameters.ax_width;

    v = [1 0 0];
    plot3(v(1)*trail,v(2)*trail,v(3)*trail,'.','Color',parameters.ax_color(1,:),'MarkerSize',ms)
    v = [0 1 0];
    plot3(v(1)*trail,v(2)*trail,v(3)*trail,'.','Color',parameters.ax_color(2,:),'MarkerSize',ms)
    v = [0 0 1];
    plot3(v(1)*trail,v(2)*trail,v(3)*trail,'.','Color',parameters.ax_color(3,:),'MarkerSize',ms)

end

if nargin > 1
    if parameters.showSPAS % spectral PAS (SPAS)
        hold on
        L = 1.5;

        if abs(spectral_PAS.filteredB_shape.FA - 1) < 1e-10
            v = spectral_PAS.filteredB_shape.V(:,1);
            plot3(v(1)*L*[-1:1],v(2)*L*[-1:1],v(3)*L*[-1:1],':','Color', parameters.SPAS_color(1,:),'Linewidth',parameters.ax_width)

        else

            trail = 0.9*linspace(-1,1,500);

            ms = 24;
            v = spectral_PAS.filteredB_shape.V(:,1);
            plot3(v(1)*L*trail,v(2)*L*trail,v(3)*L*trail,'.','Color',parameters.SPAS_color(:,1),'MarkerSize',ms)
            v = spectral_PAS.filteredB_shape.V(:,2);
            plot3(v(1)*L*trail,v(2)*L*trail,v(3)*L*trail,'.','Color',parameters.SPAS_color(:,2),'MarkerSize',ms)
            v = spectral_PAS.filteredB_shape.V(:,3);
            plot3(v(1)*L*trail,v(2)*L*trail,v(3)*L*trail,'.','Color',parameters.SPAS_color(:,3),'MarkerSize',ms)


        end
    end
end
end

