function make_axes1(parameters, spectral_PAS, L, ms, order)
% Create axes for plotting.
% Expects the "filteredB_shape" field, e.g., spectral_PAS.filteredB_shape.FA.
% Note: The "filteredB_shape" field is removed in make_axes2.


if nargin < 5
    order = [1 2 3];
end
hold on

if nargin > 1
    if parameters.showSPAS % spectral PAS (SPAS)

        if abs(spectral_PAS.filteredB_shape.FA - 1) < 1e-10
            v = spectral_PAS.filteredB_shape.V(:,1);
            plot3(v(1)*L*[-1:1],v(2)*L*[-1:1],v(3)*L*[-1:1],':','Color', parameters.SPAS_color(1,:),'Linewidth',parameters.ax_width)

        else

            trail = 0.9*linspace(-1,1,500);

            v = spectral_PAS.filteredB_shape.V(:,order(1));
            plot3(v(1)*L*trail,v(2)*L*trail,v(3)*L*trail,'.','Color',parameters.SPAS_color(:,1),'MarkerSize',ms)
            v = spectral_PAS.filteredB_shape.V(:,order(2));
            plot3(v(1)*L*trail,v(2)*L*trail,v(3)*L*trail,'.','Color',parameters.SPAS_color(:,2),'MarkerSize',ms)
            v = spectral_PAS.filteredB_shape.V(:,order(3));
            plot3(v(1)*L*trail,v(2)*L*trail,v(3)*L*trail,'.','Color',parameters.SPAS_color(:,3),'MarkerSize',ms)

        end
    end
end
end