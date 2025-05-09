function correlate_maps_struct = SPAS_make_correlate_maps_struct(correlate_maps_struct, ...
    correlate_prepend, correlate_name, correlate_append, correlate_range, correlate_gamma_cor, ...
    modulate_prepend, modulate_name, modulate_append, modulate_range, modulate_gamma_cor, color_order, title_str)

% used in ..._show_correlation_maps.m

correlate_maps_struct(end+1).correlate_prepend = correlate_prepend;
correlate_maps_struct(end).correlate_name = correlate_name;
correlate_maps_struct(end).correlate_append = correlate_append;
correlate_maps_struct(end).correlate_range = correlate_range;
correlate_maps_struct(end).correlate_gamma_cor = correlate_gamma_cor;

correlate_maps_struct(end).modulate_prepend = modulate_prepend;
correlate_maps_struct(end).modulate_name = modulate_name;
correlate_maps_struct(end).modulate_append = modulate_append;
correlate_maps_struct(end).modulate_range = modulate_range;
correlate_maps_struct(end).modulate_gamma_cor = modulate_gamma_cor;

correlate_maps_struct(end).color_order = color_order;
correlate_maps_struct(end).title_str = title_str;

end
