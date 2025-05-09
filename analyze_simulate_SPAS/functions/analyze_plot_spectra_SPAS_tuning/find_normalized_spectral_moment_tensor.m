function m = find_normalized_spectral_moment_tensor(f, PS, cum_tracePS, p);
% Find the normalized p-th spectral moment:
% - p = 2 for Vw (variance),
% - p = 1 for the centroid.


b = trace(squeeze(real(sum(PS))));
m.filteredB = real(squeeze(sum(PS .* repmat(f',1, 3, 3).^p )))/b; % assuming f >= 0

% [V L] = eig(gB);
m.filteredB_shape = tensor_shape(real(m.filteredB));

for ind = 1:3
    [azimuth,elevation,~] = cart2sph(m.filteredB_shape.V(1,ind),m.filteredB_shape.V(2,ind),m.filteredB_shape.V(3,ind));
    m.phi(ind) = mod(azimuth,pi);
    m.theta(ind) = pi/2 - elevation;
end

end
