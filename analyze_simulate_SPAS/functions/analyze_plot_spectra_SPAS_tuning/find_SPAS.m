function spectral_PAS = find_SPAS(f,PS,cum_tracePS,LF_power_thresh)
% find spectral PAS based on filtered b-tensor

spectral_PAS.thresh_index = max(find(cum_tracePS<LF_power_thresh*cum_tracePS(end)));
spectral_PAS.fthresh = f(spectral_PAS.thresh_index);
spectral_PAS.filteredB = real(squeeze(sum(PS(1:spectral_PAS.thresh_index,:,:))));

% [V L] = eig(gB);
spectral_PAS.filteredB_shape = tensor_shape(real(spectral_PAS.filteredB));

for ind = 1:3
    [azimuth,elevation,~] = cart2sph(spectral_PAS.filteredB_shape.V(1,ind),spectral_PAS.filteredB_shape.V(2,ind),spectral_PAS.filteredB_shape.V(3,ind));
    spectral_PAS.phi(ind) = mod(azimuth,pi);
    spectral_PAS.theta(ind) = pi/2 - elevation;
end

end
