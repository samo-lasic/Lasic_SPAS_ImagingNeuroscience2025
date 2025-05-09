function fcTensor = find_normalized_spectral_moments_from_projections(f,PSu, u, p);
% Compute the p-th spectral moment from spectral projections PSu along u.

PSu_norm = PSu./sum(PSu,2);

for n = 1:size(PSu,1)
    fc(n) = PSu_norm(n,:)*f';
end

fcTensor = tensor_from_projections(u,fc);

for ind = 1:3
    [azimuth,elevation,~] = cart2sph(fcTensor.V(1,ind), fcTensor.V(2,ind), fcTensor.V(3,ind));
    fcTensor.phi(ind) = mod(azimuth,pi);
    fcTensor.theta(ind) = pi/2 - elevation;
end

end

