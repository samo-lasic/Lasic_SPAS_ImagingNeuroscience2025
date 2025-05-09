
function [bt, bT, f_full, f, PS_full, PS, cPS, trPS, gmat, qmat] = wfm_power_spectra_ensure_b(gmat, b, N0, TE)
% Scale gradient waveforms (gmat) to a specified b-value and computes cross power spectra normalized to b = 1

% bt (4D array [Nx3x3xNwaves]): Time-domain cumulative b-tensor (integrated over time) for each waveform.
% bT (3x3xNwaves array): Total b-tensor for each waveform
% f_full (vector [Nfreq]): Full frequency range, including positive and negative frequencies
% f (vector [Nfreq_positive]): Positive frequencies extracted from f_full
% PS_full (4D array [Nfreqx3x3xNwaves]): Cross power spectra for each waveform.
% PS (4D array [Nfreq_positivex3x3xNwaves]): Cross power spectra for each waveform. restricted to positive frequencies
% cPS (4D array [Nfreq_positivex3x3xNwaves]): Cumulative cross power spectra
% trPS (2D array [Nfreq_positivexNwaves]): Cumulative spectral trace 
% gmat (3D array [Nx3xNwaves]): Normalized gradient waveform, scaled to match the specified b-value.
% qmat (3D array [Nx3xNwaves]): Normalized dephasing, scaled to ensure correct diffusion weighting.


gmr = 26.75e7;

N = length(gmat);
t = linspace(0,TE,N);
dt = mean(diff(t));

qmat = gmr*cumsum(gmat)*dt;

for nWave = 1:size(gmat,3)

    g = squeeze(gmat(:,:,nWave));
    q = squeeze(qmat(:,:,nWave));

    [V, L] = eig(q'*q*dt);

    g = g/sqrt(trace(L))*sqrt(b);
    q = q/sqrt(trace(L))*sqrt(b);

    for i = 1:3
        for j = 1:3
            qq(:,i,j) = q(:,i).*q(:,j)*dt;
        end
    end

    bT(:,:,nWave) = sum(qq);
    bt(:,:,:,nWave) = cumsum(qq);

    q_zeropad = [zeros(N0,3); q; zeros(N0,3)];

    [f_full, PS_full] = power_spectra(q_zeropad,dt);

    PS_all(:,:,:,nWave) = PS_full;

    gmat(:,:,nWave) = g;
    qmat(:,:,nWave) = q;

    % find high frequency index (cumulative power up to 1-1e-6)
    ind = find(f_full>=0);
    f = f_full(ind);
    PS = PS_full(ind,:,:);
    P = cumsum(PS(:,1,1) + PS(:,2,2) + PS(:,3,3))/(sum(PS(:,1,1) + PS(:,2,2) + PS(:,3,3)));
    ind = max(find(P<1-1e-6));
    maxP(nWave) = P(ind); % close to 1
    f0(nWave) = f(ind);

end

% eliminate some super high frequency points
f0 = max(f0);

ind = find(f_full>=-f0 & f_full<=f0);
PS_full = PS_all(ind,:,:,:);
f_full = f_full(ind);

ind = find(f_full>=0);
f = f_full(ind);
PS = PS_full(ind,:,:,:);

% ensure b = 1
for nWave = 1:size(gmat,3)
    trB = trace(real(squeeze(sum(PS(:,:,:,nWave))))); % total power

    PS(:,:,:,nWave) = maxP(nWave)*PS(:,:,:,nWave)/trB;

    cPS(:,:,:,nWave) = cumsum(PS(:,:,:,nWave));
    trPS(:,nWave) = cPS(:,1,1,nWave) + cPS(:,2,2,nWave) + cPS(:,3,3,nWave);

    trB = trace(real(squeeze(sum(PS_full(:,:,:,nWave))))); % total power
    PS_full(:,:,:,nWave) = maxP(nWave)*PS_full(:,:,:,nWave)/trB;

end
end
