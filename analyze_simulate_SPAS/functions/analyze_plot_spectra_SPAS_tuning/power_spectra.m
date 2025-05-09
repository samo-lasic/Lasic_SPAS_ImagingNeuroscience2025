function [f, PS] = power_spectra(F, dt, mode)
% - Compute the cross power spectra for F (NT x N) 
% - Two modes:
%   1. 'positive' mode: spectra for positive frequencies only.
%   2. Full mode: spectra for both positive and negative frequencies.


N = size(F,2);
NT = size(F,1);

f = 1/dt*(0:(NT/2))/NT;

if (exist('mode', 'var'))
    if strcmp(mode,'positive')
        mode = 1;
        PS = zeros(length(f),N,N);
    end
else
    mode = 0;
end

if ~mode
    f = [-fliplr(f(2:end)) f(1:end)];
    PS = zeros(NT+1,N,N);
end


for m = 1:N
    for n = 1:N

        ps = fft(F(:,m)).*conj(fft(F(:,n)));

        if mode
            PStmp = ps(1:length(f));
        else
            PStmp = fftshift(ps);
            PStmp = [PStmp; PStmp(1)];
        end

        PS(:,m,n) = PStmp;

    end
end

end

