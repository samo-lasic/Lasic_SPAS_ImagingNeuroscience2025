

function exGamma = exchange_gamma_projection(q,b,dt,u)
% Calculate Gamma for all projections qu (Nt x Ndirection).
% q4(t) = integral_0^tau q^2(x) * q^2(x+t) dx
% Gamma = 2 * integral_0^tau q4(t) * t dt
% This is the 1D version of Gamma4.

% References:
% - Ning L, Nilsson M, Lasic S, Westin C-F, Rathi Y. "Cumulant expansions for measuring water exchange using diffusion MRI." J Chem Phys 2018; 148.
% - Szczepankiewicz F, Westin CF, Nilsson M. "Gradient waveform design for tensor-valued encoding in diffusion MRI." J Neurosci Methods 2020; 109007.


q2 = q*u';
q2 = q2.^2;

N = length(q2);
t = [0:dt:(N-1)*dt]';


for n = 1:size(q2,2)
    c = xcorr(q2(:,n),q2(:,n))*dt;
    q4(:,n) = c(N:end);
end

exGamma = 2*sum(q4.*t)*dt./b.^2;

end

