

function exGamma = exchange_gamma_tensorial(q,b,dt,u)

% Tensorial version of Exchange Gamma:
% Calculate Gamma components that can be projected as a tensor using i, j, k, l indices of unit vectors.

% q4(t) = integral_0^tau q^2(x) * q^2(x+t) dx
% Gamma = 2 * integral_0^tau q4(t) * t dt

% References:
% - Ning L, Nilsson M, Lasic S, Westin C-F, Rathi Y. "Cumulant expansions for measuring water exchange using diffusion MRI." J Chem Phys 2018; 148.
% - Szczepankiewicz F, Westin CF, Nilsson M. "Gradient waveform design for tensor-valued encoding in diffusion MRI." J Neurosci Methods 2020; 109007.

% Exchange Gamma in Voigt notation:
% - exGamma_ij: 6x6 matrix in Voigt notation.
% - exGamma_u = sum_ijkl u_i * u_j * u_k * u_l * integral_s q_i(s) * q_j(s) * q_k(t-s) * q_l(s-t) ds.
% - exGamma_u = sum_ijkl u_i * u_j * u_k * u_l * exGamma4_ijkl.

% In 6x6 notation:
% - exGamma_u = sum_ij uu_ij * exGamma4_ij, where i, j = 11, 22, 33, 2x23, 2x13, 2x12.

% Definitions:
% - u6 = [u(1)^2, u(2)^2, u(3)^2, 2*u(2)*u(3), 2*u(1)*u(3), 2*u(1)*u(2)].
% - u6x6 = u6 * u6.


N = size(q,1);
t = [0:dt:(N-1)*dt]';

uu(:,1) = u(:,1).*u(:,1);
uu(:,2) = u(:,2).*u(:,2);
uu(:,3) = u(:,3).*u(:,3);
uu(:,4) = 2*u(:,2).*u(:,3);
uu(:,5) = 2*u(:,1).*u(:,3);
uu(:,6) = 2*u(:,1).*u(:,2);

qq(:,1) = q(:,1).*q(:,1);
qq(:,2) = q(:,2).*q(:,2);
qq(:,3) = q(:,3).*q(:,3);
qq(:,4) = q(:,2).*q(:,3);
qq(:,5) = q(:,1).*q(:,3);
qq(:,6) = q(:,1).*q(:,2);

for m = 1:6
    for n = 1:6
        qq1 = qq(:,m);
        qq2 = qq(:,n);

        c = xcorr(qq1,qq2)*dt;
        q4 = c(N:end);

        exGamma6x6(m,n) = 2*sum(q4.*t)*dt/b^2;
    end
end


Nu = size(u,1);
for n = 1:Nu
    uu1 = uu(n,:);
    exGamma(n) = sum(sum((uu1'*uu1).*exGamma6x6));
end


end
