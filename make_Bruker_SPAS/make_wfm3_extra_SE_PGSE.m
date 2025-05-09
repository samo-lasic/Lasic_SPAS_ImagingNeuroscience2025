% Generates a PGSE (pulsed gradient spin echo) waveform with a gap for the 180° RF pulse, followed by an appended polarity-flipped copy.
% Saves the resulting gradient waveform (3D vector) as a .mat file.


clear all
% close all
warning ('off','all');

gmr = 2.6751e+08;

save_path = '.../waveforms/PGSE_te21.mat'; % example path, adjust as needed

G = 0.6; % will be rescaled to match other waveforms, e.g. SPAS
N = 1051; % make sure the number of points is the same as in the original SPAS waveform (before additing 180 gap)

dt_exp = 2e-5;

t_ramp = 0.5e-3;
t_plato = 2e-3;

N_ramp = round(t_ramp/dt_exp);
N_plato = round(t_plato/dt_exp);

g = G*[linspace(0,1,N_ramp) ones(1,N_plato) linspace(1,0,N_ramp) zeros(1, N - N_plato - 2*N_ramp)]';
g = [g*0 g*0 g];

dur180Gap = 5.05e-3;
flipPolarityAfter180 = 0; % 1 or 0 to flip polarity after adding 180 gap

Ngap = round(dur180Gap/dt_exp);
display(sprintf('gap time is %.3f ms', Ngap*dt_exp*1000))
gap = zeros(Ngap,3);

g = [g; gap; -flipud(g)];

save(save_path,'g')

q = gmr*cumsum(g)*dt_exp;
b = trace(q'*q)*dt_exp;
maxG = max(abs(g(:)));
slew = max(max(abs(diff(g))))/dt_exp;

N180 = round(length(q)/2);
qc = q(N180,:);
qc = sqrt(qc*qc'); % magnitude of crushing vector

disp_str = sprintf('Gmax = %.3f, slew = %g, |qc| = %g, b = %g',maxG,slew,qc,b);
display(disp_str)

return
% PLOT RESULTS
t = [0:length(g)-1]*dt_exp;
m1 = gmr*cumsum(g.*t')*dt_exp;

figure%,clf
subplot(3,1,1)
plot(t,g,'.')
xlabel('Time [ms]')
ylabel('Gradient amplitude [mT/m]')
% title(disp_str, 'Interpreter', 'none');

subplot(3,1,2)
plot(t,m1,'.')
ylabel('m1')

subplot(3,1,3)
plot(t,q,'.')
xlabel('Time [ms]')
ylabel('q')

