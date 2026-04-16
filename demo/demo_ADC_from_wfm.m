% Demo:
% Calculates ADC for different restricted geometries and diffusion encoding waveforms.
% Uses the effective gradient waveform g(t), with raster time dt set explicitly.
% Obtains g(t), q(t), and the normalized encoding power spectrum for a target
% b-value, and computes ADC from the spectrum using the Gaussian phase approximation.

clear all
restoredefaultpath

% --- paths ---
addpath(genpath(fullfile(pwd, '..', 'analyze_simulate_SPAS', 'functions')));

waveform_dir = fullfile('..', 'waveforms','M0_SPAS_2x21ms');
waveform_names{1} = 'Bruker7T_NOW_M0_nogap_10_6_interp_g_te21_flip0';

wfm = load(fullfile(waveform_dir,[waveform_names{1} '.mat']));
g = wfm.g; % effective gradient waveform [T/m]

% --- waveform + spectrum settings ---
N0       = 1e4;       % zero padding for power spectra
b_target = 1e9;       % [s/m^2]

dt = 2e-5;            % [s] Bruker raster used for these waveforms
Nt = size(g,1);
TE = Nt * dt;         % [s]
t  = (0:Nt-1) * dt;   % [s]

N0 = 10000; % zero padding power spectra
b_target = 1e9;

% If you store the full struct with true gradient and RF sign:
% g = wfm.g .* wfm.rf;
% dt = wfm.dt;
% TE = size(g,1) * dt;

t = linspace(0,1,Nt) * TE;

% Scale g to hit b_target and compute PS_full of q(t)
[~, ~, f_full, ~, PS_full, ~, ~, ~, g, q] = ...
    wfm_power_spectra_ensure_b(g, b_target, N0, TE); % use effective gradient here 


% ---------- plots ----------
if (1)
    figure('Color','w');

    % 1) g(t)
    subplot(3,1,1);
    plot(t*1e3, g*1e3, 'LineWidth', 1.2); grid on
    xlabel('Time [ms]'); ylabel('g [mT/m]'); title('Gradient waveform g(t)');

    % 2) q(t)
    subplot(3,1,2);
    plot(t*1e3, q, 'LineWidth', 1.2); grid on
    xlabel('Time [ms]'); ylabel('q'); title('Dephasing q(t)');

    % 3) trace power spectrum vs frequency
    subplot(3,1,3);
    trPS_full = real(PS_full(:,1,1) + PS_full(:,2,2) + PS_full(:,3,3));
    plot(f_full, trPS_full, 'LineWidth', 1.2); grid on
    xlabel('Frequency [Hz]'); ylabel('(PS)'); title('Power spectrum (trace)');
    xlim([min(f_full) max(f_full)]);
end

% Package waveform info for ADC calculator
wfmInfo.f_full = f_full;
wfmInfo.PS_full = PS_full;
wfmInfo.g       = g;
wfmInfo.q       = q;

% --- substrate ---
select_substrate = 2;
D0 = 1e-9;         % [m^2/s]

% size R in [m]
switch select_substrate
    case 1  % sphere: Pin = [R, D0]
        auxP.shape = 1;
        R = 10e-6;                       
        Pin = [R, D0];                   
        
    case 2  % cylinder (axial free): Pin = [R, D0]
        auxP.shape = 2;
        R = 5e-6;                        
        Pin = [R, D0];                   
        
    case 3  % ellipsoid/spheroid: Pin = [R1, R2, D0]
        auxP.shape = 3;
        R1 = 4e-6; R2 = 8e-6;           
        Pin = [R1, R2, D0];             
        
    case 4  % stick (planar restriction axially): Pin = [R, D0]
        auxP.shape = 4;
        R = 1e-6;                         
        Pin = [R, D0];                    
        
    case 5  % free diffusion (possibly anisotropic): Pin = [D3, D2-D3, D1-D2]
        auxP.shape = 5;
        Pin = [D0, 0, 0];
        
    case 6  % asymmetric cylinder (two radial radii): Pin = [R1, R2, D0]
        auxP.shape = 6;
        R1 = 3e-6; R2 = 6e-6;            
        Pin = [R1, R2, D0];              
end
auxP.tortuosity = 0;   % alpha in Dw - use 0 for a single impermeable compartment

% Orientation angles: fADCfromWFM needs cos/sin of two angles (phi,theta).
auxP.cP = 1; auxP.sP = 0;   % cos(phi), sin(phi)
auxP.cT = 1; auxP.sT = 0;   % cos(theta), sin(theta)

ADC = fADCfromWFM(Pin, wfmInfo, auxP);

% Example: electrostatic-repulsion orientation set (UDSRTriN*)
if (0)
    load UDSRTriN15   % UDSRTriN15, UDSRTriN300
    auxP.cP = cos(UDSR.phi(:));
    auxP.sP = sin(UDSR.phi(:));
    auxP.cT = cos(UDSR.theta(:));
    auxP.sT = sin(UDSR.theta(:));
    ADC = fADCfromWFM(Pin, wfmInfo, auxP)
end






