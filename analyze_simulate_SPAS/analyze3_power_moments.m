% a helper script to check spectral moments like spectral variance (Vw) or centroid frequency
% from power spectra stored in waveform info MAT files. Requires files with the wfm structure containing f and PS fields.

clear all

restoredefaultpath
addpath(genpath(fullfile(pwd,'functions')));

% path to waveforms 
waveform_dir = fullfile('..', 'waveforms', 'M1_SPAS_2x21ms_rot');

% use _info files (with spectral info)
base_wfm_name = 'STE1_21000_5040_20_1051_252_220829_rot'; % no extension

waveform_names = {...
    strcat(base_wfm_name, '_info.mat'),...
    strcat(base_wfm_name, '_SPAS1_info.mat'), ...
    strcat(base_wfm_name, '_SPAS2_info.mat'), ...
    strcat(base_wfm_name, '_SPAS3_info.mat'), ...
    strcat(base_wfm_name, '_tunedLTE_res_info.mat'), ...
    };


calc_spectral_moments(1, waveform_dir, waveform_names);
calc_spectral_moments(2, waveform_dir, waveform_names);

% --- FUNCTIONS ---

function calc_spectral_moments(p, waveform_dir, waveform_names)

display(sprintf('-------  moment %d ------', p))

for n = 1:numel(waveform_names)
    load(fullfile(waveform_dir,waveform_names{n}))

    % get SPAS from filtered B - p-th moment (not using moments just eigenvectors, so they can be properly normalized later)
    Mp = real(squeeze(sum(wfm.PS .* repmat(wfm.f',1, 3, 3).^p ))); % assuming f >= 0

    [Vm, Lm] = eig(Mp);

    %projections along SPAS
    for k = 1:3 % assuming f >= 0
        outer_prod = Vm(:,k)*Vm(:,k)';
        % spectral projections
        su = real(sum(sum(wfm.PS .* permute(repmat(outer_prod,1,1,length(wfm.f)),[3 1 2]),2),3));
        bu(k) = sum(su);
        fp(k) = sum(wfm.f.^p .* su')/bu(k);
    end

    s = real(wfm.PS(:,1,1) + wfm.PS(:,2,2) + wfm.PS(:,3,3));
    fc = sum(wfm.f.^p .* s')/sum(s);
    
    % Comparison of two centroid frequency calculations for validation 
    %fc_weighted = nansum(fp.*bu)/sum(s);
    %display(sprintf('centroid fc_trace|fc_weighted = %.2f|%.2f Hz', fc, fc_weighted))
    
    display(sprintf('%s', waveform_names{n}))
    display(sprintf('centroid fc_trace = %.2f', fc))
    display(sprintf('eigen f = %.2f | %.2f | %.2f Hz', fp(1), fp(2), fp(3)))
end

end

