% Copies an interpolated gradient waveform, inserts a gap for the 180° RF pulse, and appends either a polarity-flipped or non-flipped copy after the gap.
% When flip_polarity_after_180 == 1, the resulting waveform is zeroth moment balanced by construction.



clear all
close all
warning ('off','all');

gmr = 2.6751e+08;


gap180 = 5.05e-3;
flip_polarity_after_180 = 1; % 1 or 0 to flip polarity after adding 180 gap

dt_exp = 2e-5; % waveform time resolution 

% Directory where the interpolated waveform files are stored
waveform_dir = '.../waveforms/Bruker7T_NOW_M0_nogap_10_interp'; % example path, adjust as needed

% List of interpolated waveform files to process
wfm_names = {'Bruker7T_NOW_M0_nogap_10_3_interp_g.mat'}; % example path, adjust as needed


Ngap = round(gap180/dt_exp);
display(sprintf('gap time is %.3f ms', Ngap*dt_exp*1000))
gap = zeros(Ngap,3);

for n = 1:numel(wfm_names)
    wfm_name = wfm_names{n};

    load(fullfile(waveform_dir,wfm_name),'g')

    new_name = extractBefore(wfm_name,'.mat');
    if flip_polarity_after_180 == 1
        new_name = [new_name '_flip1.mat'];
        g = [g; gap; -g];
    else
        new_name = [new_name '_flip0.mat'];
        g = [g; gap; g];
    end


    q = gmr*cumsum(g)*dt_exp;
    b = trace(q'*q)*dt_exp;
    maxG = max(abs(g(:)));
    slew = max(max(abs(diff(g))))/dt_exp;

    N180 = round(length(q)/2);
    qc = q(N180,:);
    qc = sqrt(qc*qc'); % magnitude of crushing vector

    disp_str = sprintf('Gmax = %.3f, slew = %g, |qc| = %g, b = %g',maxG,slew,qc,b);
    display(sprintf('%s',disp_str))

    save(fullfile(waveform_dir,new_name),'g')


    % PLOT RESULTS
    t = [0:length(g)-1]*dt_exp;
    m1 = gmr*cumsum(g.*t')*dt_exp;

    figure%,clf
    subplot(3,1,1)
    plot(t,g,'.')
    xlabel('Time [ms]')
    ylabel('Gradient amplitude [mT/m]')
    title(disp_str, 'Interpreter', 'none');

    subplot(3,1,2)
    plot(t,m1,'.')
    ylabel('m1')

    subplot(3,1,3)
    plot(t,q,'.')
    xlabel('Time [ms]')
    ylabel('q')


end

