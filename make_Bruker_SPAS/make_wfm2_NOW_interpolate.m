% Interpolates waveforms previously optimized using the Numerical Optimization of gradient Waveforms (NOW, https://github.com/jsjol/NOW).
% Compares interpolated results and saves the waveforms. 

% - Folders named with wfm_names may contain several optimization repetitions.
% - Another output folder (appended _interp) is created for interpolation.
% - Interpolated_wfm.mat has all the waveform details (separate left/right gradient waveforms).
% - g.mat has the effective gradient waveform (accounting for 180 RF pulse).
% - interpolation runs for all optimization repetitions.

clear all
close all
warning ('off','all');

% Specify the root directory containing waveform subdirectories (e.g., waveforms/wfm1, waveforms/wfm2). 
% Each subdirectory should store waveform files to be processed.
waveform_dir = '../waveforms'; % Adjust this path to point to your waveforms directory

% List of waveform subdirectory names to process, each containing .mat files with optimization results
wfm_names = {'Bruker7T_NOW_M1_21',...
    'Bruker7T_NOW_M1_23',...
    'Bruker7T_NOW_M1_25'}; % example names, adjust as needed


% Set the path to the Numerical Optimization of gradient Waveforms (NOW)
NOW_path = '...';

% Add NOW toolbox and its subfolders to MATLAB path
addpath(genpath(NOW_path));


write_interp = 1; % after interpolation

dt_exp = 2e-5; % interpolatation resolution

gmr = 2.6751e+08;

for m = 1:numel(wfm_names)
    d = dir(fullfile(waveform_dir, [wfm_names{m} '/*.mat']));
    clear res
    for n = 1:numel(d)
        load(fullfile(d(n).folder,d(n).name),'problem_result')

        name = extractBefore(d(n).name,'.mat');

        problem = problem_result.problem;
        result = problem_result.result;

        g = result.g;
        dt = result.dt;

        interp_dir = [d(n).folder '_interp'];
        mkdir(interp_dir)

        if (1)
            now_print_requested_and_real_times(problem);
        end

        %% interpolate waveform
        interpolated_wfm = wfm_interp(problem, result, dt_exp);

        if write_interp
            g = interpolated_wfm.g;
            dt = interpolated_wfm.dt;

            save(fullfile(interp_dir, [name '_interp']),'interpolated_wfm')
            save(fullfile(interp_dir, [name '_interp_g']),'g')
        end

        q = gmr*cumsum(g)*dt;
        b = trace(q'*q)*dt;
        maxG = max(abs(g(:)));
        slew = max(max(abs(diff(g))))/dt;

        N180 = round(length(q)/2);
        qc = q(N180,:);
        qc = sqrt(qc*qc'); % magnitude of crushing vector


        display(sprintf('%s : Gmax = %.3f, slew = %g, |qc| = %g, b = %g',name,maxG,slew,qc,b))

        res.maxG(n) = maxG;
        res.slew(n) = slew;
        res.b(n) = b;


        % PLOT RESULTS
        t = [0:problem.dt:problem.totalTimeActual]*1e-3;
        q = gmr*cumsum(result.g)*problem.dt*1e-6; % T/m
        q_interp = gmr*cumsum(interpolated_wfm.g)*interpolated_wfm.dt;
        m1_interp = gmr*cumsum(interpolated_wfm.g.*interpolated_wfm.t')*interpolated_wfm.dt;

        figure
        subplot(3,1,1)
        plot(t,result.g,'o',interpolated_wfm.t,interpolated_wfm.g*1e3,'.')
        xlabel('Time [ms]')
        ylabel('Gradient amplitude [mT/m]')
        title(sprintf('%s: b = %.0f, max slew = %.4g', ...
            name, trace(result.B)*1e-6, max(abs(result.slew(:)))),...
            'Interpreter', 'none');

        subplot(3,1,2)
        plot(interpolated_wfm.t,m1_interp,'.')
        ylabel('m1')

        subplot(3,1,3)
        plot(t,q,'o',interpolated_wfm.t,q_interp,'.')
        xlabel('Time [ms]')
        ylabel('q')


    end

    figure
    hold on
    plot(res.b,'o-')
    plot(res.slew/max(res.slew)*max(res.b),'*-')
    title(wfm_names{m},'interpreter','none')
end






% ------ FUNCTIONS ---------

function exp_wfm = wfm_interp(problem, result, dt_exp)
% dt_exp in seconds

tau = problem.totalTimeActual*1e-3;
tau1 = problem.durationFirstPartActual*1e-3;
tau2 = problem.durationSecondPartActual*1e-3;

dt = problem.dt*1e-3;
exp_wfm.dt = dt_exp;

t = 0:dt:tau;
exp_wfm.N = ceil(tau/dt_exp);
exp_wfm.tau = exp_wfm.N*dt_exp;
exp_wfm.time_stretch = exp_wfm.tau/tau;
t_intep = linspace(0,tau,exp_wfm.N);
exp_wfm.t = t_intep*exp_wfm.time_stretch;

exp_wfm.g = interp1(t,result.g,t_intep)*1e-3; % T/m

exp_wfm.N1 = ceil(tau1/dt_exp)+1;
exp_wfm.tau1 = exp_wfm.N1*dt_exp;

exp_wfm.N2 = ceil(tau2/dt_exp)+1;
exp_wfm.tau2 = exp_wfm.N2*dt_exp;

exp_wfm.g1 = exp_wfm.g(1:exp_wfm.N1,:);
exp_wfm.g2 = exp_wfm.g(end-exp_wfm.N2+1:end,:);
end

