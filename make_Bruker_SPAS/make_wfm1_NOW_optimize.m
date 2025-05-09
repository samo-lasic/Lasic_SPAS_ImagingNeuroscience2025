
% Generates multiple repetitions of optimized STE gradient waveforms using Numerical Optimization of gradient Waveforms (NOW, https://github.com/jsjol/NOW).
% Based on the Scripted example of Numerical Optimization of gradient Waveforms (NOW)
% https://github.com/jsjol/NOW

clear all

restoredefaultpath

% Set the path to the Numerical Optimization of gradient Waveforms (NOW)
NOW_path = '...';

% Specify the root directory containing waveform subdirectories (e.g., waveforms/wfm1, waveforms/wfm2). 
% Each subdirectory should store waveform files to be processed.
waveform_dir = '../waveforms'; % Adjust this path to point to your waveforms directory

% Add NOW toolbox and its subfolders to MATLAB path
addpath(genpath(NOW_path));

Nrep = 3; % save Nrep optimization repetitions

write_NOW_txt = 0; % If set to 1, saves optimizer output as .txt files

% Select optimization case for different waveform constraints
% Define the name of the waveform subfolder where optimization results will be saved

Case = 2; % Bruker M0
wfm_name = 'Bruker7T_NOW_M0';

Case = 3; % Bruker M1
wfm_name = 'Bruker7T_NOW_M1';

Case = 4; % Bruker M0 without gap
wfm_name = 'Bruker7T_NOW_M0_nogap';

te1 = 10; % time in ms (left/right waveform parts)
te2 = te1;

wfm_name = sprintf('%s_%d', wfm_name, te1);
mkdir(fullfile(waveform_dir, wfm_name))

switch Case

    case 1 % default

        %%  PREP
        % First, set up the optimization problem. Do this first to create a
        % structure where fields are pre-specified. Note that some fields are
        % read-only and that new fields cannot be created.
        problem = optimizationProblem;

        % Define the hardware specifications of the gradient system
        problem.gMax =  80; % Maximal gradient amplitude, in [mT/m]
        problem.sMax = 100; % Maximal gradient slew (per axis), in [T/(sm)]

        % Request encoding and pause times based on sequence timing in [ms]
        problem.durationFirstPartRequested    = 28; %32;
        problem.durationSecondPartRequested   = 22; %27;
        problem.durationZeroGradientRequested = 8;

        % Define the b-tensor shape in arbitrary units. This example uses an
        % isotropic b-tensor that results in spherical tensor encoding (STE).
        problem.targetTensor = eye(3);

        % Define the number of sample points in time. More points take longer to
        % optimize but provide a smoother waveform that can have steeper slopes.
        problem.N = 50;

        % Set the balance between energy consumption and efficacy
        problem.eta = 0.9; %In interval (0,1]

        % Set the threshold for concomitant gradients (Maxwell terms).
        % Please see https://doi.org/10.1002/mrm.27828 for more information on how
        % to set this parameter.
        problem.MaxwellIndex = 100; %In units of (mT/m)^2 ms

        % Set the desired orders of motion compensation and corresponding
        % thresholds for allowed deviations. See Szczepankiewicz et al., MRM, 2020
        % for details. maxMagnitude in units s^order / m.
        % problem.motionCompensation.order = [1, 2];
        % problem.motionCompensation.maxMagnitude = [0, 1e-4];
        % problem.motionCompensation.order = [1];
        % problem.motionCompensation.maxMagnitude = [0];

    case 2 % Bruker
        %%  PREP
        % First, set up the optimization problem. Do this first to create a
        % structure where fields are pre-specified. Note that some fields are
        % read-only and that new fields cannot be created.
        problem = optimizationProblem;

        % Define the hardware specifications of the gradient system
        problem.gMax =  600; %*sqrt(1/3)*600; % Maximal gradient amplitude, in [mT/m]
        problem.sMax = 3800; %*sqrt(1/3)*1000 % Maximal gradient slew (per axis), in [T/(sm)]

        %enforceSymmetry = true;

        % Request encoding and pause times based on sequence timing in [ms]
        problem.durationFirstPartRequested   = te1; %12-16
        problem.durationSecondPartRequested  = te2; %12-16
        problem.durationZeroGradientRequested = 5.05;

        % Define the b-tensor shape in arbitrary units. This example uses an
        % isotropic b-tensor that results in spherical tensor encoding (STE).
        problem.targetTensor = eye(3);

        % Define the number of sample points in time. More points take longer to
        % optimize but provide a smoother waveform that can have steeper slopes.
        problem.N = 100; %50;

        % Set the balance between energy consumption and efficacy
        problem.eta = 0.9; %In interval (0,1]

        % Set the threshold for concomitant gradients (Maxwell terms).
        % Please see https://doi.org/10.1002/mrm.27828 for more information on how
        % to set this parameter.
        problem.MaxwellIndex = 100; %In units of (mT/m)^2 ms

        % Set the desired orders of motion compensation and corresponding
        % thresholds for allowed deviations. See Szczepankiewicz et al., MRM, 2020
        % for details. maxMagnitude in units s^order / m.
        % problem.motionCompensation.order = [1, 2];
        %problem.motionCompensation.maxMagnitude = [0, 1e-4];
        % problem.motionCompensation.order = [1];
        % problem.motionCompensation.maxMagnitude = [0];

    case 3 % Bruker
        %%  PREP
        % First, set up the optimization problem. Do this first to create a
        % structure where fields are pre-specified. Note that some fields are
        % read-only and that new fields cannot be created.
        problem = optimizationProblem;

        % Define the hardware specifications of the gradient system
        problem.gMax =  600; %*sqrt(1/3)*600; % Maximal gradient amplitude, in [mT/m]
        problem.sMax = 3800; %*sqrt(1/3)*1000 % Maximal gradient slew (per axis), in [T/(sm)]

        %enforceSymmetry = true;

        % Request encoding and pause times based on sequence timing in [ms]
        problem.durationFirstPartRequested   = te1;
        problem.durationSecondPartRequested   = te2;
        problem.durationZeroGradientRequested = 5.05;

        % Define the b-tensor shape in arbitrary units. This example uses an
        % isotropic b-tensor that results in spherical tensor encoding (STE).
        problem.targetTensor = eye(3);

        % Define the number of sample points in time. More points take longer to
        % optimize but provide a smoother waveform that can have steeper slopes.
        problem.N = 100; %50;

        % Set the balance between energy consumption and efficacy
        problem.eta = 0.9; %In interval (0,1]

        % Set the threshold for concomitant gradients (Maxwell terms).
        % Please see https://doi.org/10.1002/mrm.27828 for more information on how
        % to set this parameter.
        problem.MaxwellIndex = 100; %In units of (mT/m)^2 ms

        % Set the desired orders of motion compensation and corresponding
        % thresholds for allowed deviations. See Szczepankiewicz et al., MRM, 2020
        % for details. maxMagnitude in units s^order / m.
        % problem.motionCompensation.order = [1, 2];
        %problem.motionCompensation.maxMagnitude = [0, 1e-4];
        problem.motionCompensation.order = [1];
        problem.motionCompensation.maxMagnitude = [0];

    case 4 % Bruker
        %%  PREP
        % First, set up the optimization problem. Do this first to create a
        % structure where fields are pre-specified. Note that some fields are
        % read-only and that new fields cannot be created.
        problem = optimizationProblem;

        % Define the hardware specifications of the gradient system
        problem.gMax =  600; %*sqrt(1/3)*600; % Maximal gradient amplitude, in [mT/m]
        problem.sMax = 3800; %*sqrt(1/3)*1000 % Maximal gradient slew (per axis), in [T/(sm)]

        %enforceSymmetry = true;

        % Request encoding and pause times based on sequence timing in [ms]
        problem.durationFirstPartRequested   = te1;
        problem.durationSecondPartRequested   = te2;
        problem.durationZeroGradientRequested = 0;

        % Define the b-tensor shape in arbitrary units. This example uses an
        % isotropic b-tensor that results in spherical tensor encoding (STE).
        problem.targetTensor = eye(3);

        % Define the number of sample points in time. More points take longer to
        % optimize but provide a smoother waveform that can have steeper slopes.
        problem.N = 100; %50;

        % Set the balance between energy consumption and efficacy
        problem.eta = 0.9; %In interval (0,1]

        % Set the threshold for concomitant gradients (Maxwell terms).
        % Please see https://doi.org/10.1002/mrm.27828 for more information on how
        % to set this parameter.
        problem.MaxwellIndex = 100; %In units of (mT/m)^2 ms

        % Set the desired orders of motion compensation and corresponding
        % thresholds for allowed deviations. See Szczepankiewicz et al., MRM, 2020
        % for details. maxMagnitude in units s^order / m.
        % problem.motionCompensation.order = [1, 2];
        %problem.motionCompensation.maxMagnitude = [0, 1e-4];
        % problem.motionCompensation.order = [1];
        % problem.motionCompensation.maxMagnitude = [0];

end

% Make a new optimizationProblem object using the updated specifications.
% This explicit call is necessary to update all private variables.
problem = optimizationProblem(problem);

%% PRINT REQUESTED AND TRUE TIMES
% Note that due to the coarse raster, the requested and actual times may
% differ slightly.
clc
now_print_requested_and_real_times(problem);

for nRep = 1:Nrep
  
    %% RUN OPTIMIZATION
    [problem_result.result, problem_result.problem] = NOW_RUN(problem);

    rep_name = sprintf('%s_%d', wfm_name, nRep);
    save(fullfile(waveform_dir, wfm_name, rep_name),'problem_result')
    % txt files from optimizer
    if write_NOW_txt
        fn = now_problem_to_name(problem_result.problem);
        fn = strrep(fn,'NOW',rep_name);
        now_write_wf(problem_result.result, problem_result.problem, fullfile(waveform_dir, wfm_name), fn);
    end

end



