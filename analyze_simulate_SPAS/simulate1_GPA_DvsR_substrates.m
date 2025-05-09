% Frequency domain Gaussian phase approximation (GPA) calculations of ADCs for various waveforms and substrates (with varying compartment size and shape). 
% Powder averaging is done by rotating eigen spectra.
% Makes the `DvsR_*.mat` file containing mean diffusivity and diffusion variance as a function of substrate radius for each waveform; 
% Requires `*_info.mat` files in `waveform_dir`; these are combined into `wfms_info.mat` saved in `model_*` folders.


clear all
close all

restoredefaultpath
addpath(genpath(fullfile(pwd,'functions')));

% ====================================
% Set save_wfm_common = 1 the FIRST TIME you run this script
% This will generate wfms_info.mat by collecting spectral traces and thresholded spectra
% After the initial run, you can set save_wfm_common = 0 to reuse the saved wfms_info.mat
% ====================================
save_wfm_common = 1; % Only needed once for a new set of waveforms – collects spectral traces and thresholded spectra


trPS_thresh = 0.999;
save_results = 1;

% path to waveforms 
waveform_dir = fullfile('..', 'waveforms', 'M1_SPAS_2x21ms_rot');

% waveform names
base_waveform_names{1} = 'STE1_21000_5040_20_1051_252_220829_rot';

baseName = 'DvsR';
select_substrate = 1; %1 - cylinder, 2 - spheroid, 3 - sphere, 4 - asymmetric cylinder, 5 - stick

D0 = 2.48e-9; % 28 C (https://dtrx.de/od/diff/index.html#arrhenius)


for n_model = 1:numel(base_waveform_names)
    waveforms_base_name = base_waveform_names{n_model};

    clear waveform_names
    waveform_names{1} = waveforms_base_name;
    waveform_names{end+1} =[waveforms_base_name '_tunedLTE_res'];
    waveform_names{end+1} =[waveforms_base_name '_optTunedLTE_res'];
    waveform_names{end+1} =[waveforms_base_name '_SPAS1'];
    waveform_names{end+1} =[waveforms_base_name '_SPAS2'];
    waveform_names{end+1} =[waveforms_base_name '_SPAS3'];

    model_name = ['model_' waveforms_base_name];
    %model_dir = fullfile(pwd,'model');
    model_dir = fullfile(fileparts(waveform_dir),model_name);
    mkdir(model_dir)


    switch select_substrate
        case 1
            addName = '_cylinder';
            auxP.shape = 2; SaveAsName = [baseName addName];

            R.Rs =  1e-6*linspace(1,20,300);

            Parameters = [R.Rs; D0*ones(1,length(R.Rs))]';

            % for rotating eigen spectra
            load UDSRTriN300 %UDSRTriN15 % UDSRTriN300
            auxP.cP = cos(UDSR.phi);
            auxP.sP = sin(UDSR.phi);
            auxP.cT = cos(UDSR.theta);
            auxP.sT = sin(UDSR.theta);

        case 2
            addName = '_spheroid';
            auxP.shape = 3; SaveAsName = [baseName addName];

            R.Rs = 1e-6*linspace(1,20,100);
            R.axes_ratio = 3;
            [R.R1s, R.R3s] = spheroid_ratio2axes_mat(R.Rs, R.axes_ratio);

            Parameters = [R.R1s(:) R.R3s(:) D0*ones(length(R.R3s(:)),1)];

            % for rotating eigen spectra
            load UDSRTriN300 %UDSRTriN15 % UDSRTriN300
            auxP.cP = cos(UDSR.phi);
            auxP.sP = sin(UDSR.phi);
            auxP.cT = cos(UDSR.theta);
            auxP.sT = sin(UDSR.theta);
        case 3
            addName = '_sphere';
            auxP.shape = 3; SaveAsName = [baseName addName];

            R.Rs = 1e-6*linspace(1,20,100);

            R.axes_ratio = 1;
            [R.R1s, R.R3s] = spheroid_ratio2axes_mat(R.Rs, R.axes_ratio);

            Parameters = [R.R1s(:) R.R3s(:) D0*ones(length(R.R3s(:)),1)];

            % for rotating eigen spectra
            load UDSRTriN15 %UDSRTriN15 % UDSRTriN300
            auxP.cP = cos(UDSR.phi);
            auxP.sP = sin(UDSR.phi);
            auxP.cT = cos(UDSR.theta);
            auxP.sT = sin(UDSR.theta);
        case 4
            addName = '_asymmetric_cylinder';
            auxP.shape = 6; SaveAsName = [baseName addName];

            R.Rs = 1e-6*linspace(1,20,64);

            R.eps = linspace(0,0.99,68);

            [R.R1s, R.R2s] = asym_cylinder2axes_mat(R.Rs, R.eps);

            Parameters = [R.R1s(:) R.R2s(:) D0*ones(length(R.R1s(:)),1)];

            % for rotating eigen spectra
            [phi,theta,psi]=dpfg5design_angles();
            ang = [phi' theta' psi'];
            % add this for compatibility with innerProductEuler
            ang = -ang;
            auxP.cA = cos(ang(:,3));
            auxP.sA = sin(ang(:,3));
            auxP.cB = cos(ang(:,2));
            auxP.sB = sin(ang(:,2));
            auxP.cG = cos(ang(:,1));
            auxP.sG = sin(ang(:,1));

            auxP.weightADCs = 1/length(psi)*ones(length(psi),1);

        case 5
            addName = '_stick';
            auxP.shape = 4; SaveAsName = [baseName addName];

            R.Rs =  1e-6*linspace(1,20,100);

            Parameters = [R.Rs; D0*ones(1,length(R.Rs))]';

            % for rotating eigen spectra
            load UDSRTriN300 %UDSRTriN15 % UDSRTriN300
            auxP.cP = cos(UDSR.phi);
            auxP.sP = sin(UDSR.phi);
            auxP.cT = cos(UDSR.theta);
            auxP.sT = sin(UDSR.theta);

    end


    auxP.tortuosity = 0;

    % collect relevant waveform info
    if save_wfm_common

        for n = 1:numel(waveform_names)
            path = fullfile(waveform_dir,[waveform_names{n} '_info.mat']);

            load(path)

            wfms(n).trPS_thresh = trPS_thresh;

            wfms(n).f0 = wfm.f(max(find(wfm.trPS<trPS_thresh)));
            wfms(n).f_max = max(wfm.f);

            ind = find(abs(wfm.f_full)<wfms(n).f0);
            wfms(n).f_full = wfm.f_full(ind);
            PS_full = wfm.PS_full(ind,:,:);
            tmp_trace = PS_full(:,1,1)+PS_full(:,2,2)+PS_full(:,3,3);
            b = sum(tmp_trace);
            wfms(n).trPS_full = tmp_trace/b; % normalize
            wfms(n).PS_full = PS_full/b;% normalize

            wfms(n).b_delta = wfm.b_delta;
            wfms(n).TE = wfm.TE;
            wfms(n).b = wfm.b;
            wfms(n).name = waveform_names{n};

        end

        save(fullfile(model_dir,'wfms_info.mat'),'wfms')
    else
        load(fullfile(model_dir,'wfms_info.mat'),'wfms')
    end


    R.D0 = D0;
    R.auxP = auxP;

    wfms = wfms(ismember({wfms.name},waveform_names));

    M = numel(wfms);
    N = size(Parameters,1);

    for m = 1:M
        R.b_delta(m) = wfms(m).b_delta;
        R.TE(m) =  wfms(m).TE;
        R.name{m} = wfms(m).name;

        for n = 1:N

            if select_substrate == 4
                ADC = fADCfromWFM_asymmetric(Parameters(n,:),wfms(m),auxP);
                R.MD(m,n) = sum(auxP.weightADCs.*ADC);
                R.V(m,n) = sum(auxP.weightADCs.*ADC.^2) - R.MD(m,n)^2;

            else

                ADC = fADCfromWFM(Parameters(n,:),wfms(m),auxP);
                R.MDpow(m,n) = mean(ADC);
                R.VApow(m,n) = mean(ADC.^2)-mean(ADC)^2;

                % sufficient for MD
                [l1 l2 l3] = fLw(Parameters(n,:),wfms(m),auxP,1);
                trD = (l1+l2+l3)/3;
                R.MD(m,n) = sum(wfms(m).trPS_full.*trD);

                for i = 1:3
                    for j = 1:3
                        Delta(i,j) =  real(sum(wfms(m).PS_full(:,i,j).*(l3-l1)));
                    end
                end

                R.VA(m,n) = 0;
                for i = 1:3
                    for j = 1:3
                        R.VA(m,n) =  R.VA(m,n) + 3*Delta(i,j)^2 - Delta(i,i)*Delta(j,j); % axial symmetry and no compartment size/shape dispersion
                    end
                end

                R.VA(m,n) = 2/45*R.VA(m,n);

            end

            R.ADCs(m,n,:) = ADC;
            disp(sprintf('wfm = %d/%d, R = %d/%d',m,M,n,N))

        end

    end

    if save_results
        save(fullfile(model_dir,[SaveAsName '.mat']),'R')
    end

end




function [a, c] = spheroid_ratio2axes_mat(r1_vec, ratio_vec)
% spheroid axis lengths a and c from radius vector r1_vec and axis ratio vector ratio_vec
c = r1_vec'*ratio_vec;
a = repmat(r1_vec,[length(ratio_vec), 1])';
end


function [phi,theta,psi]=dpfg5design_angles()
% adapted by Samo Lasic to output the Euler angles
% R = R_phi*R_theta*R_psi (ZYZ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by
% Sune Nørhøj Jespersen
% CFIN/MindLab and Dept. of Physics and Astronomy
% Aarhus Universitet
% Nørrebrogade 44, bygn 10G, 5. sal
% 8000 Århus C
% Denmark
% Phone: +45 78463334
% Cell: +45 60896642
% E-mail: sune@cfin.au.dk
% Web: http://www.cfin.au.dk/~sune
%
% This Matlab function computes 60 pairs of diffusion directions for
% rotationally invariant sampling of double PFG MR diffusion experiments,
% as described in ref. 2. This system is based on the theory of
% quadrature on the sphere, more specifically the framework proposed in
% Ref. 1. It uses 60 rotations to perform 60 rotatiosn of the input
% vectors. The input vectors can be arbitrary unit vectors. Free usage is
% permitted, but please acknowledge Refs 1 and 2.
%
% % References:
% 1) Graf, M., Potts, D., 2009. Sampling Sets and Quadrature Formulae on the Rotation Group. Numerical
% Functional Analysis and Optimization 30, 665-688
% 2) Jespersen, S., Lundell, H., Sønderby C. K., Dyrby, T.B. NMR in
% Biomedicine (2013).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Input arguments: 2 unit vectors upon which the the 60 rotations will
%  act. If no input is supplied, these 2 vectors will default to (1,0,0)
%  and (0,0,1).
% Output: 60 pairs of diffusion directions.

% if nargin==0
%     e1=[1,0,0]';
%     e2=[0,0,1]';
% elseif nargin==2
%     e1=e1(:);
%     e2=e2(:);
% else    error('Syntax error. See help');
%     return;
% end

phi=[0 2*pi/5*[1:5] (pi/5+2*pi/5*[1:5]) 0];
theta=[0 ones(1,5)*atan(2) pi-ones(1,5)*atan(2) pi];
psi=zeros(5,12);
for i=1:12
    for k=1:5
        psi(k,i)=pi/5*not(i==1)*not(i==12)+2/5*pi*(k-1)-phi(i);
    end
end
phi=repmat(phi,5,1);
theta=repmat(theta,5,1);
phi=reshape(phi,[1,60]);
theta=reshape(theta,[1,60]);
psi=reshape(psi,[1,60]);

end

function [r1, r2] = asym_cylinder2axes_mat(rV_vec, eps_vec)
% two principal axes r1 and r2 of an asymmetric cylinder from radius vector rV_vec and eccentricity vector eps_vec.
Nr = length(rV_vec);
Neps = length(eps_vec);

for nr = 1:Nr
    rV = rV_vec(nr);
    for neps = 1:Neps
        eps = eps_vec(neps);
        r1(nr,neps) = rV*(1-eps.^2).^(1/4);
        r2(nr,neps) = rV*(1-eps.^2).^(-1/4);
    end
end
end












