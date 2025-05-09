% This script generates Bruker-ready gradient waveforms.
%
% Inputs:
% - Interpolated SPAS and tuned waveforms stored in the specified `waveform_dir`.
% - A sample Bruker text file (`Bruker_SAMPLE.txt`) used as a template for output formatting.
%
% Outputs:
% - Bruker-ready waveforms including STE, SPAS-LTE, tuned LTE, and additional waveforms.
%   These are scaled, time-aligned, split into pre/post blocks, and formatted as Bruker-compatible text files.
%
% Key Features:
% - Ensures equal b-values across waveforms.
% - Splits waveforms into pre- and post-180-degree pulse segments.
% - Scales SPAS waveforms to a global maximum gradient amplitude (gmax).
% - Exports waveforms as Bruker-compatible text files.

clear all
gmr = 2.6751e+08;

sample_file = fullfile(fileparts(mfilename('fullpath')),'Bruker_SAMPLE.txt');

% Specify the directory containing interpolated SPAS and tuned waveforms.
waveform_dir = '../waveforms/Bruker7T_NOW_M0_nogap_10_interp'; % Adjust this path to point to your waveforms directory

% Output folder for Bruker-ready waveforms
out_folder = 'BrukerReady';

% A common prefix for a set of SPAS waveform files. 
SPAS_base_name = 'Bruker7T_NOW_M0_nogap_10_3_interp_g_flip1';

% Define input and output filenames for additional LTE waveforms to process and save.
additional_LTE_base_name = SPAS_base_name;
additional_LTE_in = strcat(additional_LTE_base_name,{'_optTunedLTE_res','_TunedLTE_res'})
additional_LTE_out = strcat(additional_LTE_base_name,{'_optLTE_res','_tLTE_res'}) % if you want to change names

% Any additional LTE waveforms that should have the same b-value could go here
% additional_LTE_in{end+1} = 'PGSE_te21';
% additional_LTE_out{end+1} = additional_LTE_in{end};

SPAS_ord = [1 2 3]; % e.g. [2 1 3] for [SPAS2 SPAS1 SPAS3]
SPAS_names_in = strcat(SPAS_base_name,{'_SPAS1','_SPAS2','_SPAS3'})


dt = 2e-5;

% in units of fraction of maximum value (not T/m!)
gmax = 1;

for n = 1:3
    wfm_name = strcat(SPAS_names_in{n},'.mat');
    load(fullfile(waveform_dir,wfm_name))
    % ensure the non-zero axis is taken
    gSPAS(:,n) = g(:,find(~prod(g == 0)));

end
gSPAS = gSPAS(:,SPAS_ord);

qSPAS = gmr*cumsum(gSPAS)*dt;
bSPAS = sum(qSPAS.^2)*dt;
bSPAS_max = max(bSPAS);

% plot(qSPAS)

% in case b values are different
br = sqrt(bSPAS_max./bSPAS);
gSPAS = gSPAS.*br;
qSPAS = qSPAS.*br;
bSPAS = sum(qSPAS.^2)*dt;
bSPAS = mean(bSPAS);

[gSPAS_pre, gSPAS_post] = splitg(gSPAS);
for n = 1:length(additional_LTE_in)
    wfm_name = strcat(additional_LTE_in{n},'.mat');
    load(fullfile(waveform_dir,wfm_name))
    % ensure the non-zero axis is taken
    gLTE(:,n) = g(:,find(~prod(g == 0)));
end

% ensure equal b
qLTE = gmr*cumsum(gLTE)*dt;
bLTE = sum(qLTE.^2)*dt;

br = sqrt(bSPAS./bLTE);
gLTE = gLTE.*br;
qLTE = qLTE.*br;

% figure(1),clf
% subplot(2,1,1)
% plot(gLTE)
% subplot(2,1,2)
% plot(qLTE)


[gLTE_pre, gLTE_post] = splitg(gLTE);
%plot(cumsum([gLTE_pre; gLTE_post]))

N_pre = length(gSPAS_pre);
T_pre_ms = round(N_pre*dt*1e3);

N_post = length(gSPAS_post);
T_post_ms = round(N_post*dt*1e3);

% ensure waveform parts have equal length
if length(gLTE_pre) ~= N_pre || length(gLTE_post) ~= N_post
    display('waveform lengths for SPAS and other LTEs are different')
    return
end

% get the real 180 gap
N180gap = length(gSPAS)-N_pre-N_post;
T180gap = N180gap*dt;
display(sprintf('T180gap = %g ms', T180gap*1e3))


% ensure global gmax
gr = gmax/max(abs(gSPAS(:)));
gSPAS = gSPAS*gr;
gSPAS_pre = gSPAS_pre*gr;
gSPAS_post = gSPAS_post*gr;
gLTE = gLTE*gr;
gLTE_pre = gLTE_pre*gr;
gLTE_post = gLTE_post*gr;

qSPAS = gmr*cumsum(gSPAS)*dt;
bSPAS = sum(qSPAS.^2)*dt;
% qLTE = gmr*cumsum(gLTE)*dt;
% bLTE = sum(qLTE.^2)*dt;
bSPAS = mean(bSPAS);
display(['b_max = ' num2str(bSPAS,4)])


% --------------------------- save to txt -----------------------------
out_path = fullfile(waveform_dir,out_folder);

dt_us = round(dt*1e6);
T180gap_us = round(T180gap*1e6);
T_pre_us = round(T_pre_ms*1e3);
T_post_us = round(T_post_ms*1e3);

dateTime = datestr(datetime,'yymmdd');
commonName_pre = makeCommonNameID(T_pre_us, T180gap_us, dt_us, N_pre, N180gap, dateTime);
commonName_post = makeCommonNameID(T_pre_us, T180gap_us, dt_us, N_pre, N180gap, dateTime);

% SPAS
for n = 1:3

    %gradient axis = n; %SPAS1->x, SPAS2->y, SPAS3->z
    vg = [1 2 3] == n;


    %rotation reference axis
    %vr = vg; % change with gradient
    vr = [0 0 1]; % fixed to z

    g_pre = gSPAS_pre(:,n);
    g_post = gSPAS_post(:,n);

    g_pre = repmat(g_pre,1,3);
    g_post = repmat(g_post,1,3);

    g_pre = g_pre.*[0 0 1];
    g_post = g_post.*[0 0 1];


    %reference axis = 1; %SPAS1->x, SPAS2->x, SPAS3->x
    name = sprintf('pre_SPAS%d_%s', n, commonName_pre);
    f_write2bruker(g_pre, sample_file, fullfile(out_path,name),T_pre_ms,vr);

    name = sprintf('post_SPAS%d_%s', n, commonName_post);
    f_write2bruker(-g_post, sample_file, fullfile(out_path,name),T_post_ms,vr);
end

% qiso = gmr*cumsum(gSPAS)*dt;
% b = trace(qiso'*qiso)*dt/3

% STE
axis_order = [1 2 3; 1 3 2]; % flip around x (LF) axis
vr = [0 0 1]; %rotation reference

for n = 1:size(axis_order,1)
    name = sprintf('pre_STE%d_%s', n, commonName_pre);
    f_write2bruker(gSPAS_pre(:,axis_order(n,:))/sqrt(3), sample_file, fullfile(out_path,name), T_pre_ms, vr)

    name = sprintf('post_STE%d_%s', n, commonName_post);
    f_write2bruker(-gSPAS_post(:,axis_order(n,:))/sqrt(3), sample_file, fullfile(out_path,name), T_post_ms, vr)
end

% additional LTEs
vr = [0 0 1]; %rotation reference

for n = 1:numel(additional_LTE_in)
    name = sprintf('pre_%s_%s', additional_LTE_out{n}, commonName_pre);
    f_write2bruker(repmat(gLTE_pre(:,n),1,3).*[0 0 1], sample_file, fullfile(out_path,name), T_pre_ms, vr)

    name = sprintf('post_%s_%s', additional_LTE_out{n}, commonName_post);
    f_write2bruker(-repmat(gLTE_post(:,n),1,3).*[0 0 1], sample_file, fullfile(out_path,name), T_post_ms, vr)
end



% asym3D: Generate asymmetric 3D SPAS gradients by applying weighted contributions to each axis (x, y, z)
if (0)
    %1/2+1/3+1/6 = 1
    vr = [0 0 1]; %rotation reference

    g_pre = gSPAS_pre.*[1/sqrt(2) 1/sqrt(3) 1/sqrt(6)];
    name = sprintf('pre_asym3D_%s', commonName_pre);
    f_write2bruker(g_pre, sample_file, fullfile(out_path,name),T_pre_ms,vr)

    g_post = gSPAS_post.*[1/sqrt(2) 1/sqrt(3) 1/sqrt(6)];
    name = sprintf('post_asym3D_%s', commonName_post);
    f_write2bruker(-g_post, sample_file, fullfile(out_path,name),T_post_ms,vr)
end


% ------ FUNCTIONS -------
function [g1, g2] = splitg(g)
% find 180 gap, g can be 1D, 2D, 3D etc.
SZ = size(g,2);
isZero = g(:,1) == 0;
for n = 2:size(g,2)
    isZero = isZero & g(:,n) == 0;
end
ind180 = find(isZero);
ind180 = ind180(find(diff(ind180) == 1));
g1 = g(1:min(ind180),:);
g2 = g(max(ind180)+1:end,:);
N1 = length(g1);
N2 = length(g2);
if N1 > N2
    g2 = [g2; zeros(N1-N2,SZ)];
elseif N1 < N2
    g1 = [zeros(N2-N1,SZ); g1];
end
end


function f_write2bruker(G,template_file, file_path, dur, dir)
%dir vector defines the direction in G that is rotated. Must be the same as the active axis in LTE

np=length(G);
fid  = fopen(template_file,'r');
f=fread(fid,'*char')';
fclose(fid);

[folder,name] = fileparts(file_path);

f = strrep(f,'PAR_NAME',name);
f = strrep(f,'PAR_NP',num2str(np));
f = strrep(f,'PAR_DUR',num2str(dur));
f = strrep(f,'PAR_DIR',num2str(dir));

mkdir(folder)
fid  = fopen(file_path,'w');
fprintf(fid,'%s',f);
fclose(fid);

for i=1:np
    dlmwrite(file_path, G(i,:),'delimiter',' ','-append','precision','%.6f')
end

dlmwrite(file_path,'##END','delimiter','','-append')
end

function name = makeCommonNameID(T_pre_us, T180gap_us, dt_us, N_pre, N180gap, dateTime)
name = sprintf('%d_%d_%d_%d_%d_%s', T_pre_us, T180gap_us, dt_us, N_pre, N180gap, dateTime);
end