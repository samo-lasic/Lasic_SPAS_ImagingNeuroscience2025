
function [hplot, hlegend, parameters] = plot_color_weighted_SPAS_moments(figureN, legendN, f_in, PS_in, trPS_in, ODF, g, q, bt, dt, parameters)
% Patch for plotting spectral moments.
% Contains unused features from old code, to be cleaned in future updates.


if parameters.showTuning & parameters.tuning.thresh > parameters.thresh
    thresh = parameters.tuning.thresh;
    need_rethresh = 1;
else
    thresh = parameters.thresh;
    need_rethresh = 0;
end

% cutoff for color coded spectra (after tuning)
if parameters.f0 == 0
    f0 = max(f_in(trPS_in<parameters.thresh));
else
    f0 = parameters.f0;
end

% limit frq range
ind = find(f_in<=f0);
f = f_in(ind);
PS = PS_in(ind,:,:);

% interpolate to speed up
[f, PS] = interpolate_power_spectra(f,PS, 1000);
tracePS = PS(:,1,1)+PS(:,2,2)+PS(:,3,3);


% ----------- PS projections and normalization --------
u = [ODF.x ODF.y ODF.z];
PSu = real(tensor_projection(PS,u));

% normalize to trace
PSu = PSu/sum(tracePS); % mean(sum(PSu,2)) = 1/3
tracePS = tracePS/sum(tracePS);
% figure(1),clf, hold on, plot(tracePS,'--b')


showSQ = parameters.showSQ & norm(eye(3)-3*bt/trace(bt))>0.1; % ensure that bt is not spherical


% if necessary further threshold spectra
if need_rethresh
    % rough limit frq range
    ind = find(f_in<=f0);
    f = f_in(ind);
    PS = PS_in(ind,:,:);

    % ------- interpolate to speed up
    [f, PS] = interpolate_power_spectra(f,PS, 1000);
    tracePS = PS(:,1,1)+PS(:,2,2)+PS(:,3,3);

    if showSQ == 0
        % ----------- PS projections and normalization --------
        PSu = real(tensor_projection(PS,u));
        % normalize to trace
        PSu = PSu/sum(tracePS); % mean(sum(PSu,2)) = 1/3
    end
    tracePS = tracePS/sum(tracePS);

end


if parameters.showBshape == 1

    if showSQ
        ODF = fSuperQuadTensor(bt,50,3,parameters.SQoffset,[]);

        u = ODF.verts./ODF.Pdir';
        PSu = real(tensor_projection(PS,u));

        % normalize to trace
        PSu = PSu/sum(tracePS); % mean(sum(PSu,2)) = 1/3

        if parameters.showTuning
            Fu = scatteredInterpolant(u(:,1),u(:,2),u(:,3),ODF.Pdir');
            Fu.Method = 'natural'; % 'nearest', 'linear', or 'natural'.
            Fu.ExtrapolationMethod = 'linear'; % 'nearest', 'linear', or 'none'.
        end
    else
        ODF.Pdir = tensor_projection(bt,[ODF.x ODF.y ODF.z]);
        ODF.verts = repmat(ODF.Pdir,[1 3]).*[sin(ODF.theta).*cos(ODF.phi)...
            sin(ODF.theta).*sin(ODF.phi) ...
            cos(ODF.theta)];
    end

else
    showSQ = 0;
    ODF.Pdir = ones(ODF.N,1);
    ODF.verts = repmat(ODF.Pdir,[1 3]).*[sin(ODF.theta).*cos(ODF.phi)...
        sin(ODF.theta).*sin(ODF.phi) ...
        cos(ODF.theta)];
end

% cumulative trace
cum_tracePS = cumsum(tracePS);

% normalize ODF
maxPdir = max(ODF.Pdir);
ODF.verts = ODF.verts/maxPdir;
ODF.Pdir = ODF.Pdir/maxPdir;

% ----------- D(w) for color map ----------------------------
if parameters.colormap.D2Rm4>0
    if parameters.colormap.restricted
        % restricted
        D0 = parameters.colormap.D0;
        R = (D0^2/parameters.colormap.D2Rm4)^(1/4);
        Dw = DwSpherical(2*pi*f,R,D0,0,50)';
    else
        % incoherent flow (%FT of exponential autocorrelation exp(-|t|/tau) -> 2a/(a^2+w^2), a = 1/tau, a = 1/1e-1)
        D0 = parameters.colormap.D0;
        tau = (1/parameters.colormap.D2Rm4)^(1/2);  %tau = 0.001;
        a = 1./tau;
        Dw = 2*a./(a^2+(2*pi*f').^2);
        Dw = D0*Dw/max(Dw);
    end
end


Ncolors = size(PS,1);

if isempty(parameters.binning)
    if parameters.colormap.D2Rm4>0
        % split color according to Dw
        if parameters.colormap.restricted
            Nred = max(find(Dw/D0<1/3));
            Ngreen = max(find(Dw/D0<2/3))-Nred;
            Nblue = Ncolors-Nred-Ngreen;
        else
            Nred = max(find(Dw/D0>2/3));
            Ngreen = max(find(Dw/D0>1/3))-Nred;
            Nblue = Ncolors-Nred-Ngreen;
        end

        parameters.binning(1) = Nred/Ncolors;
        parameters.binning(2) = Ngreen/Ncolors;
        parameters.split_color_weighted_spectra = 'according to even binnin of diffusion spectrum';
    else
        % split color according to power
        Nred = max(find(cum_tracePS<1/3));
        Ngreen = max(find(cum_tracePS<2/3))-Nred;
        Nblue = Ncolors-Nred-Ngreen;
        parameters.binning(1) = Nred/Ncolors;
        parameters.binning(2) = Ngreen/Ncolors;
        parameters.split_color_weighted_spectra = 'according to even binnin of encoding spectral trace';
    end
else
    % split bins
    Nred = round(Ncolors*parameters.binning(1));
    Ngreen = round(Ncolors*parameters.binning(2));
    Nblue = Ncolors-Nred-Ngreen;
    %binning = [];
    parameters.split_color_weighted_spectra = 'according to custom bins';
end
display(['binning = ' num2str(parameters.binning)])
display(['f0 = ' num2str(f0)])

parameters.f0 = f0;


% red-green-blue
col(:,1) = [1*ones(1,Nred) 0*ones(1,Ngreen) 0*ones(1,Nblue)];
col(:,2) = [0*ones(1,Nred) 1*ones(1,Ngreen) 0*ones(1,Nblue)];
col(:,3) = [0*ones(1,Nred) 0*ones(1,Ngreen) 1*ones(1,Nblue)];


% --------------    find spectral PAS
spectral_PAS = find_SPAS(f,PS,cum_tracePS, parameters.SPAS.LF_power_thresh);

hlegend = figure(legendN);
hlegend.Color = 'white';
clf

patch('Faces',[1:4],'Vertices',...
    [0 0; Nred 0; Nred 1; 0 1],...
    'EdgeColor','none','FaceColor',parameters.path_brightness*[1 0 0],'FaceAlpha',1);
patch('Faces',[1:4],'Vertices',...
    [Nred 0; Nred+Ngreen 0; ...
    Nred+Ngreen 1; Nred 1],...
    'EdgeColor','none','FaceColor',parameters.path_brightness*[0 1 0],'FaceAlpha',1);
patch('Faces',[1:4],'Vertices',...
    [Nred+Ngreen 0; Ncolors 0; ...
    Ncolors 1; Nred+Ngreen 1],...
    'EdgeColor','none','FaceColor',parameters.path_brightness*[0 0 1],'FaceAlpha',1);


if parameters.colormap.D2Rm4>0
    plot_function = Dw'/D0; %tracePS/max(tracePS).*Dw'/D0;
else
    plot_function = cum_tracePS;
end


hold on
plot(1:length(plot_function),plot_function,'k:','Linewidth',4)
plot(1:length(tracePS),tracePS/max(tracePS),'k-','Linewidth',4)

if parameters.tuning.colormap_showDw
    D0 = parameters.tuning.contour.D0;
    N = length(parameters.tuning.contour.D2Rm4s);
    for n = 1:N
        D2Rm4 = parameters.tuning.contour.D2Rm4s(n);
        ind = 1+mod(n-1,numel(parameters.tuning.contour.style));
        style = parameters.tuning.contour.style{ind};

        if parameters.tuning.contour.restricted
            R = (D0^2/D2Rm4)^(1/4);
            Dw = DwSpherical(2*pi*f,R,D0,0,50)'/D0; % spherical restriction
        else
            %FT of exponential autocorrelation exp(-|t|/tau) -> 2a/(a^2+w^2), a = 1/tau
            tau = (1/D2Rm4)^(1/2);  %tau = 0.001;
            a = 1./tau;
            Dw = 2*a./(a^2+(2*pi*f').^2);
            Dw = Dw/max(Dw);
        end
        plot(1:length(Dw),Dw,style,'color', 0.99*[1 1 1],'Linewidth',3)
    end
end

box off
pbaspect([5 1 1])
axis off


% ------------------------- COLOR CODE -----------------------

% normalize each direction
PSu_norm = PSu./sum(PSu,2);

if parameters.show_LP_color
    for n = 1:size(PSu,1)
        C(n,:) = sum(repmat(PSu_norm(n,:),3,1)'.*col);
    end
elseif parameters.show_Vw_axis % color code depending on Vw
    for n = 1:size(PSu,1)
        bVw(n) = PSu(n,:) * (2*pi*f').^2; % PSu is normalized to trace
    end
    color_range = [min(bVw) max(bVw)];
    C = color_map_RGB((bVw' - color_range(1)) / diff(color_range));
else % color code depending on centroid frequency
    for n = 1:size(PSu,1)
        fc(n) = PSu_norm(n,:)*f';
    end
    color_range = [min(fc) max(fc)];
    C = color_map_RGB((fc' - color_range(1)) / diff(color_range));
end

hplot = figure(figureN);
clf


pa = patch('Faces',ODF.tri,'Vertices',ODF.verts/max(ODF.Pdir));
axis tight off, axis square, axis equal


% --------------    find spectral PAS
spectral_PAS = find_SPAS(f,PS,cum_tracePS, parameters.SPAS.LF_power_thresh);

display('SPAS eigenvectors')
spectral_PAS.filteredB_shape.V


make_axes(parameters,spectral_PAS);


tmp = parameters;
parameters.SPAS_color = parameters.SPAS_color * .7;
parameters.ax_width = parameters.ax_width * .7;

ms = 16;
ax_length = 1.8;

if parameters.show_Vw_axis
    spectral_m = find_normalized_spectral_moment_tensor(f,PS,cum_tracePS,2);
    make_axes1(parameters,spectral_m, ax_length, ms, [3 2 1]);
else
    fcTensor = find_normalized_spectral_moments_from_projections(f,PSu,u,1);
    make_axes2(parameters,fcTensor, ax_length, ms, [3 2 1]);
end

parameters = tmp;


if (parameters.show_trajectory) % make q trajectory
    hold on
    factor = 1.2;
    tmp = q;
    tmp = parameters.trajectory_scale*tmp/max(abs(tmp(:)));
    plot3(tmp(:,1),tmp(:,2),tmp(:,3),'-k','Linewidth',parameters.trajectory_lw)
end

set(pa,'FaceColor','interp','FaceVertexCData',C,...
    'EdgeColor',parameters.edge_brightness*[1 1 1],'LineWidth',.5,'EdgeAlpha',parameters.edge_alpha,'FaceAlpha',parameters.face_alpha)

set(pa,'FaceLighting','phong',...%'FaceColor','interp',...
    'AmbientStrength',parameters.ambient,'SpecularStrength',parameters.specular)


set(gca,'XTick',[],'YTick',[],'ZTick',[],'FontSize',14)


if ~isempty(parameters.lightPos1)
    light('Position',parameters.lightPos1,'Style','infinite'), lighting gouraud
end

if ~isempty(parameters.lightPos2)
    light('Position',parameters.lightPos2,'Style','infinite'), lighting gouraud
end

h = gca;
h.View = [parameters.az,parameters.el];
h.CameraViewAngle = parameters.CameraViewAngle;

% make colorbar
if ~parameters.show_LP_color
    C = color_map_RGB(linspace(0,1,1000)');
    hlegend = figure(legendN+1);
    hlegend.Color = 'white';
    clf
    for n = 1:length(C)
        patch('Faces',[1:4],'Vertices',...
            [n-1 0; n 0; n 1; n-1 1],...
            'EdgeColor','none','FaceColor',parameters.path_brightness*C(n,:),'FaceAlpha',1);
    end
    box off
    pbaspect([10 1 1])
    axis off

    display(sprintf('legend range = %.2f-%.2f', color_range(1), color_range(2) ))

end

end
