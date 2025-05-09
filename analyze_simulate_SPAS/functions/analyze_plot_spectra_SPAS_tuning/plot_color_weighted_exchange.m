
function [hplot, hlegend] = plot_color_weighted_exchange(figureN, legendN, f_in, PS_in, trPS_in, q, bt, dt, b, ODF, parameters);
% Visualize exchange weighting for b-tensors

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

% ------- interpolate to speed up
[f, PS] = interpolate_power_spectra(f,PS, 1000);
tracePS = PS(:,1,1)+PS(:,2,2)+PS(:,3,3);


% ----------- PS projections and normalization --------
u = [ODF.x ODF.y ODF.z];


exchange_Gamma_u = exchange_gamma_tensorial(q,b,dt,u);


showSQ = parameters.showSQ & norm(eye(3)-3*bt/trace(bt))>0.1; % ensure that bt is not spherical


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

Ncolors = size(u,1);
% red to blue
col(:,1) = linspace(1,0,Ncolors);
col(:,2) = zeros(1,Ncolors);
col(:,3) = linspace(0,1,Ncolors);

% map exchange_Gamma_u to range 0-1
exchange_Gamma_min = min(exchange_Gamma_u);
exchange_Gamma_max = max(exchange_Gamma_u);
exchange_Gamma_mean = mean(exchange_Gamma_u);
display(sprintf('exchange: mean = %g, range = [%g-%g] ms', exchange_Gamma_mean*1e3, exchange_Gamma_min*1e3, exchange_Gamma_max*1e3))

exchange_map = (exchange_Gamma_u - exchange_Gamma_min)/(exchange_Gamma_max-exchange_Gamma_min);
c_ind = ceil(1+exchange_map*(Ncolors-1));
C = col(c_ind,:);
exchange_range = linspace(exchange_Gamma_min,exchange_Gamma_max,Ncolors);


% --------------    find spectral PAS
spectral_PAS = find_SPAS(f,PS,cum_tracePS, parameters.SPAS.LF_power_thresh);

hlegend = figure(legendN);
hlegend.Color = 'white';
clf
Nred = 100;
Ngreen = 200;


for n = 1:Ncolors-1
    patch('Faces',[1:4],'Vertices',...
        [exchange_range(n) 0; exchange_range(n+1) 0; exchange_range(n+1) 1; exchange_range(n) 1],...
        'EdgeColor','none','FaceColor',parameters.path_brightness*col(n,:),'FaceAlpha',1);
end


box off
pbaspect([5 1 1])
axis off




hplot = figure(figureN);
clf


pa = patch('Faces',ODF.tri,'Vertices',ODF.verts/max(ODF.Pdir));
axis tight off, axis square, axis equal

make_axes(parameters, spectral_PAS)


if (parameters.show_trajectory) % q trajectory
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

end

