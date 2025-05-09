function [hplot, hlegend, out] = plot_tuning(figureN, legendN, f, PS, trPS, ODF, g, q, dt, N0, parameters)
% - Calculate tuning metrics based on power spectra.
% - Visualize tuning using color coding, contour lines, and tuned projections.



u = [ODF.x ODF.y ODF.z];

% limit frq range
ind = find(trPS<parameters.thresh);
f = f(ind);
PS = PS(ind,:,:);

% interpolate to speed up
[f, PS] = interpolate_power_spectra(f,PS, 1000);
tracePS = PS(:,1,1)+PS(:,2,2)+PS(:,3,3);

PSu = real(tensor_projection(PS,u));

% normalize to trace
trace_norm = sum(tracePS);
PSu = PSu/trace_norm; % mean(sum(PSu,2)) = 1/3
tracePS = tracePS/trace_norm;

% tracePSu = mean(PSu);
% figure(1),clf, hold on, plot(3*tracePSu,'-r'), plot(tracePS,'--b')

D0 = parameters.contour.D0;
if parameters.contour.restricted
    R = (D0^2/parameters.contour.D2Rm4)^(1/4);
    Dw = DwSpherical(2*pi*f,R,D0,0,50)'; % spherical restriction
else
    %FT of exponential autocorrelation exp(-|t|/tau) -> 2a/(a^2+w^2), a = 1/tau

    tau = (1/parameters.contour.D2Rm4)^(1/2);
    a = 1./tau;
    Dw = 2*a./(a^2+(2*pi*f').^2);
    Dw = D0*Dw/max(Dw);
end


% figure(1),clf, hold on, plot(tracePS/sum(tracePS)), plot(Dw/sum(Dw))

trD = sum(tracePS.*Dw)/D0;
Du = sum(3*PSu'.*repmat(Dw,1,size(PSu,1)))/D0;

out.tuning.rel_range = [min((Du-trD)/trD) max((Du-trD)/trD)];
display(sprintf('range of (Du-trD)/trD = %g - %g',out.tuning.rel_range(1),out.tuning.rel_range(2)))

out.tuning.v = abs(Du-trD)/trD;
out.tuning.m = mean(out.tuning.v);
out.tuning.std = std(out.tuning.v);
out.tuning.rstd = out.tuning.std/out.tuning.m;

tuning = map_matrix_to_range(out.tuning.v,0,1);

%----  tuning index (just tuning, no other parameters)
[~,tuning_ind] = min(tuning);
tuning_ind = tuning_ind(1);
out.tunedLTE.vec = u(tuning_ind,:);
out.tunedLTE.g = sqrt(3)*g*out.tunedLTE.vec';
out.tunedLTE.q = sqrt(3)*q*out.tunedLTE.vec';
out.tunedLTE.maxG = max(abs(out.tunedLTE.g));
out.tunedLTE.maxSlew = max(abs(diff(out.tunedLTE.g)))/dt;

%-------  tuning indices below the first contour line -------------

inds = find(tuning<=parameters.contour.limits(1));
u = u(inds,:);

out.tuning_cont.amplitude1 = get_magnitude(g,u,dt);
out.tuning_cont.slew_rate1 = get_slew_rate(g,u,dt);

amplitude = map_matrix_to_range(out.tuning_cont.amplitude1.v,0,1);
slew_rate = map_matrix_to_range(out.tuning_cont.slew_rate1.v,0,1);


out.opt = 'amp x slew';
[~,ind] = min(amplitude.*slew_rate);
out.optLTE.vec = u(ind(1),:);

out.optLTE.g = sqrt(3)*g*out.optLTE.vec';
out.optLTE.q = sqrt(3)*q*out.optLTE.vec';

out.optLTE.maxG = max(abs(out.optLTE.g));
out.optLTE.maxSlew = max(abs(diff(out.optLTE.g)))/dt;


% ---------------------------- plot --------------------

bt = q'*q*dt;


% --- find countour lines for tuning
u = [ODF.x ODF.y ODF.z];
out.tuning_cont = find_contour_lines(tuning, u, parameters.contour);

if parameters.plot.showBshape == 1

    showSQ = parameters.plot.showSQ & norm(eye(3)-3*bt/trace(bt))>0.1; % ensure that bt is not spherical
    %parameters.SQoffset = 0;

    if showSQ
        ODF = fSuperQuadTensor(bt,50,3,parameters.plot.SQoffset,[]);
        u = ODF.verts./ODF.Pdir';
        PSu = real(tensor_projection(PS,u));

        % normalize to trace
        PSu = PSu/trace_norm; % mean(sum(PSu,2)) = 1/3

        Du = sum(3*PSu'.*repmat(Dw,1,size(PSu,1)))/D0;
        tuning = map_matrix_to_range(abs(Du-trD)/trD,0,1);
        % radius
        Fu = scatteredInterpolant(u(:,1),u(:,2),u(:,3),ODF.Pdir');
        Fu.Method = 'natural'; % 'nearest', 'linear', or 'natural'.
        Fu.ExtrapolationMethod = 'linear'; % 'nearest', 'linear', or 'none'.

    else
        %u = [ODF.x ODF.y ODF.z];
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

maxPdir = max(ODF.Pdir);
ODF.verts = ODF.verts/maxPdir;
ODF.Pdir = ODF.Pdir/maxPdir;


%color map from min to max Du or use a range relative to Dtrace
% ---------------

Du = Du/trD;
display(sprintf('colormap Du/trD = %g - %g',min(Du),max(Du)))

if isempty(parameters.plot.colormap_limits) % from 0 to 1
    Du = map_matrix_to_range(Du,0,1);
else % limit range
    Du = map_matrix_to_range_limit(Du, parameters.plot.colormap_limits(1), parameters.plot.colormap_limits(2), 0, 1);
end

col = make_color_map(Du,2);

% Shade the tuning band corresponding to the lowest contour threshold for emphasis.
if parameters.contour.do_shade
    %col(inds,:) = repmat([0 1 0],length(inds),1);
    mix = exp_pulse(tuning, 0, parameters.contour.shade); %__/\__
    col = repmat(1-mix',1,3).*col ...
        + repmat(mix',1,3).*repmat(parameters.contour.shade.color,length(col),1);
end


[~,sortIndex] = sort(Du);
col_sorted = col(sortIndex,:); %.*repmat(fdev,3,1)';


hlegend = figure(legendN);
hlegend.Color = 'white';

clf

hold on
for n = 1:length(col)-1
    [px,py] = patch_XY(n:n+1,0*ones(2,1),ones(2,1),1);
    p = patch(px,py,'k');
    p.FaceColor = parameters.plot.patch_brightness*col_sorted(n,:);
    p.FaceAlpha = 1;
    p.EdgeColor = 'none';
    %p.LineStyle = '-';
end


box off
pbaspect([5 1 1])
axis off


%-----


hplot = figure(figureN);
hplot.Color = 'white';
clf

pa = patch('Faces',ODF.tri,'Vertices',ODF.verts);
axis tight off, axis square, axis equal


make_axes(parameters.plot)

if (parameters.plot.show_trajectory) % q trajectory
    hold on
    factor = 1.2;
    tmp = q;
    tmp = parameters.plot.trajectory_scale*tmp/max(abs(tmp(:)));
    plot3(tmp(:,1),tmp(:,2),tmp(:,3),'-k','Linewidth',parameters.plot.trajectory_lw)
end


% ------  plot contour lines for tuning and get final waveforms ------------
L = 1.01;
for n = 1:numel(out.tuning_cont.XYZ)

    XYZ = out.tuning_cont.XYZ{n};

    amplitude = get_magnitude(g,XYZ,dt);
    slew_rate = get_slew_rate(g,XYZ,dt);


    amplitude.v = map_matrix_to_range(amplitude.v,0,1);
    slew_rate.v = map_matrix_to_range(slew_rate.v,0,1);

    c = amplitude.v.*slew_rate.v;
    c1 = amplitude.v;
    c2 = slew_rate.v;

    if parameters.plot.showBshape == 1
        if showSQ
            % interpolate
            XYZ = Fu(XYZ(:,1),XYZ(:,2),XYZ(:,3))/maxPdir.*XYZ;

        else
            XYZ = XYZ.*tensor_projection(bt,XYZ)/maxPdir;
        end
    end

    M = 50;
    strokes = round(length(c)/M);
    M = floor(length(c)/strokes);

    for mStrokes = 1:strokes
        if mStrokes == strokes
            ind = (mStrokes-1)*M+1:length(c);
        else
            ind = (mStrokes-1)*M+[1:M+1];
        end
        m = round(mean(ind));
        plot3(XYZ(ind,1)*L, XYZ(ind,2)*L, XYZ(ind,3)*L,'-','Color',...
            (1-out.tuning_cont.scaledCL{n})*[c1(m),1,c2(m)],'LineWidth',2)
    end
end

L = 1.5;
v1 = out.tunedLTE.vec;
plot3(v1(1)*L*[-1:1],v1(2)*L*[-1:1],v1(3)*L*[-1:1],':','Color',parameters.plot.tunedLTE_color,'Linewidth',2*parameters.plot.ax_width)

if parameters.plot.show_optTuned
    v2 = out.optLTE.vec;
    plot3(v2(1)*L*[-1:1],v2(2)*L*[-1:1],v2(3)*L*[-1:1],':','Color',parameters.plot.optLTE_color,'Linewidth',2*parameters.plot.ax_width)
end

if parameters.plot.showBshape == 1
    if showSQ
        % interpolate (it would be better to calculate projections directly)
        v1 = Fu(v1(:,1),v1(:,2),v1(:,3))/maxPdir.*v1;
        if parameters.plot.show_optTuned
            v2 = Fu(v2(:,1),v2(:,2),v2(:,3))/maxPdir.*v2;
        end
    else
        v1 = v1.*tensor_projection(bt,v1)/maxPdir;
        if parameters.plot.show_optTuned
            v2 = v2.*tensor_projection(bt,v2)/maxPdir;
        end
    end
end

L = 1.02;
plot3(v1(1)*L*[-1 1],v1(2)*L*[-1 1],v1(3)*L*[-1 1],'.','Color',parameters.plot.tunedLTE_color,...
    'MarkerSize',parameters.plot.LTE_MarkerSize)

if parameters.plot.show_optTuned
    plot3(v2(1)*L*[-1 1],v2(2)*L*[-1 1],v2(3)*L*[-1 1],'.','Color',parameters.plot.optLTE_color,...
        'MarkerSize',parameters.plot.LTE_MarkerSize)
end

set(pa,'FaceColor','interp','FaceVertexCData',col,...
    'EdgeColor',parameters.plot.edge_brightness*[1 1 1],'LineWidth',.5,'EdgeAlpha',...
    parameters.plot.edge_alpha,'FaceAlpha',parameters.plot.face_alpha)

set(pa,'FaceLighting','phong',...%'FaceColor','interp',...
    'AmbientStrength',parameters.plot.ambient,'SpecularStrength',parameters.plot.specular)

set(gca,'XTick',[],'YTick',[],'ZTick',[],'FontSize',14)

if ~isempty(parameters.plot.lightPos1)
    light('Position',parameters.plot.lightPos1,'Style','infinite'), lighting gouraud
end

if ~isempty(parameters.plot.lightPos2)
    light('Position',parameters.plot.lightPos2,'Style','infinite'), lighting gouraud
end

h = gca;
h.View = [parameters.plot.az,parameters.plot.el];
h.CameraViewAngle = parameters.plot.CameraViewAngle;


end

