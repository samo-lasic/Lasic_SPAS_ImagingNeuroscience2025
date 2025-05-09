
function h = plot_SQ(figN, SQ, title_str, parameters)
% Visualize a superquadric (SQ) glyph to represent directional or tensor data with customizable lighting and color mapping.

fs = 14;
width1 = .25;
height1 = width1;
width2 = .14;
height2 = width2;

lw_ax = 2.5;
scale_ax = 1.3;

if ~isfield(parameters,'ambient')
    parameters.ambient = 0.5;
end
if ~isfield(parameters,'specular')
    parameters.specular = 0.5;
end


ShowTitle = parameters.ShowTitle;

Pdir = sqrt(sum(SQ.verts.^2,2))'; %bu(:,1)';

h = figure(figN);
h.Color = 'white';
clf
pos = get(gcf,'Position');

ShowTitle = parameters.ShowTitle;

if isfield(parameters,'maxScale')
    maxL = parameters.maxScale;
else
    maxL = max(abs(Pdir));
end

edgeC = 0.0;
ODFc = Pdir';% Ndir x 1
if max(ODFc)-min(ODFc)<0.001
    ODFc = mean(ODFc)*ones(size(ODFc));
end

p = patch('Faces',SQ.tri,'Vertices',SQ.verts);

colormap bone;
caxis([0 max(ODFc)]);

set(p,'FaceVertexCData',SQ.c)
set(p,'FaceLighting','phong','FaceColor','interp','FaceAlpha',1,...
    'AmbientStrength',parameters.ambient,'SpecularStrength',parameters.specular)

shading interp
axis equal off

axis(scale_ax*maxL*[-1 1 -1 1 -1 1])

% provisional fix
try
    if parameters.showXYZ
        line(scale_ax*maxL*1*[-1 1], scale_ax*maxL*0*[-1 1],scale_ax*maxL*0*[-1 1], 'LineStyle','--','linewidth', lw_ax,'color','r');
        line(scale_ax*maxL*0*[-1 1], scale_ax*maxL*1*[-1 1],scale_ax*maxL*0*[-1 1], 'LineStyle','--','linewidth', lw_ax,'color','g');
        line(scale_ax*maxL*0*[-1 1], scale_ax*maxL*0*[-1 1],scale_ax*maxL*1*[-1 1], 'LineStyle','--','linewidth', lw_ax,'color','b');
    end
catch
    line(scale_ax*maxL*1*[-1 1], scale_ax*maxL*0*[-1 1],scale_ax*maxL*0*[-1 1], 'LineStyle','--','linewidth', lw_ax,'color','r');
    line(scale_ax*maxL*0*[-1 1], scale_ax*maxL*1*[-1 1],scale_ax*maxL*0*[-1 1], 'LineStyle','--','linewidth', lw_ax,'color','g');
    line(scale_ax*maxL*0*[-1 1], scale_ax*maxL*0*[-1 1],scale_ax*maxL*1*[-1 1], 'LineStyle','--','linewidth', lw_ax,'color','b');
end
set(gca,'XTick',[],'YTick',[],'ZTick',[],'FontSize',fs)
view(parameters.az,parameters.el)


if ~isempty(parameters.lightPos1)
    light('Position',parameters.lightPos1,'Style','infinite'), lighting gouraud
end

if ~isempty(parameters.lightPos2)
    light('Position',parameters.lightPos2,'Style','infinite'), lighting gouraud
end


if ShowTitle
    title([title_str ', [' num2str(min(Pdir),3) '-' num2str(max(Pdir),3) '], ' num2str(mean(Pdir),3) ' (' num2str(std(Pdir),3) ')'],'FontSize',fs)
end


end