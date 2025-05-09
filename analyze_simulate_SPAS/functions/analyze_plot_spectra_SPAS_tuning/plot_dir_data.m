function h = plot_dir_data(figN, Pdir, title_str, ODF, parameters)
% 3D visualization of directional data (`Pdir`) 


fs = 14;
width1 = .25;
height1 = width1;
width2 = .14;
height2 = width2;

lw_ax = 2;
scale_ax = 1.3;

ShowTitle = parameters.ShowTitle;


if isfield(parameters,'maxScale')
    maxL = parameters.maxScale;
else
    maxL = max(abs(Pdir));
end

edgeC = 0.0;
ODFc = Pdir';
if max(ODFc)-min(ODFc)<0.001
    ODFc = mean(ODFc)*ones(size(ODFc));
end

h = figure(figN);
h.Color = 'white';
clf


pos = get(gcf,'Position');

verts = repmat(Pdir',[1 3]).*[sin(ODF.theta).*cos(ODF.phi)...
    sin(ODF.theta).*sin(ODF.phi) ...
    cos(ODF.theta)];

ODFcolor = abs([sin(ODF.theta).*cos(ODF.phi)...
    sin(ODF.theta).*sin(ODF.phi) ...
    cos(ODF.theta)]);

if isfield(parameters,'monochrome')
    if ~isempty(parameters.monochrome)
        tmp = ones(length(ODFcolor),1);
        ODFcolor(:,1) = parameters.monochrome(1)*tmp;
        ODFcolor(:,2) = parameters.monochrome(2)*tmp;
        ODFcolor(:,3) = parameters.monochrome(3)*tmp;
    end
end


p = patch('Faces',ODF.tri,'Vertices',verts);
colormap bone;
caxis([0 max(ODFc)]);

set(p,'FaceColor','interp','FaceVertexCData',ODFcolor,...
    'EdgeColor',edgeC*[1 1 1],'LineWidth',.5,'EdgeAlpha',.5)

axis equal off

axis(scale_ax*maxL*[-1 1 -1 1 -1 1])


line(scale_ax*maxL*1*[-1 1], scale_ax*maxL*0*[-1 1],scale_ax*maxL*0*[-1 1], 'LineStyle','--','linewidth', lw_ax,'color','r');
line(scale_ax*maxL*0*[-1 1], scale_ax*maxL*1*[-1 1],scale_ax*maxL*0*[-1 1], 'LineStyle','--','linewidth', lw_ax,'color','g');
line(scale_ax*maxL*0*[-1 1], scale_ax*maxL*0*[-1 1],scale_ax*maxL*1*[-1 1], 'LineStyle','--','linewidth', lw_ax,'color','b');


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