
function cont = find_contour_lines(Z, u, parameters)
% Find the contour of Z using the given directions u and parameters.


% sample points
[az,el,~] = cart2sph(u(:,1),u(:,2),u(:,3));

% figure(1),clf, plot3(az,el,Z,'.'), grid on

F = scatteredInterpolant([az; 2*pi+az],[el; el],[Z Z]');
F.Method = 'natural'; % 'nearest', 'linear', or 'natural'.
F.ExtrapolationMethod = 'none'; % 'nearest', 'linear', or 'none'.
% make query grid (shift to avoid discontinuity)
[azq,elq] = ndgrid(pi*linspace(-1,3,2*parameters.Nres),pi/2*linspace(-1,1,parameters.Nres));

% figure(1),clf, plot(azq(:,1),elq(1,:),'.')

tuningq = F(azq(:),elq(:));
tuningq = reshape(tuningq,size(azq,1),size(azq,2));

ind = round(size(azq,2)/2) + [1:size(azq,2)];
azq = azq(ind,:);
elq = elq(ind,:);
tuningq = tuningq(ind,:);

% figure(1),clf
% plot3(pi+az,el,Z,'.',azq,elq,tuningq,'o'), grid on
% title('Linear Interpolation')
% xlabel('x'), ylabel('y'), zlabel('Values')
% legend('Sample data','Interpolated query data','Location','Best')

azq = azq + pi; % just to move the gap out of sight

% ------------- contour levels ---------------------

if parameters.log == 1
    contLevels = logspace(log10(parameters.limits(1)),...
        log10(parameters.limits(2)),...
        parameters.Nlevels);
else
    contLevels = linspace(parameters.limits(1),...
        parameters.limits(2),...
        parameters.Nlevels);
end
if parameters.Nlevels==1
    contLevels = [1 1]*contLevels;
end
% ---------------------------------------------------


C = contourc(elq(1,:)',azq(:,1),tuningq,contLevels);


cont.XYang = {};
cont.nCL = {};
cont.CL = {};
cont.scaledCL = {};% 0-1
cont.NCL = parameters.Nlevels;
minCL = min(contLevels);
maxCL = max(contLevels);

for nCL = 1:cont.NCL
    CL = contLevels(nCL);
    ind = find(C(1,:) == CL);
    if ~isempty(ind)
        for n = ind
            cont.XYang{end+1} = C(:,n+[1:C(2,n)]);
            cont.nCL{end+1} = nCL;
            cont.CL{end+1} = CL;
            if minCL == maxCL
                cont.scaledCL{end+1} = 0;
            else
                cont.scaledCL{end+1} = (CL-minCL)/(maxCL-minCL);
            end
        end
    end
end

for n = 1:numel(cont.XYang)
    az = cont.XYang{n}(2,:); % note the change of order
    el = cont.XYang{n}(1,:);
    [x,y,z] = sph2cart(az+pi,el,1);
    cont.XYZ{n}(:,1) = x;
    cont.XYZ{n}(:,2) = y;
    cont.XYZ{n}(:,3) = z;
end

end

