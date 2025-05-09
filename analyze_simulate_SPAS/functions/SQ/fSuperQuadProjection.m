function P = fSuperQuadProjection(T,u,N,gamma,lambda_offset)
% lambda_offset is given in fraction of max(lambda)
% if no RGBval (vector 3x1, e.g. 0.7*[1 1 1]) is given or RGBval = [], color are calculated form vertices
% Original function provided by Prof. Daniel Topgaard, Division of Physical Chemistry, Lund University, Sweden (https://www.physchem.lu.se/people/seniors/topgaard/)
% Adapted by Samo Lasic


T = (T+T')/2; % ensure symmetric matrix, so that eigen-values and -vectors are real

[V, L] = eig(T);

[l,ind] = sort([L(1,1) L(2,2) L(3,3)],'descend');
maxl = max(l);
if maxl == 0
    l = 1e-32*[1 1 1];
else
    l = l + maxl*lambda_offset;
end

R(:,1) = V(:,ind(1));
R(:,2) = V(:,ind(2));
R(:,3) = V(:,ind(3));
%R = eye(3);

% L = [l(1) 0 0; 0 l(2) 0 ; 0 0 l(3)];
% R*L*inv(R)-T

DT.cl = (l(1)-l(2))/(l(1)+l(2)+l(3));
DT.cp = 2*(l(2)-l(3))/(l(1)+l(2)+l(3));
DT.cs = 3*l(3)/(l(1)+l(2)+l(3));

% eval(['load UDSRTriN' num2str(N)])
% SQ = UDSR;
% TR = TriRep(SQ.tri, SQ.x, SQ.y, SQ.z);
% Nsubdiv = 3;
% TR=SubdivideSphericalMesh(TR,Nsubdiv);
% SQ.tri = TR.Triangulation;
% SQ.verts = TR.X;
% SQ.x = TR.X(:,1);
% SQ.y = TR.X(:,2);
% SQ.z = TR.X(:,3);
% SQ.N = numel(SQ.x);
% SQ.theta = acos(SQ.z);
% SQ.phi = atan2(SQ.y,SQ.x);

% [az,el,~] = cart2sph(u(:,1),u(:,2),u(:,3));
% SQ.phi = az;
% SQ.theta = pi/2 - el;

SQ.theta = acos(u(:,3));
SQ.phi = atan2(u(:,2),u(:,1));


SQ.gamma = gamma;

if DT.cl >= DT.cp
    SQ.alpha = (1-DT.cp)^SQ.gamma;
    SQ.beta = (1-DT.cl)^SQ.gamma;
    q.x = l(1)*sign(cos(SQ.theta)).*abs(cos(SQ.theta)).^SQ.beta;
    q.y = l(2)*(-1)*sign(sin(SQ.phi)).*abs(sin(SQ.phi)).^SQ.alpha.*sign(sin(SQ.theta)).*abs(sin(SQ.theta)).^SQ.beta;
    q.z = l(3)*sign(cos(SQ.phi)).*abs(cos(SQ.phi)).^SQ.alpha.*sign(sin(SQ.theta)).*abs(sin(SQ.theta)).^SQ.beta;
elseif DT.cl < DT.cp
    SQ.alpha = (1-DT.cl)^SQ.gamma;
    SQ.beta = (1-DT.cp)^SQ.gamma;
    q.x = l(1)*sign(cos(SQ.phi)).*abs(cos(SQ.phi)).^SQ.alpha.*sign(sin(SQ.theta)).*abs(sin(SQ.theta)).^SQ.beta;
    q.y = l(2)*sign(sin(SQ.phi)).*abs(sin(SQ.phi)).^SQ.alpha.*sign(sin(SQ.theta)).*abs(sin(SQ.theta)).^SQ.beta;
    q.z = l(3)*sign(cos(SQ.theta)).*abs(cos(SQ.theta)).^SQ.beta;
end


SQ_PAS.x = q.x;
SQ_PAS.y = q.y;
SQ_PAS.z = q.z;


SQ_LF.x = R(1,1)*SQ_PAS.x + R(1,2)*SQ_PAS.y + R(1,3)*SQ_PAS.z;
SQ_LF.y = R(2,1)*SQ_PAS.x + R(2,2)*SQ_PAS.y + R(2,3)*SQ_PAS.z;
SQ_LF.z = R(3,1)*SQ_PAS.x + R(3,2)*SQ_PAS.y + R(3,3)*SQ_PAS.z;

SQ.verts = [SQ_LF.x SQ_LF.y SQ_LF.z];

% TR = triangulation(SQ.tri, SQ.verts);
% SQ.norms = vertexNormal(TR,(1:SQ.N)');
% 
% SQ.c = abs(SQ.verts)/max(abs(SQ.verts(:)));
% 
% if nargin==5 && ~isempty(RGBval)
%     %SQ.c = zeros(SQ.N,3);
%     SQ.c(:,1) = RGBval(1)*ones(SQ.N,1);
%     SQ.c(:,2) = RGBval(2)*ones(SQ.N,1);
%     SQ.c(:,3) = RGBval(3)*ones(SQ.N,1);
% end

P = sqrt(sum(SQ.verts.^2,2)); %bu(:,1)';
% used to get a common axis scaling when plotting multiple SQs
% SQ.Pdir = sqrt(sum(SQ.verts.^2,2))'; %bu(:,1)';
% SQ.maxScale = max(abs(SQ.Pdir));
% SQs(n) = SQ;
