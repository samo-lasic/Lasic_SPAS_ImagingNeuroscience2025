function rotL = rotate_diagonal_tensors_around_XYZ(L,ang)
% rotate diagonal tensors L around XYZ
%size(L) = [m,3]
%size(ang) = [n,3]
%size(L11) = [m,n]

M = size(L,1);
N = size(ang,1);

cAx = repmat(cos(ang(:,1)),1,M)';
sAx = repmat(sin(ang(:,1)),1,M)';

cAy = repmat(cos(ang(:,2)),1,M)';
sAy = repmat(sin(ang(:,2)),1,M)';

cAz = repmat(cos(ang(:,3)),1,M)';
sAz = repmat(sin(ang(:,3)),1,M)';


L11 = repmat(L(:,1),1,N);
L22 = repmat(L(:,2),1,N);
L33 = repmat(L(:,3),1,N);

rotL = zeros(M,N,3,3);

rotL(:,:,1,1) = L22.*(cAx.*sAz - cAz.*sAx.*sAy).^2 + L33.*(sAx.*sAz + cAx.*cAz.*sAy).^2 + L11.*cAy.^2.*cAz.^2
rotL(:,:,1,2) = L11.*cAy.^2.*cAz.*sAz - L33.*(sAx.*sAz + cAx.*cAz.*sAy).*(cAz.*sAx - cAx.*sAy.*sAz) - L22.*(cAx.*cAz + sAx.*sAy.*sAz).*(cAx.*sAz - cAz.*sAx.*sAy);
rotL(:,:,1,3) = L33.*cAx.*cAy.*(sAx.*sAz + cAx.*cAz.*sAy) - L22.*cAy.*sAx.*(cAx.*sAz - cAz.*sAx.*sAy) - L11.*cAy.*cAz.*sAy;
rotL(:,:,2,2) = L22.*(cAx.*cAz + sAx.*sAy.*sAz).^2 + L33.*(cAz.*sAx - cAx.*sAy.*sAz).^2 + L11.*cAy.^2.*sAz.^2;
rotL(:,:,2,3) = L22.*cAy.*sAx.*(cAx.*cAz + sAx.*sAy.*sAz) - L33.*cAx.*cAy.*(cAz.*sAx - cAx.*sAy.*sAz) - L11.*cAy.*sAy.*sAz;
rotL(:,:,3,3) = L11.*sAy.^2 + L33.*cAx.^2.*cAy.^2 + L22.*cAy.^2.*sAx.^2;

rotL(:,:,2,1) = rotL(:,:,1,2);
rotL(:,:,3,1) = rotL(:,:,1,3);
rotL(:,:,3,2) = rotL(:,:,2,3);
