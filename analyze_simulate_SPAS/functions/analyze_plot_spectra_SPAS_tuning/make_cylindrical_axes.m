
function make_cylindrical_axes(L, r, n_circ, n_seg, r_seg, col)
% L - length of axes
% r - axes width
% n_circ - number of points on the circle
% n_seg - number of segments (divisions of line)
% r_seg - ratio length/spacing of segments
% color (3x3)

if nargin < 6
    col = eye(3);
end

L_spacing = L/n_seg;
L_seg = r_seg*L_spacing;
spacing = linspace(-L+L_seg/2,L-L_seg/2,2*n_seg)+L_seg;

[Xx, Yx, Zx, Xy, Yy, Zy, Xz, Yz, Zz] = makeCylindricalAxes1(L_seg,r,n_circ);



for i = 1:2:length(spacing)
    Xx_tmp = Xx+spacing(i);
    sx(1,i) = fill3(Xx_tmp(1,:),Yx(1,:),Zx(1,:),col(1,:),'LineStyle','none');
    sx(2,i) = fill3(Xx_tmp(2,:),Yx(2,:),Zx(2,:),col(1,:),'LineStyle','none');
    sx(3,i) = surf(Xx_tmp,Yx,Zx,'facecolor',col(1,:),'LineStyle','none');

    Yy_tmp = Xx+spacing(i);
    sy(1,i) = fill3(Xy(1,:),Yy_tmp(1,:),Zy(1,:),col(2,:),'LineStyle','none');
    sy(2,i) = fill3(Xy(2,:),Yy_tmp(2,:),Zy(2,:),col(2,:),'LineStyle','none');
    sy(3,i) = surf(Xy,Yy_tmp,Zy,'facecolor',col(2,:),'LineStyle','none');

    Zz_tmp = Xx+spacing(i);
    sz(1,i) = fill3(Xz(1,:),Yz(1,:),Zz_tmp(1,:),col(3,:),'LineStyle','none');
    sz(2,i) = fill3(Xz(2,:),Yz(2,:),Zz_tmp(2,:),col(3,:),'LineStyle','none');
    sz(3,i) = surf(Xz,Yz,Zz_tmp,'facecolor',col(3,:),'LineStyle','none');
end




    function [Xx, Yx, Zx, Xy, Yy, Zy, Xz, Yz, Zz] = makeCylindricalAxes1(L1,r,n)

        [Xz,Yz,Zz] = cylinder(r,n);

        Zz = 2*L1*(Zz-0.5);

        Rx = [1 0 0; 0 0 -1; 0 1 0];
        Ry = [0 0 1; 0 1 0; -1 0 0];


        XYZ1 = Ry*[Xz(1,:); Yz(1,:); Zz(1,:)];
        XYZ2 = Ry*[Xz(2,:); Yz(2,:); Zz(2,:)];
        Xx(1,:) = XYZ1(1,:); Yx(1,:) = XYZ1(2,:); Zx(1,:) = XYZ1(3,:);
        Xx(2,:) = XYZ2(1,:); Yx(2,:) = XYZ2(2,:); Zx(2,:) = XYZ2(3,:);

        XYZ1 = Rx*[Xz(1,:); Yz(1,:); Zz(1,:)];
        XYZ2 = Rx*[Xz(2,:); Yz(2,:); Zz(2,:)];
        Xy(1,:) = XYZ1(1,:); Yy(1,:) = XYZ1(2,:); Zy(1,:) = XYZ1(3,:);
        Xy(2,:) = XYZ2(1,:); Yy(2,:) = XYZ2(2,:); Zy(2,:) = XYZ2(3,:);

    end
end

