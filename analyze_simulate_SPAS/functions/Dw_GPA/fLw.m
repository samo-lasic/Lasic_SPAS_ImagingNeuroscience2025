
function [l1 l2 l3] = fLw(Pin,wfmInfo,auxP,full_range)
% can be free or restricted
% wfm has waveform info, auxP etc.
% Pin - size and D0

if full_range == 1
    f = wfmInfo.f_full;
else
    f = wfmInfo.f;
end
Nf = length(f); % auxP.Nf

switch auxP.shape
    case 1
        % sphere
        l1 = DwSpherical(2*pi*f,Pin(1),Pin(2),auxP.tortuosity,50)';
        l2 = l1;
        l3 = l1;

    case 2
        % cylinder
        l1 = DwCylindrical(2*pi*f,Pin(1),Pin(2),auxP.tortuosity,50)';
        l2 = l1;
        l3 = Pin(2)*ones(Nf,1);

    case 3
        % ellipsoid
        l1 = DwSpherical(2*pi*f,Pin(1),Pin(3),auxP.tortuosity,50)';
        l2 = l1;
        l3 = DwSpherical(2*pi*f,Pin(2),Pin(3),auxP.tortuosity,50)';

    case 4
        % stick
        l1 = auxP.tortuosity*Pin(2)*ones(Nf,1); %radial
        l2 = l1; %radial
        l3 = DwPlanar(2*pi*f,2*Pin(1),Pin(2),auxP.tortuosity,50); % axial

    case 5 % Free diffusion
        % Lg = eye(3);
        l3 = Pin(1); % Lg(3,3)
        l2 = Pin(1) + Pin(2); % Lg(2,2)
        l1 = Pin(1) + Pin(2) + Pin(3); % Lg(1,1)

    case 6 % asymmetric cylinder
        l1 = DwCylindrical(2*pi*f,Pin(1),Pin(3),auxP.tortuosity,50)';
        l2 = DwCylindrical(2*pi*f,Pin(2),Pin(3),auxP.tortuosity,50)';
        l3 = Pin(3)*ones(Nf,1);

end


end