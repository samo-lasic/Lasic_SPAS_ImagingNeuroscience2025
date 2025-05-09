
function ADC = fADCfromWFM_asymmetric(Pin, wfmInfo, auxP)
% calculate ADC for a given waveform, diffusion eigenspectra and their rotations
% eigenspectra can be asymmetric and rotations are given in Euler angles alpha, beta, gamma as R = Ra * Rb * Rg;
% requires cos and sin of angles as auxP.cA, auxP.sA, auxP.cB, auxP.sB, auxP.cG, auxP.sG
% diffusion can be free or restricted
% wfm has waveform info, auxP etc.
% Pin - size and D0


[l1 l2 l3] = fLw(Pin, wfmInfo, auxP, 1);

s = wfmInfo.PS_full; % rank 3


% ---------------- attenuation -----------------------

ADC = real(sum(s(:,1,1).*l1)*(auxP.sA.*auxP.sG - auxP.cA.*auxP.cB.*auxP.cG).^2 + sum(s(:,1,1).*l2)*(auxP.cA.*auxP.sG + auxP.cB.*auxP.cG.*auxP.sA).^2 + sum(s(:,1,1).*l3)*auxP.cG.^2.*auxP.sB.^2 ...
    + sum(s(:,2,2).*l1)*(auxP.cG.*auxP.sA + auxP.cA.*auxP.cB.*auxP.sG).^2 + sum(s(:,2,2).*l2)*(auxP.cA.*auxP.cG - auxP.cB.*auxP.sA.*auxP.sG).^2 + sum(s(:,2,2).*l3)*auxP.sB.^2.*auxP.sG.^2 ...
    + sum(s(:,3,3).*l1)*auxP.cA.^2.*auxP.sB.^2 + sum(s(:,3,3).*l3)*auxP.cB.^2 + sum(s(:,3,3).*l2)*auxP.sA.^2.*auxP.sB.^2 ...
    + 2*(sum(s(:,1,3).*l3)*auxP.cB.*auxP.cG.*auxP.sB + sum(s(:,1,3).*l1)*auxP.cA.*auxP.sB.*(auxP.sA.*auxP.sG - auxP.cA.*auxP.cB.*auxP.cG) + sum(s(:,1,3).*l2)*auxP.sA.*auxP.sB.*(-auxP.cA.*auxP.sG - auxP.cB.*auxP.cG.*auxP.sA)) ...
    + 2*(sum(s(:,1,2).*l3)*auxP.cG.*auxP.sG.*auxP.sB.^2 - sum(s(:,1,2).*l1)*(auxP.cG.*auxP.sA + auxP.cA.*auxP.cB.*auxP.sG).*(auxP.sA.*auxP.sG - auxP.cA.*auxP.cB.*auxP.cG) - sum(s(:,1,2).*l2)*(auxP.cA.*auxP.sG + auxP.cB.*auxP.cG.*auxP.sA).*(auxP.cA.*auxP.cG - auxP.cB.*auxP.sA.*auxP.sG)) ...
    + 2*(sum(s(:,2,3).*l3)*auxP.cB.*auxP.sB.*auxP.sG - sum(s(:,2,3).*l1)*auxP.cA.*auxP.sB.*(auxP.cG.*auxP.sA + auxP.cA.*auxP.cB.*auxP.sG) + sum(s(:,2,3).*l2)*auxP.sA.*auxP.sB.*(auxP.cA.*auxP.cG - auxP.cB.*auxP.sA.*auxP.sG)));


end