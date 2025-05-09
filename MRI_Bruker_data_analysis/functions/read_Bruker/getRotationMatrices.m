function [DwR, DwGAmpRot, DwGAmp] = getRotationMatrices(scan)
for m = 1:3
    for n = 1:3
        % normalized
        field = sprintf('DwR%d%d',m-1,n-1);
        DwR(:,m,n) = scan.method.(field);

        % gradient-scaled (%) rotation matrix
        % used to obtain scales for each volume
        field = sprintf('DwGAmpRot%d%d',m-1,n-1);
        DwGAmpRot(:,m,n) = scan.method.(field);
    end
end

for n = 1:size(DwGAmpRot,1)
    DwGAmp(n) = norm(squeeze(DwGAmpRot(n,:,:)));
end