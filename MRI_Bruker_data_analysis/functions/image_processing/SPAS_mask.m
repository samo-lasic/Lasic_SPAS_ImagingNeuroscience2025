function out_nii_fn = SPAS_mask(nii_fn,opt)
% ---------- mask (based on contours for each slice) -------------

if nargin < 2
    opt.mask.thresh = 0.1;
    opt.mask.Ncontours = 10;
    opt.density_std = 3;
    opt.kurtosis_pow = .5;
    opt.show_mask = 0;
end

xps = mdm_xps_load(mdm_fn_nii2xps(nii_fn));
[sig, h_sig] = mdm_nii_read(nii_fn);

sig_sum = sum(sig,4);
sz = size(sig_sum);
mask = zeros(sz);



% threshold
x = 1:sz(1);
y = 1:sz(2);
[x1,y1] = ind2sub(sz(1:2),1:prod(sz(1:2)));

for ns = 1:sz(3)
    slice = squeeze(sig_sum(:,:,ns));

    % figure(1),clf
    % plot(sum(slice,2)), return

    cont = fFindContours(x,y,slice,opt.mask.Ncontours);

    CL_inds = find([cont.nCL{:}] > opt.mask.thresh*opt.mask.Ncontours);
    in = [];
    for nCL = 1:length(CL_inds)
        CL_ind = CL_inds(nCL);
        XY = cont.XY{CL_ind};
        in = [in find(inpolygon(x1,y1,XY(1,:),XY(2,:)))];
    end

    %         [X,Y] = meshgrid(x,y);
    %         figure(1),clf
    %         hold on
    %         [M,c] = contour(X,Y,slice',opt.mask.Ncontours);
    %         plot(x1,y1, 'r.');

    mask1 = mask(:,:,ns);
    mask1(in) = 1;

    %         figure(2),clf
    %         hold on
    %         imagesc(slice.*mask1)
    %         plot(x1,y1, 'r.');

    mask(:,:,ns) = mask1;
end

% fill mask holes
mask = mio_mask_fill(mask);

%M = mio_mask_radial(sig_sum);
%s = mdm_s_mask(s, @mio_mask_radial,[],opt);
%s = mdm_s_mask(s, @mio_mask_mic, mask_path, opt);
%s = mdm_s_mask(s, @mio_mask_pca, mask_path, opt);
%s = mdm_s_mask(s, @mio_mask_simple, [], opt);
%s = mdm_s_mask(s, @mio_mask_threshold, [], opt);


% more masking
sig_sum2 = sig_sum .* mask;

for ns = 1:sz(3)
    slice = squeeze(sig_sum2(:,:,ns));

    s_proj1 = sum(slice,1);
    s_proj2 = sum(slice,2);
    k(1) = kurtosis(s_proj1); 
    k(2) = kurtosis(s_proj2); 
    axis_ratio = k.^(-opt.kurtosis_pow);

    %     figure(1),clf
    %     hold on
    %     plot(s_proj1)
    %     plot(s_proj2)
    %     return


    r1 = 0;
    r2 = 0;
    for x = 1:sz(1)
        for y = 1:sz(2)
            r1 = r1 + [x y]*slice(x,y);
            r2 = r2 + [x y].^2*slice(x,y);
        end
    end
    r1 = r1 / sum(slice(:));
    r2 = r2 / sum(slice(:));

    r = opt.density_std*sqrt(r2 - r1.^2);
    %r(1) = r(1)*axis_ratio;
    r = r .* axis_ratio;

    for x = 1:sz(1)
        for y = 1:sz(2)
            dist = ([x y]-r1) ./ r;

            dist2 = sum(dist.^2);
            %mask1(x,y) = dist < 1;
            mask2(x,y, ns) = dist2 < 1;
        end
    end

    %     figure(1),clf
    %     hold on
    %     imagesc(slice .* mask1)
    %     plot(r1(2),r1(1),'ro')

end

mask = mask & mask2;

if opt.show_mask
    sig_sum = tanh(5 * sig_sum / max(sig_sum(:))) .* mask;
    stack = sig_sum * 255;
    stack_sz = size(stack);
    figure
    stack = imresize(stack,stack_sz(1:2) .* h_sig.pixdim(2:3)');
    clf, montage(uint8(stack));
end

out_nii_fn = append_nii_fn(nii_fn, 'mask');
mdm_nii_write(int8(mask),out_nii_fn,h_sig);


% -------------- FUNCTIONS -------------------------------------------------


function cont = fFindContours(x,y,z,NCL)
z = double(z);

minCL = min(z(:));
maxCL = max(z(:));

contLevels = linspace(minCL,maxCL,NCL);
C = contourc(x,y,z',contLevels);

cont.XY = {};
cont.nCL = {};
cont.CL = {};
cont.scaledCL = {};
cont.NCL = NCL;


for nCL = 1:cont.NCL
    CL = contLevels(nCL);
    ind = find(C(1,:) == CL);
    if ~isempty(ind)
        for n = ind
            cont.XY{end+1} = C(:,n+[1:C(2,n)]);
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


