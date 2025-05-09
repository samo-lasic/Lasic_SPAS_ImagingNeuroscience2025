function out_nii_fn = smooth_slices(nii_fn,opt)
% Gaussian smoothing per slice
if (nargin < 2), opt.smooth_sz = 1; end

[I, h] = mdm_nii_read(nii_fn);


for nsl = 1:size(I,3)
    for nvol = 1:size(I,4)
        I1 = double(squeeze(I(:,:,nsl,nvol)));
        %I1 = imfilter(I1, filter);
        Iblur = single(imgaussfilt(I1,opt.smooth_sz));
        
        I(:,:,nsl,nvol) = Iblur;
        
        %                 figure(1),clf
        %                 subplot(2,1,1)
        %                 imagesc(I1)
        %                 subplot(2,1,2)
        %                 imagesc(Iblur)
        %montage({I1,Iblur})
        
    end
end

out_nii_fn = append_nii_fn(nii_fn, 's');
mdm_nii_write(I,out_nii_fn,h);
xps = mdm_xps_load(mdm_fn_nii2xps(nii_fn));
mdm_xps_save(xps, mdm_fn_nii2xps(out_nii_fn));


end
