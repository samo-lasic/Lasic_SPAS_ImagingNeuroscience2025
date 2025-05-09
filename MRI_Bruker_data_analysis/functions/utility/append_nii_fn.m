function out_nii_fn = append_nii_fn(nii_fn, name)
out_nii_fn = [extractBefore(nii_fn,'.nii') '_' name '.nii.gz'];