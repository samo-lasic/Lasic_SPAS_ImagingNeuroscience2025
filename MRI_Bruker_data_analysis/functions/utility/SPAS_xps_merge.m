function xps = SPAS_xps_merge(nii_fn_cell)
% merge xps inluding wfm source ----------
% xps = mdm_xps_merge(mdm_fn_nii2xps(nii_fn)); % this doesn't work if xps has strings like wfm_src

xps_cell = {};
wfm_src = {};
for n_fn = 1:numel(nii_fn_cell)
    xps_tmp = mdm_xps_load(mdm_fn_nii2xps(nii_fn_cell{n_fn}));
    if isfield(xps_tmp,'wfm_src')
        wfm_src{n_fn} = xps_tmp.wfm_src;
        xps_cell{n_fn} = rmfield(xps_tmp,'wfm_src');
    else
        wfm_src{n_fn} = '';
        xps_cell{n_fn} = xps_tmp;
    end
end
xps = mdm_xps_merge(xps_cell);
xps.wfm_src = wfm_src;


% --- add waveform names
wfm_names_display = [];
for n = 1:numel(nii_fn_cell)
    [~, wfm_name_tmp] = fileparts(extractBefore(nii_fn_cell{n},'.nii.gz'));
    wfm_names{n} = wfm_name_tmp;
    wfm_names_display = [wfm_names_display ', ' wfm_name_tmp ];
end
display(sprintf('%s', wfm_names_display(3:end)))

xps.wfm_names = wfm_names;

