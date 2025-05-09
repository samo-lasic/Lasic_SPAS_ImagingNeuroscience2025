function SPAS_wfms = SPAS_order(xps)
% find wfm indices for 1-SPAS1, 2-SPAS2, 3-SPAS3, 4-STE, 5-tLTE, 6-geoSPAS in this order
% if waveform not present set index to 0
% also give flags to know if waveforms or their combinations are contained in the xps.wfm_names list

contains_SPAS = 0; % contains SPAS1 and SPAS3 (SPAS2 not required)
contains_STE = 0;
contains_geoSPAS = 0;
contains_tLTE = 0;

ordered_wfm_ind = [];
SPAS_wfms.names = {};
for m = 1:3
    wfm_name = sprintf('SPAS%d',m);
    ind = find(contains(xps.wfm_names,wfm_name));
    if length(ind) == 1
        ordered_wfm_ind(end+1) = ind;
        SPAS_wfms.names{end+1} = wfm_name;
    elseif length(ind) == 0
        ordered_wfm_ind(end+1) = 0;
        SPAS_wfms.names{end+1} = '';
    else
        display(sprintf('more than one %s found!',wfm_name))
        return
    end
end
%contains_SPAS = numel(ordered_wfm_ind) > 1;
contains_SPAS = ordered_wfm_ind(1) & ordered_wfm_ind(3);

wfm_name = 'STE';
ind = find(contains(xps.wfm_names,wfm_name));
if length(ind) == 1
    ordered_wfm_ind(end+1) = ind;
    SPAS_wfms.names{end+1} = wfm_name;
    contains_STE = 1;
elseif length(ind) == 0
    ordered_wfm_ind(end+1) = 0;
    SPAS_wfms.names{end+1} = '';
else
    display(sprintf('more than one %s found!',wfm_name))
    return
end

wfm_name = 'tLTE';
ind = find(contains(xps.wfm_names,wfm_name));
if length(ind) == 1
    ordered_wfm_ind(end+1) = ind;
    contains_tLTE = 1;
    SPAS_wfms.names{end+1} = wfm_name;
elseif length(ind) == 0
    ordered_wfm_ind(end+1) = 0;
    SPAS_wfms.names{end+1} = '';
else
    display(sprintf('more than one %s found!',wfm_name))
    return
end

wfm_name = 'geoSPAS';
ind = find(contains(xps.wfm_names,wfm_name));
if length(ind) == 1
    ordered_wfm_ind(end+1) = ind;
    contains_geoSPAS = 1;
    SPAS_wfms.names{end+1} = wfm_name;
elseif length(ind) == 0
    ordered_wfm_ind(end+1) = 0;
    SPAS_wfms.names{end+1} = '';
else
    display(sprintf('more than one %s found!',wfm_name))
    return
end



SPAS_wfms.contains_SPAS = contains_SPAS;
SPAS_wfms.contains_STE = contains_STE;
SPAS_wfms.contains_geoSPAS = contains_geoSPAS;
SPAS_wfms.contains_tLTE = contains_tLTE;
SPAS_wfms.ordered_wfm_ind = ordered_wfm_ind;

