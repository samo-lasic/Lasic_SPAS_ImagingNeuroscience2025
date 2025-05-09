function [files, names] = SPAS_find_files(all_files, merged_name, map_names_prepend, map_names, map_names_append)
% used in SPAS_make_correlation_maps

if isempty(map_names)
    names = strcat(map_names_prepend, map_names_append);
elseif isempty(map_names_append)
    names = strcat(map_names_prepend, map_names);
else
    names = strcat(map_names_prepend, map_names, '_', map_names_append);
end

names = strcat('_', names,'.nii.gz');

% ensure cell array
if ~iscell(names)
    names = {names};
end

% need loop to keep order of files and prevent alphabetical sorting

for n = 1:numel(names)
    tmp = all_files(contains({all_files.name}, names{n}));
    if numel(tmp) ~= 1
        display(['missing ' names{n}])
        files = [];
        break
    else
        files(n) = tmp;
    end
end

%files = all_files(contains({files.name}, names));
end