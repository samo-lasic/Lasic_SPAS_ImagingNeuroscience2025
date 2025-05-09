function help_data_path_struct2XLSX(data_path_struct, filename)
% make excel table from data_path_struct (adds an 'LPS' field)

for c = 1:numel(data_path_struct)
    data = data_path_struct(c);

    data = rmfield(data,{'root_data_path','roi_names'});
    data.exp_folder_name = extractBefore(data.exp_folder_name,'_');

    data.LPS = 0;
    data.comment = '';

    tmp = sprintf('%d, ', data.select_subfolders);
    data.select_subfolders = join(tmp(1:end-2));

    T = struct2table(data, 'AsArray', true);

    T.Properties.VariableNames = {'exp_folder_name', 'subfolders', 'LPS', 'comment'};
    if c == 1
        writetable(T, filename, 'WriteMode', 'overwrite');
    else
        writetable(T, filename, 'WriteMode', 'append');
    end
end


