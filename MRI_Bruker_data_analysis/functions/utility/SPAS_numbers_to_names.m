function str_array = SPAS_numbers_to_names(merged_folder_numbers, str_prepend);
str_array = {};
for n = 1:length(merged_folder_numbers)
    str_array{end+1} = [str_prepend num2str(merged_folder_numbers(n))];
end
end