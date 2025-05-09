function print_mat(M, str_format)
% print a matrix of values
fprintf(strrep([repmat([str_format '\t'], 1, size(M, 2)) '\n'],'\t\n','\n'), M');