
function mat = map_matrix_to_range_limit(mat, in_min, in_max, min_val, max_val)
mat = (mat-in_min)/(in_max-in_min);
mat = min_val + (max_val - min_val)*mat;

mat(mat>max_val) = max_val;
mat(mat<min_val) = min_val;
end

