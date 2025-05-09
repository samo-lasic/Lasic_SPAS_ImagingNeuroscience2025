

function mat = map_matrix_to_range(mat,min_val,max_val)
mat = (mat-min(mat(:)))/(max(mat(:))-min(mat(:)));
mat = min_val + (max_val - min_val)*mat;
end
