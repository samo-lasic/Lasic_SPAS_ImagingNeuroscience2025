function [struct_without_cells, struct_with_cells] = struct_split(struct_in)
f = fields(struct_in);
for n = 1:length(f)
    if iscell(struct_in.(f{n}))
        struct_with_cells.(f{n}) = struct_in.(f{n});
    else
        struct_without_cells.(f{n}) = struct_in.(f{n});
    end
end
end