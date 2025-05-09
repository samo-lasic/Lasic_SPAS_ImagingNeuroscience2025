
function Tu = tensor_projection(T, u)
% T(ind1, ind2, ind3):
% - ind1: Consecutive tensor index for projecting time or frequency.
% - ind2, ind3: Indices of the tensor components.

% u(n, i):
% - n: Number of directions.
% - i = 1...3: Components of each direction.


if (numel(size(T)) == 2) % ensure order 3
    t(1,:,:) = T(:,:);
    T = t;
end

for ind1 = 1:size(T,1)
    for n = 1:size(u,1)
        Tu(n,ind1) = 0;
        for ind2 = 1:3
            for ind3 = 1:3
                Tu(n,ind1) = Tu(n,ind1)+T(ind1,ind2,ind3)*u(n,ind2)*u(n,ind3);
            end
        end
    end
end

Tu = squeeze(Tu);
end